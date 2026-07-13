"""
emin_core.py — Core E_min calculation engine, refactored for UI use.

This is the same physics as emin.py (single-point least-work calculation for
Li+Na+Mg brine separation via PHREEQC/Pitzer), restructured so nothing calls
sys.exit() or print()s its way through — every failure mode returns a
CalcResult with success=False and a human-readable message, so a Streamlit
(or any other) front end can display it instead of crashing.

Physics notes (carried over from emin.py / prior validation):
  - Water activity is pulled via a USER_PUNCH ACT("H2O") call, NOT
    phreeqpython's species_activities/activity() helpers, which return a
    constant nonphysical placeholder for H2O in this wrapper version.
  - Saturation indices are pulled the same way (true_si) — phreeqpython's
    .si() omits the water-activity term for hydrated phases.
  - Solutions are charge-balanced on pH, not on Cl-, to avoid spurious
    MgOH+ formation eating into the Mg mass balance at high ionic strength.
  - The general E_min equation (aqueous species + RT/N * n_j * ln K_sp,j
    for any precipitated solids) is evaluated directly from equilibrated
    activities, so no database parsing is needed and K_sp is automatically
    consistent with whatever .dat file is loaded.
"""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import phreeqpython

# ----------------------------------------------------------------------------
R = 8.314462618e-3        # kJ mol^-1 K^-1
T_DEFAULT_C = 25.0
KGW_MOLES = 55.50843      # mol H2O per kg water

SI_TOL = 1e-5
MAX_PRECIP_ITER = 10
EBAL_TOL = 5e-3           # relative tolerance on the ΣE_i = E final check

# Fixed, non-configurable per the source material (not exposed to the UI):
#   - 25 degC: Currently hardcoded value. But CAN allow for variation and provide an advanced setting section to allow users to change.
#   - 20 mol/kgw: the "I > 20 M?" guard from the flowchart, as given.
FIXED_TEMP_C = 25.0
I_LIMIT = 20.0

CANDIDATE_PHASES = {
    "Halite":     {"ions": {"Na+": 1, "Cl-": 1},  "h2o": 0},
    "Bischofite": {"ions": {"Mg+2": 1, "Cl-": 2}, "h2o": 6},
}

IONS = {"Li+": "Li", "Na+": "Na", "Mg+2": "Mg", "Cl-": "Cl"}
MW = {"Li": 6.941, "Na": 22.98977, "Mg": 24.305, "Cl": 35.453}

CONC_UNITS = ("mol/kgw", "mol/L", "ppm")


@dataclass
class CalcInputs:
    Li_val: float
    Na_val: float
    Mg_val: float
    conc_unit: str            # one of CONC_UNITS
    Y: float
    CF: float
    S_Na: float
    S_Mg: float
    database: str = "pitzer.dat"
    database_dir: str | None = None   # None => package default database dir


@dataclass
class CalcResult:
    success: bool
    error: str | None = None
    warnings: list[str] = field(default_factory=list)
    E_min: float | None = None
    E_min_LW: float | None = None     # simplified LiCl-only reference value
    E_components: dict[str, float] = field(default_factory=dict)
    solids: dict[str, float] = field(default_factory=dict)     # mol solid / mol Li recovered
    I_feed: float | None = None
    I_product: float | None = None
    n_simulations: int = 0
    streams: dict[str, dict[str, float]] = field(default_factory=dict)  # molal, for display


def _molal_feed(inp: CalcInputs) -> dict[str, float]:
    """Convert user-entered feed concentrations to molality (mol/kgw) —
    the unit PHREEQC and the E_min equations both use natively.

    Unit handling, exactly as offered in the UI:
      - "mol/kgw" (molality): used as-is. Exact, no assumption involved.
      - "mol/L" (molarity): used numerically as if it were mol/kgw. This
        assumes the feed's density is ~1 kg solution per liter (i.e. 1 L
        of feed contains ~1 kg of water) — the same dilute-solution
        shortcut implicit in Seeley's original mg/L-based code. Exact
        only for dilute feeds; introduces error as the feed gets
        concentrated (their density departs from 1 kg/L).
      - "ppm" (mg solute per kg solution): converted via molar mass to
        mol per kg solution, then used as mol/kgw. At true ppm-level
        concentrations the solute mass is negligible next to the
        solvent's, so kg-solution ~ kg-water is an excellent
        approximation here (much better than for the mol/L case above).
    """
    if inp.conc_unit == "ppm":
        return {
            "Li": inp.Li_val / 1000.0 / MW["Li"],
            "Na": inp.Na_val / 1000.0 / MW["Na"],
            "Mg": inp.Mg_val / 1000.0 / MW["Mg"],
        }
    elif inp.conc_unit in ("mol/kgw", "mol/L"):
        return {"Li": inp.Li_val, "Na": inp.Na_val, "Mg": inp.Mg_val}
    else:
        raise ValueError(f"Unknown concentration unit: {inp.conc_unit!r}")


def _ionic_strength(molal: dict[str, float]) -> float:
    z = {"Li": 1, "Na": 1, "Mg": 2, "Cl": -1}
    return 0.5 * sum(m * z[el] ** 2 for el, m in molal.items())


def _compute_streams(inp: CalcInputs):
    mF = _molal_feed(inp)
    mP = {"Li": inp.CF * mF["Li"]}
    mP["Na"] = inp.CF / inp.S_Na * mF["Na"] if inp.S_Na > 0 else 0.0
    mP["Mg"] = inp.CF / inp.S_Mg * mF["Mg"] if inp.S_Mg > 0 else 0.0
    mF["Cl"] = mF["Li"] + mF["Na"] + 2 * mF["Mg"]
    mP["Cl"] = mP["Li"] + mP["Na"] + 2 * mP["Mg"]

    kgw = {"F": 1.0, "P": inp.Y / inp.CF, "R": 1.0 - inp.Y / inp.CF}
    N_Li_P = inp.Y * mF["Li"]

    mR = {}
    for el in ("Li", "Na", "Mg", "Cl"):
        nR = mF[el] * kgw["F"] - mP[el] * kgw["P"]
        if nR < -1e-12:
            return None, None, None, None, (
                f"Negative retentate {el} (mass balance violated) — "
                f"lower Y, raise CF less aggressively, or raise S_{el}."
            )
        mR[el] = max(nR, 0.0) / kgw["R"]
    return mF, mP, mR, kgw, None


def _true_si(pp, sol, phase: str) -> float:
    ip = pp.ip
    ip.accumulate_line(
        f'SELECTED_OUTPUT 1\n    -reset false\n'
        f'USER_PUNCH 1\n    -headings si\n    10 PUNCH SI("{phase}")\n'
        f'RUN_CELLS\n    -cells {sol.number}\nEND')
    ip.dll.RunAccumulated(ip.id_)
    return float(ip.get_selected_output_array()[1][0])


def _water_activity(pp, sol) -> float:
    ip = pp.ip
    ip.accumulate_line(
        f'SELECTED_OUTPUT 1\n    -reset false\n'
        f'USER_PUNCH 1\n    -headings a_w\n    10 PUNCH ACT("H2O")\n'
        f'RUN_CELLS\n    -cells {sol.number}\nEND')
    ip.dll.RunAccumulated(ip.id_)
    return float(ip.get_selected_output_array()[1][0])


def run_calculation(inp: CalcInputs) -> CalcResult:
    warnings: list[str] = []

    if inp.Li_val <= 0:
        return CalcResult(success=False, error="Feed Li\u207a concentration must be greater than 0.")

    mF, mP, mR, kgw, err = _compute_streams(inp)
    if err:
        return CalcResult(success=False, error=err)

    I_F, I_P = _ionic_strength(mF), _ionic_strength(mP)
    if I_F > I_LIMIT or I_P > I_LIMIT:
        return CalcResult(
            success=False,
            error=(f'Overlimit — feed I = {I_F:.2f}, product I = {I_P:.2f} '
                   f'mol/kgw exceeds the guard of {I_LIMIT} mol/kgw from the flowchart.'),
        )

    try:
        if inp.database_dir:
            pp = phreeqpython.PhreeqPython(
                database=inp.database, database_directory=Path(inp.database_dir))
        else:
            pp = phreeqpython.PhreeqPython(database=inp.database)
    except Exception as e:
        return CalcResult(success=False, error=f"Could not load database '{inp.database}': {e}")

    def make(comp):
        return pp.add_solution({
            "units": "mol/kgw", "temp": FIXED_TEMP_C, "pH": "7 charge",
            "Li": comp["Li"], "Na": comp["Na"], "Mg": comp["Mg"], "Cl": comp["Cl"],
        })

    try:
        solF, solP, solR = make(mF), make(mP), make(mR)
    except Exception as e:
        return CalcResult(success=False, error=f"PHREEQC could not equilibrate the input streams: {e}")

    solP_initial = {sp: solP.molality(sp, "mol") * kgw["P"] for sp in IONS}

    n_sim = 1
    for it in range(1, MAX_PRECIP_ITER + 1):
        try:
            oversat = [ph for ph in CANDIDATE_PHASES if _true_si(pp, solP, ph) > SI_TOL]
        except Exception as e:
            return CalcResult(success=False, error=f"SI lookup failed: {e}")
        if not oversat:
            break
        warnings.append(f"Simulation {it + 1}: {', '.join(oversat)} oversaturated — precipitating.")
        solP.equalize(oversat, [0.0] * len(oversat), [0.0] * len(oversat))
        n_sim += 1
    else:
        return CalcResult(success=False, error="Precipitation loop did not converge within "
                                                 f"{MAX_PRECIP_ITER} iterations.")

    n, a = {}, {}
    try:
        for tag, sol, kg in (("F", solF, kgw["F"]), ("P", solP, kgw["P"]), ("R", solR, kgw["R"])):
            cell_kgw = sol.mass
            n[tag] = {sp: sol.molality(sp, "mol") * cell_kgw * kg for sp in IONS}
            n[tag]["H2O"] = KGW_MOLES * cell_kgw * kg
            a[tag] = {sp: sol.activity(sp, "mol") for sp in IONS}
            a[tag]["H2O"] = _water_activity(pp, sol)
    except Exception as e:
        return CalcResult(success=False, error=f"Could not read back activities: {e}")

    solids, lnKsp = {}, {}
    for ph, info in CANDIDATE_PHASES.items():
        marker, nu = next(iter(info["ions"].items()))
        dn = solP_initial[marker] - n["P"][marker]
        n_solid = dn / nu if dn > 1e-12 else 0.0
        if n_solid > 0:
            solids[ph] = n_solid
            lnKsp[ph] = sum(nu_i * math.log(a["P"][ion]) for ion, nu_i in info["ions"].items()) \
                        + info["h2o"] * math.log(a["P"]["H2O"])

    N = n["P"]["Li+"] if "Li+" in n["P"] else mP["Li"] * kgw["P"]
    N = mP["Li"] * kgw["P"]  # = Y * mF['Li'], the recovered-Li basis

    def nlna(nn, aa):
        return nn * math.log(aa) if nn > 1e-15 else 0.0

    RT = R * T_from_C(FIXED_TEMP_C)
    E_i: dict[str, float] = {}
    for sp in list(IONS) + ["H2O"]:
        E_i[sp] = RT / N * (
            nlna(n["P"][sp], a["P"][sp]) + nlna(n["R"][sp], a["R"][sp]) - nlna(n["F"][sp], a["F"][sp])
        )
    for ph, n_s in solids.items():
        E_i[ph] = RT / N * n_s * lnKsp[ph]

    E = sum(E_i.values())

    if E <= 0:
        warnings.append(f"E_min came out non-positive ({E:.3f} kJ/mol) — check inputs; "
                         "this usually means Y, CF, S are inconsistent with each other.")
    if abs(sum(E_i.values()) - E) > EBAL_TOL * max(abs(E), 1.0):
        warnings.append("Energy balance check (ΣE_i = E) is off by more than tolerance — "
                         "treat this result with caution.")

    display_species = {"Li+": "Li\u207a", "Na+": "Na\u207a", "Mg+2": "Mg\u00b2\u207a",
                        "Cl-": "Cl\u207b", "H2O": "H\u2082O"}
    E_components = {display_species.get(k, k): v for k, v in E_i.items()}

    #Below is the simplified E_min, LW taken from Seeley's "Theoretical
    # Limit" formulation. It basically treats the brine as pure LiCl and ignores
    #the Na/Mg. Wont equal the E_min exactly.

    E_min_LW = RT * math.log((a["P"]["Li+"] / a["F"]["Li+"]) * (a["P"]["Cl-"] / a["F"]["Cl-"]))

    # Water is added to each stream dict two ways: "H2O" as its molality
    # (always ~55.51 mol/kgw by definition — included for completeness,
    # but it's a constant, not informative on its own) and "H2O_kgw" as
    # the water MASS relative to the feed's (1.0, Y/CF, 1-Y/CF) — this is
    # the quantity that actually varies across streams and is genuinely
    # informative about how concentrated each stream is.
    for stream_dict, tag in ((mF, "F"), (mP, "P"), (mR, "R")):
        stream_dict["H2O"] = KGW_MOLES
        stream_dict["H2O_kgw"] = kgw[tag]

    return CalcResult(
        success=True,
        warnings=warnings,
        E_min=E,
        E_min_LW=E_min_LW,
        E_components=E_components,
        solids={ph: v / N for ph, v in solids.items()},
        I_feed=I_F,
        I_product=I_P,
        n_simulations=n_sim,
        streams={"Feed (mol/kgw)": mF, "Product (mol/kgw)": mP, "Retentate (mol/kgw)": mR},
    )


def T_from_C(temp_C: float) -> float:
    return temp_C + 273.15