"""
emin.py — Single-point calculation of E_min for a Li + Na + Mg (+ Ca) system
=============================================================================
Preliminary implementation (v0.1) of the flowchart:

    Start
      └─ Input: C_i^F, Y, S, CF                      [Excel  → read_inputs]
      └─ I > I_LIMIT?  → "overlimit, change input"   [guard  → check_ionic_strength]
      └─ Input: database (modified Pitzer)           [PHREEQC → PhreeqPython(database=...)]
      └─ Calculate C_i^P, C_i^R                      [Excel  → compute_streams]
      └─ Call PHREEQC, input C^F, C^P, C^R
      └─ Run 1st simulation                          [PHREEQC → build_solutions]
      └─ Any SI > 0? ── Y → equilibrium-phase fn,
      │                     run (n+1)th simulation   [PHREEQC → precipitation_loop]
      │              ── N ↓
      └─ Read a_i^F, a_i^P, a_i^R, m(solid)          [extract_state]
      └─ Yuren eqs → E_min, E_min,i                  [compute_emin]
      └─ Final check: E > 0? ΣE_i = E?               [final_check]
      └─ Print E_min, E_min,i, m(solid)  /  "Calculation Error"

Physics implemented (from the SI, Feng, Abu-Obaid & Yip 2026):

  General:      E_min = (RT / N_Li^P) * [ Σ_i ( n_i^P ln a_i^P
                                              + n_i^R ln a_i^R
                                              - n_i^F ln a_i^F )
                                          - Σ_j n_j^S ln K_sp,j ]         (eq. S-general
                                                                           with solids)
  where i runs over aqueous species incl. H2O, j over precipitated salts
  (halite, bischofite, ...). At precipitation equilibrium SI_j = 0, so
  ln K_sp,j = ln IAP_j is evaluated from the equilibrated product activities
  — no database parsing needed.

  Definitions:  Y  = N_Li^P / N_Li^F           (recovery)
                CF = m_Li^P / m_Li^F           (concentration factor)
                r_i  = m_i^F / m_Li^F          (feed impurity/Li molar ratio)
                S_i  = (m_Li^P/m_i^P)/(m_Li^F/m_i^F)   (Li/i selectivity)
                → m_i^P = (CF / S_i) * m_i^F

Dependencies:  pip install phreeqpython openpyxl pandas
Author: preliminary draft — Claude, for R. Caldentey / Yip Lab pipeline
"""

    from __future__ import annotations

    import math
    import sys
    from dataclasses import dataclass, field

    import phreeqpython

    # ----------------------------------------------------------------------------
    # Constants & configuration
    # ----------------------------------------------------------------------------
    R = 8.314462618e-3        # kJ mol^-1 K^-1
    T = 298.15                # K
    RT = R * T                # kJ mol^-1
    KGW_MOLES = 55.50843      # mol H2O per kg water

    I_LIMIT = 20.0            # mol/kgw guard from the flowchart ("I > 20 M?")
                            # (the extended-Pitzer DB is parameterised to ~40 m)

    SI_TOL = 1e-5             # "any SI > 0?" tolerance (PHREEQC converges to ~1e-8..1e-6)
    MAX_PRECIP_ITER = 10      # safety cap on the n+1 simulation loop
    EBAL_TOL = 1e-6           # |ΣE_i − E| tolerance for the final check

    # Candidate solid phases to monitor (extend when the modified DB adds
    # Li-carnallite, anhydrite, gypsum, ...)
    CANDIDATE_PHASES = {
        #  phase name in DB : {ion: stoichiometric number}, water per formula unit
        "Halite":     {"ions": {"Na+": 1, "Cl-": 1},  "h2o": 0},
        "Bischofite": {"ions": {"Mg+2": 1, "Cl-": 2}, "h2o": 6},
        # "Li-Carnallite": {"ions": {"Li+": 1, "Mg+2": 1, "Cl-": 3}, "h2o": 7},
        # "Anhydrite": {"ions": {"Ca+2": 1, "S(6)": 1}, "h2o": 0},
    }

    # Aqueous master species tracked in the energy balance
    IONS = {"Li+": ("Li", 1), "Na+": ("Na", 1), "Mg+2": ("Mg", 2), "Cl-": ("Cl", -1)}

    MW = {"Li": 6.941, "Na": 22.98977, "Mg": 24.305, "Cl": 35.453}  # g/mol


    # ----------------------------------------------------------------------------
    # 1. Inputs   (flowchart: "Input: C_i^F, Y, S, CF")
    # ----------------------------------------------------------------------------
    @dataclass
    class Inputs:
        feed_mgL: dict[str, float]      # {"Li": ..., "Na": ..., "Mg": ...} in mg/L
        Y: float                        # Li recovery, 0 < Y < 1
        CF: float                       # concentration factor, CF > 1
        S: dict[str, float]             # selectivity S_Li/i for each impurity
        database: str = "pitzer.dat"    # TODO: point to pitzer_mod.dat (Lassin params)
        rho_feed: float = 1.0           # kg solution/L — placeholder mg/L→molal conversion

        def molal_feed(self) -> dict[str, float]:
            """mg/L → mol/kgw.  PRELIMINARY: assumes density ≈ rho_feed and dilute
            solution (kg water ≈ kg solution). For brines, replace with a proper
            density model or run a PHREEQC pre-pass to convert units."""
            return {el: (c / 1000.0 / MW[el]) / self.rho_feed
                    for el, c in self.feed_mgL.items()}


    def read_inputs(xlsx_path: str | None = None) -> Inputs:
        """Read the Excel input sheet (columns: parameter, value). Falls back to a
    hard-coded demo case (Atacama-like brine) when no file is given."""
    if xlsx_path:
        import pandas as pd
        df = pd.read_excel(xlsx_path, sheet_name="inputs", index_col=0)
        v = df["value"]
        return Inputs(
            feed_mgL={"Li": v["Li_F_mgL"], "Na": v["Na_F_mgL"], "Mg": v["Mg_F_mgL"]},
            Y=float(v["Y"]), CF=float(v["CF"]),
            S={"Na": float(v["S_Li_Na"]), "Mg": float(v["S_Li_Mg"])},
            database=str(v.get("database", "pitzer.dat")),
        )
    # demo: Atacama-like feed (Table S2), Y = 0.9, CF = 30, S = 100
    return Inputs(
        feed_mgL={"Li": 1500.0, "Na": 60000.0, "Mg": 10000.0},
        Y=0.9, CF=30.0, S={"Na": 100.0, "Mg": 100.0},
    )


# ----------------------------------------------------------------------------
# 2. Guard   (flowchart: "I > 20 M?")
# ----------------------------------------------------------------------------
def ionic_strength(molal: dict[str, float]) -> float:
    z = {"Li": 1, "Na": 1, "Mg": 2, "Cl": -1}
    return 0.5 * sum(m * z[el] ** 2 for el, m in molal.items())


def check_ionic_strength(streams: dict[str, dict[str, float]]) -> None:
    for name, comp in streams.items():
        I = ionic_strength(comp)
        if I > I_LIMIT:
            sys.exit(f'"overlimit, change input" — stream {name}: '
                     f"I = {I:.2f} mol/kgw > {I_LIMIT}")


# ----------------------------------------------------------------------------
# 3. Stream compositions   (flowchart: "Calculate C_i^P, C_i^R")
# ----------------------------------------------------------------------------
@dataclass
class Streams:
    """Molal compositions (mol/kgw) and water masses (kg) on a 1-kg-feed-water
    basis. n_i = m_i * kgw for each stream."""
    feed: dict[str, float]
    prod: dict[str, float]
    ret: dict[str, float]
    kgw: dict[str, float] = field(default_factory=dict)   # {"F":1, "P":..., "R":...}
    N_Li_P: float = 0.0                                    # mol Li in product (basis)


def compute_streams(inp: Inputs) -> Streams:
    mF = inp.molal_feed()

    # product cations from CF, Y, S:  m_Li^P = CF·m_Li^F ;  m_i^P = (CF/S_i)·m_i^F
    mP = {"Li": inp.CF * mF["Li"]}
    for el in ("Na", "Mg"):
        mP[el] = inp.CF / inp.S[el] * mF[el]

    # Cl- by electroneutrality in each stream
    mF["Cl"] = mF["Li"] + mF["Na"] + 2 * mF["Mg"]
    mP["Cl"] = mP["Li"] + mP["Na"] + 2 * mP["Mg"]

    # water split: basis 1 kg feed water; kgw_P = N_Li^P / m_Li^P = Y/CF
    kgw = {"F": 1.0, "P": inp.Y / inp.CF, "R": 1.0 - inp.Y / inp.CF}
    N_Li_P = inp.Y * mF["Li"]

    # retentate by mass balance:  n_i^R = n_i^F − n_i^P
    mR = {}
    for el in ("Li", "Na", "Mg", "Cl"):
        nR = mF[el] * kgw["F"] - mP[el] * kgw["P"]
        if nR < 0:
            sys.exit(f'"Calculation Error" — negative retentate {el}: '
                     "check Y/CF/S combination (mass balance violated)")
        mR[el] = nR / kgw["R"]
    # NOTE: in the least-work limit (infinitesimal extraction) m^R → m^F and
    # a_i^R → a_i^F (eq. S-12); the exact balance above converges to that
    # limit and stays valid for finite Y.

    return Streams(feed=mF, prod=mP, ret=mR, kgw=kgw, N_Li_P=N_Li_P)


# ----------------------------------------------------------------------------
# 4–6. PHREEQC simulations + precipitation loop
# ----------------------------------------------------------------------------
def build_solutions(pp: phreeqpython.PhreeqPython, st: Streams):
    def make(comp):
        return pp.add_solution({
            "units": "mol/kgw", "temp": T - 273.15, "pH": 7.0,
            "Li": comp["Li"], "Na": comp["Na"], "Mg": comp["Mg"],
            "Cl": f'{comp["Cl"]} charge',      # charge-balance on Cl-
        })
    return make(st.feed), make(st.prod), make(st.ret)     # 1st simulation


def true_si(pp: phreeqpython.Phre0-=00qPython, sol, phase: str) -> float:
    """PHREEQC's SI("phase") via USER_PUNCH.

    GOTCHA: phreeqpython's Solution.si()/phases omit the water-activity term
    for hydrated phases (e.g. Bischofite = MgCl2·6H2O), returning log(IAP')
    without a_w^6. Always use the punched SI for saturation decisions."""
    ip = pp.ip
    ip.accumulate_line(
        f'SELECTED_OUTPUT 1\n    -reset false\n'
        f'USER_PUNCH 1\n    -headings si\n    10 PUNCH SI("{phase}")\n'
        f'RUN_CELLS\n    -cells {sol.number}\nEND')
    ip.dll.RunAccumulated(ip.id_)
    return float(ip.get_selected_output_array()[1][0])


def precipitation_loop(pp, sol_prod):
    """Flowchart right-hand loop: while any SI > 0, call the equilibrium-phase
    function and run the (n+1)th simulation. Precipitation-only (in_phase=0):
    the solids reservoir starts empty, so nothing can dissolve back in."""
    for it in range(1, MAX_PRECIP_ITER + 1):
        oversat = [ph for ph in CANDIDATE_PHASES
                   if true_si(pp, sol_prod, ph) > SI_TOL]
        if not oversat:
            return sol_prod, it
        print(f"  simulation {it + 1}: SI > 0 for {oversat} → EQUILIBRIUM_PHASES")
        sol_prod.equalize(oversat, [0.0] * len(oversat), [0.0] * len(oversat))
    sys.exit('"Calculation Error" — precipitation loop did not converge')


# ----------------------------------------------------------------------------
# 7. Read activities & solid amounts   (flowchart: "Read a_i, m(solid)")
# ----------------------------------------------------------------------------
def water_activity(pp: phreeqpython.PhreeqPython, sol) -> float:
    """phreeqpython's .activity('H2O') is unreliable; punch ACT("H2O")."""
    ip = pp.ip
    ip.accumulate_line(
        f'SELECTED_OUTPUT 1\n    -reset false\n'
        f'USER_PUNCH 1\n    -headings a_w\n    10 PUNCH ACT("H2O")\n'
        f'RUN_CELLS\n    -cells {sol.number}\nEND')
    ip.dll.RunAccumulated(ip.id_)
    return float(ip.get_selected_output_array()[1][0])


@dataclass
class State:
    n: dict[str, dict[str, float]]      # n[stream][species] moles (incl. "H2O")
    a: dict[str, dict[str, float]]      # activities
    solids: dict[str, float]            # n_j^S, mol
    lnKsp: dict[str, float]             # ln K_sp,j (= ln IAP at equilibrium)


def extract_state(pp, st: Streams, solF, solP_eq, solR, solP_initial_moles) -> State:
    n, a = {}, {}
    for tag, sol, kgw in (("F", solF, st.kgw["F"]),
                          ("P", solP_eq, st.kgw["P"]),
                          ("R", solR, st.kgw["R"])):
        # Each PHREEQC cell was built with 1 kg water and scaled to the stream
        # basis by `kgw`. Hydrated solids (bischofite: 6 H2O) consume water on
        # precipitation, so the cell's remaining water mass sol.mass (kg) must
        # multiply both the ion moles and the water moles.
        cell_kgw = sol.mass
        n[tag] = {sp: sol.molality(sp, "mol") * cell_kgw * kgw for sp in IONS}
        n[tag]["H2O"] = KGW_MOLES * cell_kgw * kgw
        a[tag] = {sp: sol.activity(sp, "mol") for sp in IONS}
        a[tag]["H2O"] = water_activity(pp, sol)

    # m(solid) by ion mass balance across the equilibration step, and
    # ln K_sp from equilibrated activities (SI = 0 ⇒ IAP = K_sp)
    solids, lnKsp = {}, {}
    for ph, info in CANDIDATE_PHASES.items():
        # pick a marker ion unique to this phase to quantify precipitation
        marker, nu = next(iter(info["ions"].items()))
        dn = solP_initial_moles[marker] - n["P"][marker]
        n_solid = dn / nu if dn > 1e-12 else 0.0
        if n_solid > 0:
            solids[ph] = n_solid
            lnKsp[ph] = sum(nu_i * math.log(a["P"][ion])
                            for ion, nu_i in info["ions"].items()) \
                        + info["h2o"] * math.log(a["P"]["H2O"])
    # TODO: for phases sharing marker ions (e.g. halite + Li-carnallite both
    # tie up Cl-) switch to reading EQUILIBRIUM_PHASES deltas from PHREEQC
    # directly instead of single-ion mass balance.
    return State(n=n, a=a, solids=solids, lnKsp=lnKsp)


# ----------------------------------------------------------------------------
# 8. Yuren's equations   (flowchart: "Use Yuren eqs to Calculate E_min, E_min,i")
# ----------------------------------------------------------------------------
def compute_emin(st: Streams, s: State) -> tuple[float, dict[str, float]]:
    """E_min and exact species decomposition. Each aqueous species i contributes
        E_i = RT/N * ( n_i^P ln a_i^P + n_i^R ln a_i^R − n_i^F ln a_i^F )
    and each precipitated solid j contributes
        E_j = + RT/N * n_j^S ln K_sp,j .

    SIGN NOTE (⚠ check with Yuren): the SI manuscript writes the solid term
    with a MINUS sign (−n_j ln K_sp). Deriving from μ°_solid = Σν μ°_i +
    RT ln K_sp (K_sp = DISSOLUTION constant, PHREEQC's convention: halite
    log_k = +1.570) gives a PLUS sign, and only the plus sign makes E_min
    continuous across the saturation boundary (moving dn moles from a
    saturated solution into the solid changes E by dn(ln K_sp − ln IAP) = 0
    at SI = 0). The manuscript's minus sign is consistent only if its K is
    the precipitation-reaction constant (1/K_diss)."""
    N = st.N_Li_P

    def nlna(n, a):
        """n ln a with the n→0 limit (n ln a → 0 for absent species)."""
        return n * math.log(a) if n > 1e-15 else 0.0

    E_i: dict[str, float] = {}
    for sp in list(IONS) + ["H2O"]:
        E_i[sp] = RT / N * (nlna(s.n["P"][sp], s.a["P"][sp])
                            + nlna(s.n["R"][sp], s.a["R"][sp])
                            - nlna(s.n["F"][sp], s.a["F"][sp]))
    for ph, n_s in s.solids.items():
        E_i[ph] = RT / N * n_s * s.lnKsp[ph]
    return sum(E_i.values()), E_i


# ----------------------------------------------------------------------------
# 9. Final check & output
# ----------------------------------------------------------------------------
def final_check(E: float, E_i: dict[str, float]) -> None:
    if not (E > 0 and abs(sum(E_i.values()) - E) < EBAL_TOL):
        sys.exit(f'"Calculation Error" — E = {E:.4f} kJ/mol, '
                 f"ΣE_i − E = {sum(E_i.values()) - E:.2e}")


def main(xlsx: str | None = None):
    inp = read_inputs(xlsx)
    st = compute_streams(inp)
    check_ionic_strength({"F": st.feed, "P": st.prod, "R": st.ret})

    print(f"database: {inp.database}")
    pp = phreeqpython.PhreeqPython(database=inp.database)

    solF, solP, solR = build_solutions(pp, st)            # 1st simulation
    print(f"1st simulation done.  I(F) = {solF.I:.2f}, I(P) = {solP.I:.2f} mol/kgw")

    solP_initial = {sp: solP.molality(sp, "mol") * st.kgw["P"] for sp in IONS}
    solP_eq, n_sim = precipitation_loop(pp, solP)

    s = extract_state(pp, st, solF, solP_eq, solR, solP_initial)
    E, E_i = compute_emin(st, s)
    final_check(E, E_i)

    print(f"\nE_min = {E:.3f} kJ (mol Li+)^-1   [{n_sim} product simulation(s)]")
    print("species contributions E_min,i (kJ/mol Li):")
    for sp, v in sorted(E_i.items(), key=lambda kv: -abs(kv[1])):
        print(f"   {sp:<12s} {v:>10.3f}   ({100 * v / E:5.1f} %)")
    if s.solids:
        print("precipitated solids, per mol Li+ recovered:")
        for ph, n_s in s.solids.items():
            print(f"   {ph:<12s} {n_s / st.N_Li_P:.4f} mol/mol Li")
    else:
        print("no oversaturated solids (all SI ≤ 0).")
    return E, E_i, s


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)