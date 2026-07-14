"""
app.py — Streamlit front end for the E_min (minimum separation energy) calculator.

Run locally with:
    pip install streamlit phreeqpython plotly
    streamlit run app.py

This wraps emin_core.py (must be in the same folder) — a PHREEQC/Pitzer-based
implementation of Yuren's single-point E_min calculation for a Li+Na+Mg brine,
including the oversaturation / precipitation loop from the flowchart.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st

from emin_core import CalcInputs, run_calculation, MW

APP_EXPECTS_SCHEMA = 3
try:
    from emin_core import SCHEMA_VERSION
except ImportError:
    SCHEMA_VERSION = 0  # emin_core.py predates the version marker entirely

st.set_page_config(page_title="E_min Calculator", layout="wide")

if SCHEMA_VERSION != APP_EXPECTS_SCHEMA:
    st.error(
        f"**emin_core.py is out of date** (schema {SCHEMA_VERSION}, this app.py needs "
        f"schema {APP_EXPECTS_SCHEMA}). The file next to app.py wasn't replaced with the "
        "latest version — re-download emin_core.py and overwrite the old one, then rerun. "
        "(Run `python3 -c \"import emin_core; print(emin_core.__file__)\"` in your terminal "
        "if you're unsure which file is actually being loaded.)"
    )
    st.stop()

st.markdown("""
<style>

/* Main app background */
[data-testid="stAppViewContainer"]{
    background-color:#1e5473;
    color:white;
}

/* Sidebar */
[data-testid="stSidebar"]{
    background-color:#749ef2;
    color:white;
}

/* Main text */
[data-testid="stMarkdownContainer"],
p,
label,
span,
div,
h1,
h2,
h3,
h4,
h5,
h6{
    color:white !important;
}

/* Metric values */
[data-testid="stMetricValue"]{
    color:white !important;
}

/* Metric labels */
[data-testid="stMetricLabel"]{
    color:white !important;
}

/* Tabs */
button[data-baseweb="tab"]{
    color:white !important;
}

</style>
""", unsafe_allow_html=True)


# ----------------------------------------------------------------------------
# Sidebar — inputs
# ----------------------------------------------------------------------------
logo_path = Path(__file__).parent / "YIP_LAB_Logo.png"
if logo_path.exists():
    st.sidebar.image(str(logo_path), width="stretch")
else:
    st.sidebar.caption("(YIP_LAB_Logo.png not found next to app.py — place it there to show the logo.)")

st.sidebar.title("Inputs")

st.sidebar.subheader("Feed composition")

# --- Unit-conversion plumbing ------------------------------------------
# The person's entered concentrations are stored ONCE, canonically, in
# mol/kgw (in st.session_state) regardless of which unit is on screen.
# Switching units re-derives the displayed number from that canonical
# value (real conversion); editing a number re-derives the canonical
# value from what was typed (also a real conversion) — so no data is
# silently reinterpreted under a different unit.
UNIT_LABELS = {"mol/kgw (molality)": "mol/kgw", "mol/L (molarity)": "mol/L", "ppm (mg/kg)": "ppm"}
IONS_UI = ("Li", "Na", "Mg")

st.session_state.setdefault("canon_Li", 0.01)   # mol/kgw
st.session_state.setdefault("canon_Na", 0.0)
st.session_state.setdefault("canon_Mg", 0.0)
st.session_state.setdefault("conc_unit_select", "mol/kgw (molality)")


def _to_molkgw(value: float, unit: str, ion: str) -> float:
    if unit == "ppm":
        return value / 1000.0 / MW[ion]
    return value  # mol/kgw and mol/L both used numerically as-is (see help text)


def _from_molkgw(value: float, unit: str, ion: str) -> float:
    if unit == "ppm":
        return value * 1000.0 * MW[ion]
    return value


def _sync_display_from_canon():
    """Selectbox on_change: unit just changed — recompute the three
    displayed numbers from the canonical mol/kgw values so they show the
    equivalent concentration in the newly selected unit."""
    unit = UNIT_LABELS[st.session_state["conc_unit_select"]]
    for ion in IONS_UI:
        st.session_state[f"{ion}_display"] = round(
            _from_molkgw(st.session_state[f"canon_{ion}"], unit, ion), 6)


def _sync_canon_from_display(ion: str):
    """Number-input on_change: the person edited a value — convert what
    they typed (in whatever unit is currently shown) back to the
    canonical mol/kgw so it survives future unit switches correctly."""
    def _cb():
        unit = UNIT_LABELS[st.session_state["conc_unit_select"]]
        st.session_state[f"canon_{ion}"] = _to_molkgw(st.session_state[f"{ion}_display"], unit, ion)
    return _cb


# First run: seed the display values from the (zero) canonical values.
for _ion in IONS_UI:
    st.session_state.setdefault(
        f"{_ion}_display",
        _from_molkgw(st.session_state[f"canon_{_ion}"],
                     UNIT_LABELS[st.session_state["conc_unit_select"]], _ion),
    )

unit_col1, unit_col2 = st.sidebar.columns([2, 3], vertical_alignment="center")
with unit_col1:
    st.markdown("Units")
with unit_col2:
    conc_unit = st.selectbox(
        "Units", list(UNIT_LABELS.keys()),
        label_visibility="collapsed", key="conc_unit_select",
        on_change=_sync_display_from_canon,
        help="mol/kgw is used directly with no conversion. mol/L and ppm are both "
             "converted assuming solution mass \u2248 water mass \u2014 exact for ppm-level "
             "concentrations, an approximation for mol/L at high concentration. Switching "
             "units converts the displayed numbers; it does not reinterpret them.",
    )
unit_key = UNIT_LABELS[conc_unit]
unit_step = {"mol/kgw": 0.01, "mol/L": 0.01, "ppm": 50.0}[unit_key]

Li_val_disp = st.sidebar.number_input("Li\u207a", min_value=0.0, step=unit_step, format="%.4f",
                                       key="Li_display", on_change=_sync_canon_from_display("Li"))
Na_val_disp = st.sidebar.number_input("Na\u207a", min_value=0.0, step=unit_step, format="%.4f",
                                       key="Na_display", on_change=_sync_canon_from_display("Na"))
Mg_val_disp = st.sidebar.number_input("Mg\u00b2\u207a", min_value=0.0, step=unit_step, format="%.4f",
                                       key="Mg_display", on_change=_sync_canon_from_display("Mg"))
st.sidebar.caption("Cl\u207b is computed automatically from electroneutrality.")

# Canonical mol/kgw values (kept in sync by the callbacks above) are what
# actually get sent to the calculation, regardless of which unit is showing.
Li_val, Na_val, Mg_val = (st.session_state["canon_Li"], st.session_state["canon_Na"],
                           st.session_state["canon_Mg"])

st.sidebar.subheader("Process parameters")
Y = st.sidebar.slider("Li recovery, Y", min_value=0.01, max_value=0.99, value=0.90, step=0.01)
CF = st.sidebar.number_input("Concentration factor, CF", min_value=1.01, value=30.0, step=1.0)
S_Na = st.sidebar.number_input("Selectivity S (Li/Na)", min_value=0.01, value=100.0, step=1.0,
                                help="Larger = better rejection of Na from the product.")
S_Mg = st.sidebar.number_input("Selectivity S (Li/Mg)", min_value=0.01, value=100.0, step=1.0,
                                help="Larger = better rejection of Mg from the product.")

st.sidebar.subheader("Database")
db_choice = st.sidebar.radio("PHREEQC database", ["Built-in pitzer.dat", "Upload custom .dat"],
                              label_visibility="collapsed")
db_name, db_dir = "pitzer.dat", None
if db_choice == "Upload custom .dat":
    uploaded_db = st.sidebar.file_uploader("Upload database file (.dat)", type=["dat"])
    if uploaded_db is not None:
        tmp_dir = tempfile.mkdtemp()
        db_path = Path(tmp_dir) / uploaded_db.name
        db_path.write_bytes(uploaded_db.getvalue())
        db_name, db_dir = uploaded_db.name, tmp_dir
        st.sidebar.success(f"Using {uploaded_db.name}")
    else:
        st.sidebar.info("No file uploaded yet — falling back to the built-in database.")

run_clicked = st.sidebar.button("Run calculation", type="primary", width="stretch")

st.sidebar.markdown("<br>", unsafe_allow_html=True)
st.sidebar.markdown(
    """
    <span style="color:#1e5473">
    The following calculator is provided by the Yip Lab's Dr. Yuren Feng, Anna Seeley
    McGillis [PhD], Lucas Caldentey, and Dr. Ngai Yin Yip
    </span>
    """,
    unsafe_allow_html=True,
)

# ----------------------------------------------------------------------------
# Main panel
# ----------------------------------------------------------------------------
st.title("Minimum Separation Energy ($E_{\\mathrm{min}}$) Calculator")
st.caption(
    "This site serves as a single-point calculation of the thermodynamic minimum energy "
    "to concentrate lithium from a Li\u207a\u2013Na\u207a\u2013Mg\u00b2\u207a brine, via "
    "PHREEQC/Pitzer activities and various energy-balance equations. Any oversaturated "
    "salt is precipitated and the system is re-equilibrated automatically."
)

if not run_clicked:
    st.info("Set your inputs in the sidebar and click **Run calculation**.")
    st.stop()

inp = CalcInputs(
    Li_val=Li_val, Na_val=Na_val, Mg_val=Mg_val, conc_unit="mol/kgw",
    Y=Y, CF=CF, S_Na=S_Na, S_Mg=S_Mg,
    database=db_name, database_dir=db_dir,
)

with st.spinner("Running PHREEQC..."):
    result = run_calculation(inp)

if not result.success:
    st.error(result.error)
    st.stop()

for w in result.warnings:
    st.warning(w)

# --- headline metrics (emphasized) ---------------------------------------
c1, c2 = st.columns(2)
c1.metric("$E_{\\mathrm{min}}$", f"{result.E_min:.3f} kJ/mol Li")
c2.metric("Product ionic strength", f"{result.I_product:.2f} mol/kgw")

st.divider()

# --- species contribution chart -----------------------------------------
left, right = st.columns([3, 2])

with left:
    st.subheader("Species contributions to $E_{\\mathrm{min}}$")
    items = sorted(result.E_components.items(), key=lambda kv: -abs(kv[1]))
    labels = [k for k, _ in items]
    values = [v for _, v in items]
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in values]
    fig = go.Figure(go.Bar(x=values, y=labels, orientation="h", marker_color=colors,
                            text=[f"{v:+.2f}" for v in values], textposition="outside"))
    fig.update_layout(xaxis_title="kJ / mol Li recovered", height=320,
                       margin=dict(l=10, r=10, t=10, b=10))
    st.plotly_chart(fig, width="stretch")
    st.caption("Positive (red) = energy cost; negative (blue) = energy credit. "
               "Bars should sum to $E_{\\mathrm{min}}$ above (the pipeline checks this internally).")

with right:
    st.subheader("Precipitated solids")
    if result.solids:
        rows = "\n".join(f"| {ph} | {v:.4f} |" for ph, v in result.solids.items())
        st.markdown("| Phase | mol / mol Li recovered |\n|---|---|\n" + rows)
    else:
        st.write("None \u2014 all candidate phases stayed under-saturated (SI \u2264 0).")

    st.markdown("**Product stream composition**")
    prod = result.streams["Product (mol/kgw)"]
    ion_rows = "\n".join(
        f"| {ion} | {prod.get(ion, 0.0):.4f} |"
        for ion in ("Li", "Na", "Mg", "Cl") if ion in prod
    )
    st.markdown("| Ion | mol/kgw |\n|---|---|\n" + ion_rows)
    st.caption("Same product-stream values used in the plot and chart above, "
               "shown here as a quick reference for how much of each ion is present.")

# --- secondary metrics (demoted, shown below the plot) -------------------
d1, d2, d3 = st.columns(3)
d1.metric("$E_{\\mathrm{min,LW}}$ (LiCl-only ref.)", f"{result.E_min_LW:.3f} kJ/mol Li",
          help="Simplified reference value treating the brine as pure LiCl (ignores Na, Mg, "
               "and water). Not expected to match $E_{\\mathrm{min}}$ above except for a pure-LiCl feed.")
d2.metric("Feed ionic strength", f"{result.I_feed:.2f} mol/kgw")
d3.metric("PHREEQC simulations run", result.n_simulations)

st.divider()

# --- stream compositions ---------------------------------------------------
st.subheader("Stream compositions (molal, mol/kgw)")
ion_order = [c for c in ("Li", "Na", "Mg", "Cl")
             if any(c in comp for comp in result.streams.values())]
header = "| Stream | " + " | ".join(ion_order) + " | H2O (kg / kg feed water) |"
sep = "|" + "---|" * (len(ion_order) + 2)
body = "\n".join(
    "| " + name + " | " + " | ".join(f"{comp.get(ion, 0.0):.4f}" for ion in ion_order)
    + f" | {comp.get('H2O_kgw', 0.0):.4f} |"
    for name, comp in result.streams.items()
)
st.markdown("\n".join([header, sep, body]))
st.caption("Water is shown as mass relative to the feed's (1.0 for the Feed itself, "
           "Y/CF for the Product, 1\u2212Y/CF for the Retentate) \u2014 how the feed's water "
           "splits across the three streams.")

# --- download ---------------------------------------------------------
st.divider()
import csv
import io

buf = io.StringIO()
writer = csv.writer(buf)
writer.writerow(["Component", "E_i (kJ/mol Li)"])
for k, v in result.E_components.items():
    writer.writerow([k, f"{v:.6f}"])
writer.writerow(["TOTAL E_min", f"{result.E_min:.6f}"])
st.download_button(
    "Download results as CSV",
    buf.getvalue().encode("utf-8"),
    file_name="emin_results.csv",
    mime="text/csv",
)

with st.expander("What is $E_{\\mathrm{min}}$, and what do these inputs mean?"):
    st.markdown("""
- **$E_{\\mathrm{min}}$** is the least energy that thermodynamically allows for separating lithium
  from this particular brine at the specified recovery and concentration \u2014 it serves more as a
  benchmark rather than a prediction of what any real process will achieve.
- **$E_{\\mathrm{min,LW}}$** is a simplified variation on $E_{\\mathrm{min}}$ that assumes the brine
  is pure LiCl, ignoring Na, Mg, and the water term.
- **Y (recovery)** is the fraction of feed lithium captured in the product.
- **CF (concentration factor)** is how many times more concentrated the product's
  lithium is versus the feed.
- **S (selectivity)** describes how well each impurity is excluded from the product
  relative to lithium. Larger S means better rejection.
- **Units**: mol/kgw (molality) is used exactly as entered. mol/L and ppm are both
  converted assuming the feed's water mass \u2248 its solution mass.
- If any candidate salt (currently halite, bischofite) becomes oversaturated in the
  product stream at the requested CF, it is precipitated out and the system is
  re-equilibrated \u2014 shown above as an extra simulation and a warning banner.
    """)
