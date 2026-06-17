import streamlit as st
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="Water Uptake Model", layout="wide")

st.title("Water Uptake Correction Model")
st.markdown("**Trial: CMVN + EG, Membrane 1 (CMVN+EG-1)**")

# ── Data source explanation ──────────────────────────────────────────────────
with st.expander("Data source", expanded=True):
    st.markdown("""
| Variable | Value | Source |
|---|---|---|
| **m_dry_initial** | 0.0867 g | "CMVN vs FKS50" sheet, Row: :CMVN+EG-1:, Column *m dry* (original state) |
| **m_wet** | 0.0961 g | CMVN vs FKS50" sheet, Row: :CMVN+EG-1, column *m wet* (after sorption) |
| **m_dry** | 0.0857 g | CMVN vs FKS50" sheet, Row: :CMVN+EG-1, column *m dry* (after sorption) |
| **m_solute** | 0.0000719 g (default) | CMVN vs FKS50" sheet, Row: :CMVN+EG-1, *EG/PEGs mass in membrane* = 0.0719 mg |

The updated WU formula is:

**WU = [(m_wet − m_dry) − m_solute × (1 − f_residual)] / m_dry_initial**
    """)

# ── Fixed constants from data ────────────────────────────────────────────────
M_DRY_INITIAL = 0.0867
M_WET         = 0.0961
M_DRY         = 0.0857
DEFAULT_M_SOLUTE_MG = 0.0719  # mg

# ── Sidebar controls ──────────────────────────────────────────────────────────
st.sidebar.header("Controls")

# m_solute input
st.sidebar.subheader("m_solute")
st.sidebar.markdown("Default value from CMVN EG M1 desorption TOC data.")
m_solute_mg = st.sidebar.number_input(
    "m_solute (mg)",
    min_value=0.0,
    max_value=100.0,
    value=DEFAULT_M_SOLUTE_MG,
    step=0.001,
    format="%.4f",
    help="Mass of solute inside the membrane. Default = 0.0719 mg from CMVN EG M1 desorption TOC."
)
m_solute = m_solute_mg / 1000  # convert mg to g

# f_residual slider
st.sidebar.subheader("f_residual")
f_residual = st.sidebar.slider(
    "f_residual",
    min_value=0.0,
    max_value=1.0,
    value=0.0,
    step=0.01,
)

# ── Calculations ──────────────────────────────────────────────────────────────
wu_original        = (M_WET - M_DRY) / M_DRY_INITIAL
m_solute_evap      = m_solute * (1 - f_residual)
m_water_corrected  = (M_WET - M_DRY) - m_solute_evap
wu_corrected       = m_water_corrected / M_DRY_INITIAL
delta_wu           = wu_corrected - wu_original

# ── Metrics ───────────────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)
col1.metric(
    "Original WU (uncorrected)",
    f"{wu_original*100:.4f}%",
    help="(m_wet − m_dry) / m_dry_initial — no evaporation correction applied"
)
col2.metric(
    "Corrected WU",
    f"{wu_corrected*100:.4f}%",
    delta=f"{delta_wu*100:+.4f}%",
    help="Full corrected formula applying f_residual evaporation adjustment"
)
col3.metric(
    "m_solute evaporated",
    f"{m_solute_evap*1000:.5f} mg",
    help="m_solute × (1 − f_residual) — portion of solute lost during vacuum drying"
)
col4.metric(
    "f_residual selected",
    f"{f_residual:.2f}"
)

st.markdown("---")

# ── Full curve: WU vs f_residual ──────────────────────────────────────────────
f_values  = np.linspace(0, 1, 500)
wu_values = [((M_WET - M_DRY) - m_solute * (1 - f)) / M_DRY_INITIAL * 100
             for f in f_values]

fig = go.Figure()

# Corrected WU curve
fig.add_trace(go.Scatter(
    x=f_values, y=wu_values,
    mode="lines",
    name="Corrected WU",
    line=dict(color="#1D9E75", width=2.5),
    hovertemplate="f_residual: %{x:.2f}<br>WU: %{y:.5f}%<extra></extra>"
))

# Original uncorrected WU horizontal reference line
fig.add_hline(
    y=wu_original * 100,
    line_dash="dash",
    line_color="#D85A30",
    line_width=2,
    annotation_text=f"Original WU = {wu_original*100:.4f}%",
    annotation_position="top right",
    annotation_font=dict(color="#D85A30", size=12)
)

# Marker for current slider position
fig.add_trace(go.Scatter(
    x=[f_residual],
    y=[wu_corrected * 100],
    mode="markers",
    name=f"Current f_residual = {f_residual:.2f}",
    marker=dict(color="#D85A30", size=12, symbol="circle",
                line=dict(color="white", width=2)),
    hovertemplate=f"f_residual = {f_residual:.2f}<br>WU = {wu_corrected*100:.5f}%<extra></extra>"
))

# Reference lines for known solute extremes
fig.add_vline(x=0.0, line_dash="dot", line_color="#888780", line_width=1.5,
              annotation_text="EG (f=0)", annotation_position="top left",
              annotation_font=dict(size=11, color="#888780"))
fig.add_vline(x=1.0, line_dash="dot", line_color="#888780", line_width=1.5,
              annotation_text="PEG 1000 (f≈1)", annotation_position="top right",
              annotation_font=dict(size=11, color="#888780"))

fig.update_layout(
    title=dict(
        text="Corrected Water Uptake vs f_residual — CMVN+EG-1",
        font=dict(size=16)
    ),
    xaxis=dict(
        title="f_residual (fraction of solute remaining after vacuum drying)",
        tickformat=".2f",
        range=[-0.02, 1.02],
        showgrid=True, gridcolor="#EEEEEE"
    ),
    yaxis=dict(
        title="Water Uptake (%)",
        tickformat=".4f",
        showgrid=True, gridcolor="#EEEEEE"
    ),
    legend=dict(
        x=0.02, y=0.98,
        bgcolor="rgba(255,255,255,0.85)",
        bordercolor="#DDDDDD", borderwidth=1
    ),
    hovermode="x unified",
    height=480,
    plot_bgcolor="black",
    paper_bgcolor="black",
    margin=dict(t=50, b=50, l=60, r=30)
)

st.plotly_chart(fig, use_container_width=True)

# ── Bar chart: f=0 vs f=1 comparison ─────────────────────────────────────────

wu_at_f0 = ((M_WET - M_DRY) - m_solute * 1.0) / M_DRY_INITIAL * 100
wu_at_f1 = ((M_WET - M_DRY) - m_solute * 0.0) / M_DRY_INITIAL * 100
diff = abs(wu_at_f1 - wu_at_f0)


# ── Step-by-step calculation ───────────────────────────────────────────────────
st.markdown("---")
st.subheader("Step-by-Step Calculation at Selected f_residual")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("**Fixed values from CMVN+EG-1 trial**")
    st.latex(r"m_{dry\_initial} = 0.0867 \text{ g}")
    st.latex(r"m_{wet} = 0.0961 \text{ g}")
    st.latex(r"m_{dry} = 0.0857 \text{ g}")
    st.latex(rf"m_{{solute}} = {m_solute_mg:.4f} \text{{ mg}} = {m_solute:.7f} \text{{ g}}")

with col_b:
    st.markdown("**Calculation steps**")
    st.latex(rf"f_{{residual}} = {f_residual:.2f}")
    st.latex(rf"m_{{solute\_evap}} = {m_solute_mg:.4f} \times (1 - {f_residual:.2f}) = {m_solute_evap*1000:.5f} \text{{ mg}}")
    st.latex(rf"m_{{water}} = (0.0961 - 0.0857) - {m_solute_evap*1000:.5f}/1000 = {m_water_corrected*1000:.5f} \text{{ mg}}")
    st.latex(rf"WU = \frac{{{m_water_corrected*1000:.5f}/1000}}{{0.0867}} = {wu_corrected*100:.5f}\%")

# ── Full range summary ─────────────────────────────────────────────────────────
st.markdown("---")
st.subheader("Full Range Summary")

c1, c2, c3 = st.columns(3)
c1.metric("WU at f_residual = 0", f"{wu_at_f0:.5f}%",
          help="Maximum correction — all solute evaporates during drying")
c2.metric("WU at f_residual = 1", f"{wu_at_f1:.5f}%",
          help="No correction — all solute stays in membrane after drying")
c3.metric("Total range", f"{diff:.5f} % points",
          help="Maximum possible WU variation due to evaporation correction")
