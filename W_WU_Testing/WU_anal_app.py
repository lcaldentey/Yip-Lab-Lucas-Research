import streamlit as st
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="Water Uptake Model", layout="wide")

st.title("Water Uptake Correction Model")
st.markdown("**Trial: CMVN + EG, Membrane 1 (CMVN+EG-1)**")

# ── Data source explanation ──────────────────────────────────────────────────
with st.expander("Data source — how each value was found", expanded=True):
    st.markdown("""
| Variable | Value | Source |
|---|---|---|
| **m_dry_initial** | 0.0867 g | Water uptake sheet, CMVN+EG-1 row, column *m dry* (original state) |
| **m_wet** | 0.0961 g | Water uptake sheet, CMVN+EG-1 row, column *m wet* (after sorption) |
| **m_dry** | 0.0857 g | Water uptake sheet, CMVN+EG-1 row, column *m dry* (after sorption) |
| **m_solute** | 0.0000719 g | CMVN vs FKS50 sheet, CMVN EG M1 row, *EG/PEGs mass in membrane* = 0.0719 mg |

All four values come from the same membrane piece and the same experimental trial.
The corrected WU formula is:

**WU = [(m_wet − m_dry) − m_solute × (1 − f_residual)] / m_dry_initial**
    """)

# ── Fixed constants from data ────────────────────────────────────────────────
M_DRY_INITIAL = 0.0867   # g
M_WET         = 0.0961   # g  (post-sorption)
M_DRY         = 0.0857   # g  (post-sorption)
M_SOLUTE      = 0.0000719 # g  (0.0719 mg)

# ── Uncorrected WU (original thesis formula, no correction) ─────────────────
wu_original = (M_WET - M_DRY) / M_DRY_INITIAL

# ── Sidebar — single f_residual slider ───────────────────────────────────────
st.sidebar.header("Controls")
f_residual = st.sidebar.slider(
    "f_residual  (fraction of solute remaining after vacuum drying)",
    label_visibility="collapsed",
    min_value=0.0,
    max_value=1.0,
    value=0.0,        # EG is fully volatile → f_residual = 0
    step=0.01,
    help="0 = fully volatile (all solute evaporates, e.g. EG)  |  1 = fully non-volatile (solute stays, e.g. PEG 1000)"
)

# ── Corrected WU at selected f_residual ─────────────────────────────────────
m_solute_evaporated = M_SOLUTE * (1 - f_residual)
m_water_corrected   = (M_WET - M_DRY) - m_solute_evaporated
wu_corrected        = m_water_corrected / M_DRY_INITIAL
delta_wu            = wu_corrected - wu_original

# ── Metrics row ──────────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)
col1.metric("Original WU (uncorrected)", f"{wu_original*100:.4f} %",
            help="(m_wet − m_dry) / m_dry_initial — no evaporation correction")
col2.metric("Corrected WU", f"{wu_corrected*100:.4f} %",
            delta=f"{delta_wu*100:+.4f} %",
            help="Applies f_residual correction to numerator")
col3.metric("m_solute evaporated", f"{m_solute_evaporated*1000:.5f} mg",
            help="m_solute × (1 − f_residual)")
col4.metric("f_residual selected", f"{f_residual:.2f}")

st.markdown("---")

# ── Full curve: WU vs f_residual ─────────────────────────────────────────────
f_values  = np.linspace(0, 1, 500)
wu_values = [((M_WET - M_DRY) - M_SOLUTE * (1 - f)) / M_DRY_INITIAL * 100
             for f in f_values]

fig = go.Figure()

# Corrected WU curve
fig.add_trace(go.Scatter(
    x=f_values, y=wu_values,
    mode="lines",
    name="Corrected WU",
    line=dict(color="#1D9E75", width=2.5)
))

# Original (uncorrected) WU horizontal line
fig.add_hline(
    y=wu_original * 100,
    line_dash="dash",
    line_color="#D85A30",
    annotation_text=f"Original WU = {wu_original*100:.4f}%",
    annotation_position="top right",
    annotation_font_color="#D85A30"
)

# Marker for current slider value
fig.add_trace(go.Scatter(
    x=[f_residual],
    y=[wu_corrected * 100],
    mode="markers",
    name=f"Current f_residual = {f_residual:.2f}",
    marker=dict(color="#D85A30", size=12, symbol="circle")
))

# Vertical reference lines for known solutes
fig.add_vline(x=0.0, line_dash="dot", line_color="#888780",
              annotation_text="EG (f=0)", annotation_position="top left",
              annotation_font_size=11)
fig.add_vline(x=1.0, line_dash="dot", line_color="#888780",
              annotation_text="PEG 1000 (f≈1)", annotation_position="top right",
              annotation_font_size=11)

fig.update_layout(
    title=dict(
        text="Water Uptake vs f_residual — CMVN+EG-1",
        font=dict(size=16)
    ),
    xaxis=dict(
        title="f_residual (fraction of solute remaining after drying)",
        tickformat=".2f",
        range=[-0.02, 1.02]
    ),
    yaxis=dict(
        title="Water Uptake (%)",
        tickformat=".4f"
    ),
    legend=dict(x=0.02, y=0.98),
    hovermode="x unified",
    height=480,
    plot_bgcolor="black",
    paper_bgcolor="black"
)
fig.update_xaxes(showgrid=True, gridcolor="#eeeeee")
fig.update_yaxes(showgrid=True, gridcolor="#eeeeee")

st.plotly_chart(fig, use_container_width=True)

# ── Step-by-step calculation breakdown ───────────────────────────────────────
st.subheader("Step-by-step calculation at selected f_residual")

col_a, col_b = st.columns(2)

with col_a:
    st.markdown("**Fixed values from CMVN+EG-1 trial**")
    st.latex(r"m_{dry\_initial} = 0.0867 \text{ g}")
    st.latex(r"m_{wet} = 0.0961 \text{ g}")
    st.latex(r"m_{dry} = 0.0857 \text{ g}")
    st.latex(r"m_{solute} = 0.0000719 \text{ g}")

with col_b:
    st.markdown("**Calculation steps**")
    st.latex(rf"f_{{residual}} = {f_residual:.2f}")
    st.latex(rf"m_{{solute\_evap}} = 0.0000719 \times (1 - {f_residual:.2f}) = {m_solute_evaporated*1000:.5f} \text{{ mg}}")
    st.latex(rf"m_{{water}} = (0.0961 - 0.0857) - {m_solute_evaporated*1000:.5f}/1000 = {m_water_corrected*1000:.5f} \text{{ mg}}")
    st.latex(rf"WU = \frac{{{m_water_corrected*1000:.5f}/1000}}{{0.0867}} = {wu_corrected*100:.4f}\%")

# ── WU range summary ─────────────────────────────────────────────────────────
st.markdown("---")
st.subheader("Full range summary")

wu_at_0 = ((M_WET - M_DRY) - M_SOLUTE * 1.0) / M_DRY_INITIAL * 100
wu_at_1 = ((M_WET - M_DRY) - M_SOLUTE * 0.0) / M_DRY_INITIAL * 100

col_x, col_y, col_z = st.columns(3)
col_x.metric("WU at f_residual = 0 (EG, fully volatile)",
             f"{wu_at_0:.4f} %")
col_y.metric("WU at f_residual = 1 (PEG 1000, non-volatile)",
             f"{wu_at_1:.4f} %")
col_z.metric("Total WU range across all f_residual",
             f"{abs(wu_at_1 - wu_at_0)*100:.4f} % points")
