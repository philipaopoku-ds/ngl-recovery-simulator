"""Streamlit app: interactive turbo expander NGL recovery simulator."""
import altair as alt
import pandas as pd
import streamlit as st
from ngl_plant import simulate, Inputs, DEFAULT_FEED, LABELS

st.set_page_config(page_title="NGL Recovery Simulator", page_icon="🏭", layout="wide")

st.title("Turbo Expander NGL Recovery Simulator")
st.caption("A simplified steady-state model of a generic gas plant, built in Python with the "
           "Peng-Robinson equation of state. Move the sliders to see how operating choices "
           "change product recovery and energy use. Educational model; not a design tool.")

with st.sidebar:
    st.header("Operating conditions")
    feed = st.slider("Feed gas rate (MMscf/d)", 25, 300, 100, 5)
    P_in = st.slider("Inlet pressure (bar)", 40, 90, 65, 1)
    T_cold = st.slider("Cold separator temperature (°C)", -55, -10, -35, 1)
    P_col = st.slider("Demethanizer pressure (bar)", 12, 38, 22, 1,
                      help="The expander lets the gas down to this pressure.")
    eff = st.slider("Expander efficiency (%)", 60, 90, 82, 1) / 100
    st.header("Feed composition (mol %)")
    st.caption("Sweet, dry gas. Values are normalised to 100%.")
    feed_pct = {}
    cols = st.columns(2)
    for i, k in enumerate(LABELS):
        feed_pct[k] = cols[i % 2].number_input(k, 0.0, 100.0, float(DEFAULT_FEED[k]), 0.1, key=k)

if P_col >= P_in - 2:
    st.error("The demethanizer pressure must be below the inlet pressure so the expander can do work. "
             "Lower the demethanizer pressure or raise the inlet pressure.")
    st.stop()


@st.cache_data(show_spinner="Solving the flowsheet...")
def run(feed, P_in, T_cold, P_col, eff, feed_items):
    r = simulate(Inputs(feed_mmscfd=feed, P_inlet_bar=P_in, T_cold_C=T_cold,
                        P_column_bar=P_col, expander_eff=eff, P_sales_bar=P_in,
                        feed_mol_pct=dict(feed_items)))
    r.pop("inputs")
    return r


r = run(feed, P_in, T_cold, P_col, eff, tuple(feed_pct.items()))

c1, c2, c3, c4 = st.columns(4)
c1.metric("Ethane recovery", f"{r['recovery_pct']['C2']:.1f} %")
c2.metric("Propane recovery", f"{r['recovery_pct']['C3']:.1f} %")
c3.metric("NGL product", f"{r['ngl_tpd']:,.0f} t/d")
c4.metric("Net compression", f"{r['net_compression_kW']/1000:.2f} MW")

left, right = st.columns(2)
with left:
    st.subheader("Recovery by component")
    rec = pd.DataFrame({"Component": LABELS[2:], "Recovery (%)": [r["recovery_pct"][k] for k in LABELS[2:]]})
    chart = (alt.Chart(rec).mark_bar(color="#3B5FD0", cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
             .encode(x=alt.X("Component", sort=None, axis=alt.Axis(labelAngle=0)),
                     y=alt.Y("Recovery (%)", scale=alt.Scale(domain=[0, 100])),
                     tooltip=["Component", alt.Tooltip("Recovery (%)", format=".1f")]))
    st.altair_chart(chart, use_container_width=True)
with right:
    st.subheader("Key temperatures and energy")
    st.dataframe(pd.DataFrame({
        "Value": [f"{r['T_expander_out_C']:.1f} °C", f"{r['T_column_top_C']:.1f} °C",
                  f"{r['T_reboiler_C']:.1f} °C", f"{r['expander_power_kW']/1000:.2f} MW",
                  f"{r['recompression_power_kW']/1000:.2f} MW", f"{r['inlet_cooling_kW']/1000:.2f} MW",
                  f"{r['reboiler_duty_kW']/1000:.2f} MW", f"{r['residue_mmscfd']:.1f} MMscf/d"]},
        index=["Expander outlet", "Demethanizer top", "Reboiler (NGL bubble point)",
               "Expander power recovered", "Residue gas recompression", "Inlet cooling duty",
               "Reboiler duty", "Residue gas to sales"]), use_container_width=True)

st.subheader("Product compositions (mol %)")
st.dataframe(pd.DataFrame({"Feed": [feed_pct[k] / sum(feed_pct.values()) * 100 for k in LABELS],
                           "Residue gas": [r["residue_comp"][k] for k in LABELS],
                           "NGL product": [r["ngl_comp"][k] for k in LABELS]}, index=LABELS).round(2),
             use_container_width=True)

with st.expander("How the model works and its limits"):
    st.markdown("""
- **Thermodynamics:** Peng-Robinson equation of state with standard binary interaction parameters (`thermo` library).
- **Cold separator:** isothermal flash after inlet cooling (0.5 bar pressure drop).
- **Turbo expander:** isentropic expansion corrected by the expander efficiency; its power drives the booster compressor.
- **Demethanizer (simplified):** the expander and JT valve outlets meet on a cold top stage. Top vapour leaves as residue gas; the liquid is stripped so the NGL meets a methane-to-ethane ratio of 0.02.
- **Recompression:** residue gas warmed to 25 °C, then compressed back to inlet pressure (78% efficiency).
- **Limits:** no rigorous tray-by-tray column, reflux or side reboilers, and no hydrate or CO₂ freezing checks. A full design would use a rigorous simulator such as Aspen HYSYS. Results show trends and trade-offs, not design values.
""")
st.caption("Built by Philipa Opoku · Generic feed data · Not based on any real plant")
