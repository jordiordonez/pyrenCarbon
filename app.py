import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from carboni.creixement import carregar_models, simular
from carboni.equacions import CO2_PER_C, ESPECIES, HABITATS, NOMS_COMUNS
from carboni.intervencio import TIPUS, aclarir, area_basimetrica
from carboni.inventari import (
    carregar_inf1and,
    pinars,
    resum_parceles,
    rodal_analog,
    rodal_de_csv,
    rodal_de_parcela,
)

ARREL = Path(__file__).resolve().parent
DATA = ARREL / "data"

COLOR_SENSE = "#2F6B4F"
COLOR_AMB = "#C7772B"
COLOR_EXTRET = "#B9573A"
COLOR_INICIAL = "#9AA59F"
VARIABLES = {
    "co2": ("CO₂ emmagatzemat (aeri + subterrani)", "t CO₂/ha"),
    "carboni": ("Carboni total", "t C/ha"),
    "cat": ("Carboni aeri", "t C/ha"),
    "increment_brut": ("Creixement anual en carboni", "t C/ha·any"),
    "ab": ("Àrea basimètrica", "m²/ha"),
    "densitat": ("Densitat", "peus/ha"),
    "dq": ("Diàmetre quadràtic mig", "cm"),
    "ho": ("Alçada dominant", "m"),
}
NOMS_ESPECIE = {f"{NOMS_COMUNS[c]} ({ESPECIES[c]})": c for c in ESPECIES}

st.set_page_config(page_title="Carboni a les pinedes d'Andorra", page_icon="🌲", layout="wide")
st.markdown(
    """
    <style>
    .block-container {padding-top: 2rem; max-width: 1250px;}
    h1 {font-weight: 700; letter-spacing: -0.02em;}
    div[data-testid="stMetric"] {background: #FFFFFF; border: 1px solid #E3DED1; border-radius: 12px; padding: 14px 16px;}
    div[data-testid="stMetricLabel"] p {font-size: 0.85rem; color: #5B665F;}
    .subtitol {color: #5B665F; font-size: 1.05rem; margin-top: -0.6rem; margin-bottom: 1.2rem;}
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_resource
def carregar():
    parceles, arbres = carregar_inf1and(DATA / "inf1and")
    resum = resum_parceles(parceles, arbres)
    models = carregar_models(DATA / "calibracio.json")
    calibracio = json.loads((DATA / "calibracio.json").read_text(encoding="utf-8"))
    return resum, arbres, models, calibracio


@st.cache_data
def estoc_actual_pinars() -> pd.DataFrame:
    resum, arbres, models, _ = carregar()
    seleccio = pinars(resum).copy()
    seleccio["co2"] = [
        simular(rodal_de_parcela(resum, arbres, c), models, 0).serie["co2"].iloc[0] for c in seleccio["parcela"]
    ]
    return seleccio


@st.cache_data(show_spinner="Simulant el creixement…")
def simular_cache(arbres: pd.DataFrame, habitat: str, altitud: float, descripcio: str, anys: int):
    from carboni.inventari import Rodal

    _, _, models, _ = carregar()
    return simular(Rodal(arbres, habitat, altitud, descripcio), models, anys)


resum, arbres_inf, models, calibracio = carregar()


def fmt(valor: float, decimals: int = 1) -> str:
    return f"{valor:,.{decimals}f}".replace(",", " ").replace(".", ",")


def etiqueta_parcela(fila) -> str:
    return (f"{fila.parcela} · {NOMS_COMUNS[fila.especie_dominant]} · {fila.habitat} · {fila.altitud:.0f} m · "
            f"{fila.n:.0f} peus/ha · {fila.ab:.0f} m²/ha")


# ---------- Capçalera ----------
st.title("🌲 Carboni a les pinedes d'Andorra")
st.markdown(
    '<p class="subtitol">Estoc i captació de CO₂ de pinedes de pi roig i pi negre, amb creixement arbre a arbre '
    "calibrat amb l'Inventari Forestal Nacional i les parcel·les del 1r Inventari Nacional Forestal d'Andorra.</p>",
    unsafe_allow_html=True,
)

# ---------- Barra lateral ----------
with st.sidebar:
    st.header("1 · Rodal de partida")
    origen = st.radio(
        "Origen de les dades",
        ["Parcel·la de l'INF d'Andorra", "Dades del rodal", "Inventari propi (CSV)"],
        label_visibility="collapsed",
    )
    rodal = None
    error = None
    codi_seleccionat = None
    if origen == "Parcel·la de l'INF d'Andorra":
        filtre = st.segmented_control("Espècie", ["Totes", "Pi negre", "Pi roig"], default="Totes") or "Totes"
        sp_filtre = {"Pi negre": "Pu", "Pi roig": "Ps"}.get(filtre)
        opcions = pinars(resum, sp_filtre).sort_values("parcela")
        fila = st.selectbox(
            "Parcel·la", list(opcions.itertuples(index=False)), format_func=etiqueta_parcela,
            help="Pinars (≥70 % de l'àrea basimètrica de pi roig o pi negre) del 1r INF d'Andorra, 2022.",
        )
        codi_seleccionat = fila.parcela
        rodal = rodal_de_parcela(resum, arbres_inf, codi_seleccionat)
    elif origen == "Dades del rodal":
        sp = NOMS_ESPECIE[st.selectbox("Espècie", list(NOMS_ESPECIE))]
        habitat = st.radio("Hàbitat", list(HABITATS), format_func=HABITATS.get, horizontal=True)
        altitud = st.number_input("Altitud (m)", 900, 2500, 1900 if sp == "Pu" else 1500, step=50)
        c1, c2 = st.columns(2)
        n = c1.number_input("Densitat (peus/ha)", 50, 5000, 900, step=50)
        ab = c2.number_input("Àrea basimètrica (m²/ha)", 2.0, 90.0, 35.0, step=1.0)
        rodal = rodal_analog(resum, arbres_inf, sp, habitat, float(altitud), float(n), float(ab))
    else:
        plantilla = pd.DataFrame({"especie": ["Pu", "Pu", "Ps"], "dbh": [18.5, 32.0, 24.2], "ht": [11.0, None, 14.5]})
        st.download_button("Descarregar plantilla", plantilla.to_csv(index=False).encode(), "plantilla_inventari.csv",
                           "text/csv")
        fitxer = st.file_uploader("Inventari (CSV)", type="csv",
                                  help="Columnes: especie (Ps/Pu o nom), dbh en cm, ht opcional, peus_ha opcional.")
        habitat = st.radio("Hàbitat", list(HABITATS), format_func=HABITATS.get, horizontal=True)
        altitud = st.number_input("Altitud (m)", 900, 2500, 1900, step=50)
        superficie = st.number_input("Superfície de la parcel·la (m²)", 50.0, 5000.0, 452.0, step=10.0,
                                     help="Només si el fitxer no té la columna peus_ha. 452 m² = radi de 12 m.")
        if fitxer is not None:
            try:
                rodal = rodal_de_csv(pd.read_csv(fitxer, sep=None, engine="python"), habitat, float(altitud),
                                     superficie)
            except ValueError as e:
                error = str(e)

    st.header("2 · Intervenció")
    amb_intervencio = st.toggle("Comparar amb una aclarida", value=True)
    tipus = st.selectbox("Tipus d'aclarida", list(TIPUS), format_func=TIPUS.get, disabled=not amb_intervencio)
    pes = st.slider("Àrea basimètrica extreta (%)", 5, 60, 30, step=5, disabled=not amb_intervencio)

    st.header("3 · Horitzó")
    anys = st.slider("Anys de projecció", 5, 80, 30, step=5)

if rodal is None:
    st.info(error or "Puja un inventari en CSV per començar (pots descarregar la plantilla a la barra lateral).")
    st.stop()

# ---------- Simulació ----------
sim = simular_cache(rodal.arbres, rodal.habitat, rodal.altitud, rodal.descripcio, anys)
serie = sim.serie
sim_amb, extrets = None, None
if amb_intervencio:
    rodal_aclarit, extrets = aclarir(rodal, tipus, pes / 100.0)
    sim_amb = simular_cache(rodal_aclarit.arbres, rodal_aclarit.habitat, rodal_aclarit.altitud,
                            rodal_aclarit.descripcio, anys)

st.caption(f"**{rodal.descripcio}** · hàbitat {HABITATS[rodal.habitat]} · {rodal.altitud:.0f} m · "
           f"{len(rodal.arbres)} arbres inventariats")

# Avisos de rang de calibratge.
pins = rodal.arbres[rodal.arbres["sp"].isin(ESPECIES)]
sp_dom = pins.assign(g=pins["dbh"] ** 2 * pins["f"]).groupby("sp")["g"].sum().idxmax()
cal = calibracio["especies"][ESPECIES[sp_dom]]["calibratge"]
ab_total = area_basimetrica(rodal.arbres)
altres = 1 - area_basimetrica(pins) / ab_total
if not cal["altitud_min"] <= rodal.altitud <= cal["altitud_max"]:
    st.warning(f"L'altitud queda fora del rang de calibratge del {NOMS_COMUNS[sp_dom].lower()} "
               f"({cal['altitud_min']:.0f}–{cal['altitud_max']:.0f} m).")
if ab_total > cal["ab_max"]:
    st.warning(f"L'àrea basimètrica ({ab_total:.0f} m²/ha) supera la màxima de les parcel·les de calibratge "
               f"({cal['ab_max']:.0f} m²/ha): el creixement s'extrapola.")
if altres > 0.05:
    st.info(f"El {altres * 100:.0f} % de l'àrea basimètrica és d'altres espècies: compten per a la competència "
            "però no creixen ni s'inclouen en el carboni.")

# ---------- Indicadors ----------
co2_ara = serie["co2"].iloc[0]
captacio = serie["co2"].iloc[-1] - co2_ara
k1, k2, k3, k4 = st.columns(4)
k1.metric("CO₂ emmagatzemat avui", f"{fmt(co2_ara, 0)} t/ha", help="Arbres vius de pi: part aèria i subterrània.")
k2.metric(f"Captació neta en {anys} anys", f"{fmt(captacio, 0)} t CO₂/ha",
          f"{fmt(captacio / anys, 1)} t CO₂/ha·any", delta_color="off")
if sim_amb is not None:
    co2_extret = co2_ara - sim_amb.serie["co2"].iloc[0]
    captacio_amb = sim_amb.serie["co2"].iloc[-1] - sim_amb.serie["co2"].iloc[0]
    diferencia = captacio_amb - captacio
    k3.metric("Captació neta amb aclarida", f"{fmt(captacio_amb, 0)} t CO₂/ha",
              f"{'+' if diferencia >= 0 else '−'}{fmt(abs(diferencia), 0)} t vs. sense aclarida")
    k4.metric("CO₂ en la fusta extreta", f"{fmt(co2_extret, 0)} t/ha",
              help="Carboni dels arbres tallats. Que torni o no a l'atmosfera depèn de l'ús de la fusta.")
    balanc = sim_amb.serie["co2"].iloc[-1] - serie["co2"].iloc[-1]
    if diferencia > 0:
        st.success(
            f"L'aclarida **augmenta la captació** en {fmt(diferencia, 0)} t CO₂/ha en {anys} anys. "
            f"Tot i així, l'estoc en peu al final és {fmt(abs(balanc), 0)} t CO₂/ha "
            f"{'inferior' if balanc < 0 else 'superior'} al del rodal sense aclarir."
        )
    else:
        st.error(
            f"L'aclarida **redueix la captació** en {fmt(-diferencia, 0)} t CO₂/ha en {anys} anys. "
            f"L'estoc en peu al final és {fmt(abs(balanc), 0)} t CO₂/ha {'inferior' if balanc < 0 else 'superior'}."
        )
else:
    k3.metric("Carboni aeri avui", f"{fmt(serie['cat'].iloc[0], 0)} t C/ha")
    k4.metric("Àrea basimètrica", f"{fmt(serie['ab'].iloc[0], 1)} m²/ha",
              f"{fmt(serie['ab'].iloc[-1] - serie['ab'].iloc[0], 1)} en {anys} anys", delta_color="off")


# ---------- Gràfics ----------
def estil(fig, titol_y, alcada=420):
    fig.update_layout(template="simple_white", height=alcada, margin=dict(l=10, r=10, t=30, b=10),
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0), yaxis_title=titol_y,
                      hovermode="x unified", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    fig.update_xaxes(showgrid=True, gridcolor="#ECE8DD")
    fig.update_yaxes(showgrid=True, gridcolor="#ECE8DD")
    return fig


tabs = st.tabs(["📈 Evolució", "📊 Diàmetres", "🗺️ Mapa", "📋 Taula", "🔬 Model i validació", "📖 Com funciona"])

with tabs[0]:
    variable = st.segmented_control("Variable", list(VARIABLES), default="co2",
                                    format_func=lambda v: VARIABLES[v][0], label_visibility="collapsed") or "co2"
    nom, unitats = VARIABLES[variable]
    fig = go.Figure()
    inici = 1 if variable == "increment_brut" else 0
    fig.add_trace(go.Scatter(x=serie["any"][inici:], y=serie[variable][inici:], name="Sense aclarida",
                             line=dict(color=COLOR_SENSE, width=3)))
    if sim_amb is not None:
        fig.add_trace(go.Scatter(x=sim_amb.serie["any"][inici:], y=sim_amb.serie[variable][inici:],
                                 name="Amb aclarida", line=dict(color=COLOR_AMB, width=3)))
    fig.update_xaxes(title="Anys des d'avui")
    st.plotly_chart(estil(fig, f"{nom} ({unitats})"), width="stretch")

with tabs[1]:
    classes = np.arange(7.5, max(80.0, sim.arbres_finals["dbh"].max() + 5), 5.0)
    centres = (classes[:-1] + classes[1:]) / 2

    def histograma(a: pd.DataFrame) -> np.ndarray:
        a = a[a["sp"].isin(ESPECIES)]
        return np.histogram(a["dbh"], bins=classes, weights=a["f"])[0]

    fig = go.Figure()
    fig.add_bar(x=centres, y=histograma(sim.arbres_inicials), name="Avui", marker_color=COLOR_INICIAL)
    fig.add_bar(x=centres, y=histograma(sim.arbres_finals), name=f"D'aquí a {anys} anys, sense aclarida",
                marker_color=COLOR_SENSE)
    if sim_amb is not None:
        fig.add_bar(x=centres, y=histograma(extrets), name="Peus extrets", marker_color=COLOR_EXTRET)
        fig.add_bar(x=centres, y=histograma(sim_amb.arbres_finals), name=f"D'aquí a {anys} anys, amb aclarida",
                    marker_color=COLOR_AMB)
    fig.update_layout(barmode="group", hovermode="x")
    fig.update_xaxes(title="Classe diamètrica (cm)", dtick=5)
    st.plotly_chart(estil(fig, "Peus/ha"), width="stretch")

with tabs[2]:
    mapa = estoc_actual_pinars()
    fig = go.Figure(go.Scattermap(
        lat=mapa["lat"], lon=mapa["lon"], mode="markers",
        marker=dict(size=9, color=mapa["co2"], colorscale="Greens", cmin=0, showscale=True,
                    colorbar=dict(title="t CO₂/ha")),
        text=[etiqueta_parcela(f) for f in mapa.itertuples(index=False)], hovertemplate="%{text}<extra></extra>",
    ))
    seleccionada = mapa[mapa["parcela"] == codi_seleccionat]
    if not seleccionada.empty:
        fig.add_trace(go.Scattermap(lat=seleccionada["lat"], lon=seleccionada["lon"], mode="markers",
                                    marker=dict(size=18, color=COLOR_AMB), hoverinfo="skip"))
    fig.update_layout(map=dict(style="open-street-map", center=dict(lat=42.54, lon=1.57), zoom=9.6), height=560,
                      margin=dict(l=0, r=0, t=0, b=0), showlegend=False)
    st.plotly_chart(fig, width="stretch")
    st.caption(f"{len(mapa)} pinars del 1r INF d'Andorra (2022), acolorits pel CO₂ emmagatzemat. "
               "En taronja, la parcel·la seleccionada.")

with tabs[3]:
    columnes = ["any", "densitat", "ab", "dq", "ho", "bat", "bst", "cat", "cst", "co2", "increment_brut"]
    taula = serie[columnes].assign(escenari="Sense aclarida")
    if sim_amb is not None:
        taula = pd.concat([taula, sim_amb.serie[columnes].assign(escenari="Amb aclarida")], ignore_index=True)
    taula = taula[["escenari"] + columnes].round(2)
    st.dataframe(taula, hide_index=True, width="stretch", column_config={
        "escenari": "Escenari", "any": "Any", "densitat": "N (peus/ha)", "ab": "AB (m²/ha)", "dq": "Dq (cm)",
        "ho": "Ho (m)", "bat": "BAT (t/ha)", "bst": "BST (t/ha)", "cat": "CAT (t C/ha)", "cst": "CST (t C/ha)",
        "co2": "CO₂ (t/ha)", "increment_brut": "Creixement (t C/ha·any)",
    })
    st.download_button("Descarregar CSV", taula.to_csv(index=False).encode("utf-8"), "projeccio_carboni.csv",
                       "text/csv")

with tabs[4]:
    st.markdown(
        "Els models es van calibrar amb **15.500 pins** remesurats entre l'IFN2 (1989–90) i l'IFN3 (2000–01) a "
        "Lleida, i s'han validat **per validació creuada per parcel·les**: cada parcel·la es prediu amb models "
        "ajustats sense ella."
    )
    files = []
    for especie, p in calibracio["especies"].items():
        v = p["validacio"]
        files.append({
            "Espècie": especie,
            "Arbres / parcel·les": f"{p['calibratge']['arbres']} / {p['calibratge']['parceles']}",
            "AB a 11 anys: R²": round(v["parceles"]["r2"], 3),
            "AB a 11 anys: biaix (m²/ha)": round(v["parceles"]["biaix"], 2),
            "Increment d'AB: R²": round(v["increment_ab_parceles"]["r2"], 2),
            "Creixement diametral: RMSE (cm/any)": round(v["creixement"]["rmse"], 3),
            "Morts predits / observats": f"{v['mortalitat']['morts_predits']:.0f} / {v['mortalitat']['morts_observats']}",
            "SDI màxim": round(p["sdi_max"]),
        })
    st.dataframe(pd.DataFrame(files), hide_index=True, width="stretch")

    c1, c2 = st.columns(2)
    validacio = pd.read_csv(DATA / "validacio_parceles.csv")
    fig = go.Figure()
    for especie, color in zip(ESPECIES.values(), [COLOR_SENSE, COLOR_AMB]):
        v = validacio[validacio["especie"] == especie]
        fig.add_trace(go.Scatter(x=v["ab_2001_observada"], y=v["ab_2001_predita"], mode="markers", name=especie,
                                 marker=dict(color=color, size=6, opacity=0.6)))
    maxim = validacio[["ab_2001_observada", "ab_2001_predita"]].max().max()
    fig.add_trace(go.Scatter(x=[0, maxim], y=[0, maxim], mode="lines", line=dict(color="#999", dash="dash"),
                             showlegend=False))
    fig.update_xaxes(title="AB observada a l'IFN3 (m²/ha)")
    fig.update_layout(title=dict(text="Parcel·les de Lleida: AB predita vs observada", font=dict(size=14)),
                      hovermode="closest")
    c1.plotly_chart(estil(fig, "AB predita (m²/ha)", 400), width="stretch")

    fig = go.Figure()
    for especie, color in zip(ESPECIES.values(), [COLOR_SENSE, COLOR_AMB]):
        comp = pd.DataFrame(calibracio["especies"][especie]["validacio"]["comparacio_informe"])
        etiquetes = [f"{a:.0f}–{b:.0f}" if b < 100 else f">{a:.0f}" for a, b in zip(comp["ab_min"], comp["ab_max"])]
        fig.add_trace(go.Scatter(x=etiquetes, y=comp["model_sobre_observat"], name=f"Aquest model · {especie}",
                                 line=dict(color=color, width=3)))
        fig.add_trace(go.Scatter(x=etiquetes, y=comp["informe_sobre_observat"], name=f"Equacions IAVC · {especie}",
                                 line=dict(color=color, width=2, dash="dot")))
    fig.add_hline(y=1, line=dict(color="#999", dash="dash"))
    fig.update_xaxes(title="Àrea basimètrica de la parcel·la (m²/ha)")
    fig.update_layout(title=dict(text="Creixement en volum predit / observat", font=dict(size=14)))
    c2.plotly_chart(estil(fig, "Predit / observat", 400), width="stretch")
    st.caption(
        "A la dreta, el creixement en volum predit dividit per l'observat a les parcel·les de Lleida, per classes de "
        "densitat. Les equacions IAVC de l'IFN2 (les que fa servir l'informe de l'INF1-AND per a l'embornal) no "
        "tenen en compte la competència i sobreestimen el creixement de les masses denses; aquest model el reprodueix "
        "en totes les classes."
    )
    with st.expander("Coeficients dels models"):
        for especie, p in calibracio["especies"].items():
            st.markdown(f"**{especie}** · correcció de biaix {p['correccio_biaix']:.3f}")
            st.json({"creixement": p["creixement"], "mortalitat": p["mortalitat"]}, expanded=False)

with tabs[5]:
    st.markdown((ARREL / "docs" / "metode.md").read_text(encoding="utf-8"))
