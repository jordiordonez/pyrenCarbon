import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from carboni.creixement import (
    EPISODIS,
    MODES_PERTORBACIO,
    SENSE_PERTORBACIONS,
    TIPUS_PERTORBACIO,
    carregar_models,
    carregar_pertorbacions,
    pertorbacio_episodi,
    pertorbacio_pessimista,
    pertorbacio_risc_mitja,
    simular,
)
from carboni.equacions import ESPECIES, HABITATS, NOMS_COMUNS
from carboni.intervencio import TIPUS, aclarir, area_basimetrica
from carboni.inventari import (
    Rodal,
    carregar_inf1and,
    pinars,
    resum_parceles,
    rodal_analog,
    rodal_de_csv,
    rodal_de_parcela,
)

ARREL = Path(__file__).resolve().parent
DATA = ARREL / "data"
VERSIO = "0.2 · octubre 2026"

BOSC = "#173F35"
VERD = "#2D6A55"
BLAU = "#2E5F8A"
TINTA = "#202D28"
GRIS = "#59675F"
LINIA = "#D9E1D9"
GRIS_REFERENCIA = "#A7B3AB"
TERRA = "#8C7A62"
FONT = "Source Sans 3, Source Sans Pro, sans-serif"

SENSE = "Sense aclarida"
AMB = "Amb aclarida"
VARIABLES = {
    "co2": ("CO₂", "CO₂ emmagatzemat", "t CO₂/ha"),
    "carboni": ("Carboni", "Carboni total", "t C/ha"),
    "cat": ("C aeri", "Carboni aeri", "t C/ha"),
    "increment_brut": ("Creixement", "Creixement anual en carboni", "t C/ha·any"),
    "ab": ("AB", "Àrea basimètrica", "m²/ha"),
    "densitat": ("Densitat", "Densitat", "peus/ha"),
    "dq": ("Dq", "Diàmetre quadràtic mitjà", "cm"),
    "ho": ("Alçada", "Alçada dominant", "m"),
}
NOMS_ESPECIE = {f"{NOMS_COMUNS[c]} · {ESPECIES[c]}": c for c in ESPECIES}
NOMS_TIPUS_CURTS = {"baix": "per baix", "dalt": "per dalt", "sistematica": "sistemàtica"}

# Corbes de nivell abstractes per a la franja de capçalera (recurs gràfic, no cartografia).
TEXTURA = (
    "data:image/svg+xml;utf8,"
    "<svg xmlns='http://www.w3.org/2000/svg' width='600' height='220' viewBox='0 0 600 220'>"
    "<g fill='none' stroke='%23FFFFFF' stroke-opacity='0.07' stroke-width='1.2'>"
    "<path d='M-20 170 C 90 120, 180 200, 300 150 S 520 90, 640 140'/>"
    "<path d='M-20 140 C 100 90, 190 170, 310 120 S 520 60, 640 110'/>"
    "<path d='M-20 110 C 110 60, 200 140, 320 90 S 520 30, 640 80'/>"
    "<path d='M-20 80 C 120 30, 210 110, 330 60 S 520 0, 640 50'/>"
    "<path d='M-20 200 C 80 150, 170 230, 290 180 S 520 120, 640 170'/>"
    "<path d='M380 40 C 420 20, 470 30, 480 60 S 430 100, 400 80 S 360 55, 380 40Z'/>"
    "<path d='M395 50 C 420 38, 455 45, 460 62 S 430 85, 410 72 S 385 58, 395 50Z'/>"
    "</g></svg>"
)

st.set_page_config(page_title="pyrenCarbon · Carboni a les pinedes d'Andorra", page_icon="🌲", layout="wide")
st.markdown(
    f"""
    <style>
    .block-container {{padding-top: 3.75rem; padding-bottom: 3rem; max-width: 1400px;}}
    .capcalera {{
        background: {BOSC} url("{TEXTURA}") right center / cover no-repeat;
        border-radius: 12px; padding: 20px 32px 26px; margin-bottom: 20px; color: #FFFFFF;
    }}
    .capcalera .signatura {{
        display: flex; justify-content: space-between; gap: 16px; flex-wrap: wrap;
        font-size: 14px; color: #CFE0D6; border-bottom: 1px solid rgba(255,255,255,0.14);
        padding-bottom: 12px; margin-bottom: 18px;
    }}
    .capcalera .signatura strong {{color: #FFFFFF; font-weight: 600;}}
    .capcalera h1 {{color: #FFFFFF; font-size: 38px; line-height: 1.15; margin: 0 0 6px; padding: 0;}}
    .capcalera .lema {{color: #FFFFFF; font-size: 22px; font-weight: 600; margin: 0 0 8px;}}
    .capcalera p {{color: #D7E5DC; font-size: 17px; margin: 0; max-width: 62ch;}}
    .context {{color: {GRIS}; font-size: 15px; margin: 4px 0 14px;}}
    .context strong {{color: {TINTA}; font-weight: 600;}}
    .escenari {{display: inline-block; width: 14px; height: 4px; border-radius: 2px; vertical-align: middle;
               margin: 0 6px 2px 2px;}}
    .targeta {{
        background: #FFFFFF; border: 1px solid {LINIA}; border-left: 4px solid var(--accent);
        border-radius: 12px; padding: 16px 20px; min-height: 138px; margin-bottom: 16px;
    }}
    .targeta .pregunta {{font-size: 15px; color: {GRIS}; margin-bottom: 4px;}}
    .targeta .xifra {{font-size: 36px; font-weight: 600; color: {TINTA}; line-height: 1.15;
                      font-variant-numeric: tabular-nums;}}
    .targeta .unitat {{font-size: 16px; font-weight: 400; color: {GRIS}; margin-left: 4px;}}
    .targeta .detall {{font-size: 14px; color: {GRIS}; margin-top: 6px;}}
    .lectura {{background: #FFFFFF; border: 1px solid {LINIA}; border-radius: 12px; padding: 14px 20px;
               font-size: 16px; color: {TINTA}; margin-bottom: 16px; line-height: 1.5;}}
    .avis {{background: #FFF4DC; border: 1px solid #E8D2A6; border-radius: 10px; padding: 12px 18px;
            color: #8A5A16; margin-bottom: 4px;}}
    .avis strong {{font-weight: 600;}}
    .lectura-metode {{max-width: 74ch; line-height: 1.55; margin-bottom: 16px;}}
    .peu {{border-top: 1px solid {LINIA}; margin-top: 40px; padding-top: 16px; font-size: 14px; color: {GRIS};
           line-height: 1.6;}}
    section[data-testid="stSidebar"] h3 {{font-size: 18px; margin-top: 8px;}}
    @media (max-width: 640px) {{
        .capcalera {{padding: 16px 18px 20px;}}
        .capcalera h1 {{font-size: 29px;}}
        .targeta .xifra {{font-size: 30px;}}
    }}
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
    pertorbacions = carregar_pertorbacions(DATA / "calibracio.json")
    return resum, arbres, models, calibracio, parceles.set_index("parcela")["url_img"], pertorbacions


@st.cache_data
def estoc_actual_pinars() -> pd.DataFrame:
    resum, arbres, models, _, _, _ = carregar()
    seleccio = pinars(resum).copy()
    seleccio["co2"] = [
        simular(rodal_de_parcela(resum, arbres, c), models, 0).serie["co2"].iloc[0] for c in seleccio["parcela"]
    ]
    return seleccio


@st.cache_data(show_spinner="Actualitzant resultats…", max_entries=500)
def simular_cache(arbres: pd.DataFrame, habitat: str, altitud: float, descripcio: str, anys: int,
                  mode_pert: str, tipus_pert: str, any_pert: int, severitat_pert: float):
    _, _, models, _, _, parametres = carregar()
    if mode_pert == "risc_mitja":
        pertorbacio = pertorbacio_risc_mitja(parametres)
    elif mode_pert == "pessimista":
        pertorbacio = pertorbacio_pessimista(parametres, anys)
    elif mode_pert == "episodi":
        pertorbacio = pertorbacio_episodi(tipus_pert, any_pert, severitat_pert)
    else:
        pertorbacio = SENSE_PERTORBACIONS
    return simular(Rodal(arbres, habitat, altitud, descripcio), models, anys, pertorbacio)


resum, arbres_inf, models, calibracio, fotos, parametres_pert = carregar()


def especie_dominant(arbres: pd.DataFrame) -> str:
    pins = arbres[arbres["sp"].isin(ESPECIES)]
    return pins.assign(g=pins["dbh"] ** 2 * pins["f"]).groupby("sp")["g"].sum().idxmax()


def fmt(valor: float, decimals: int = 0) -> str:
    text = f"{abs(valor):,.{decimals}f}".replace(",", " ").replace(".", ",")
    return ("−" if valor < 0 and round(abs(valor), decimals) > 0 else "") + text


def fmt_signe(valor: float, decimals: int = 0) -> str:
    return ("+" if valor > 0 and round(valor, decimals) > 0 else "") + fmt(valor, decimals)


def etiqueta_parcela(fila) -> str:
    return (f"{fila.parcela} · {NOMS_COMUNS[fila.especie_dominant]} · {fila.habitat} · {fila.altitud:.0f} m · "
            f"{fila.n:.0f} peus/ha · {fila.ab:.0f} m²/ha".replace(".", ","))


def etiqueta_parcela_curta(fila) -> str:
    return f"{fila.parcela} · {NOMS_COMUNS[fila.especie_dominant]} · {fmt(fila.altitud)} m"


def targeta(pregunta: str, xifra: str, unitat: str, detall: str, accent: str) -> str:
    return (
        f'<div class="targeta" style="--accent:{accent}">'
        f'<div class="pregunta">{html.escape(pregunta, quote=False)}</div>'
        f'<div class="xifra">{xifra}<span class="unitat">{html.escape(unitat, quote=False)}</span></div>'
        f'<div class="detall">{detall}</div></div>'
    )


def marca(color: str) -> str:
    return f'<span class="escenari" style="background:{color}"></span>'


def transparent(color: str, opacitat: float) -> str:
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{opacitat})"


NOTA_BANDA = ("La franja va de l'escenari sense pertorbacions (a dalt) al cas pessimista (a baix): la pèrdua per "
              "pertorbacions en l'horitzó triat que només supera 1 de cada 10 rodals (percentil 90 d'una simulació "
              "Monte Carlo dels episodis observats al Pirineu de Lleida), repartida igual cada any.")


# ---------- Capçalera ----------
st.markdown(
    """
    <div class="capcalera">
      <div class="signatura">
        <span><strong>Eina d'estimació de carboni forestal</strong></span>
        <span>Dades del 1r Inventari Nacional Forestal d'Andorra · Andorra Recerca + Innovació</span>
      </div>
      <h1>pyrenCarbon</h1>
      <div class="lema">Carboni a les pinedes d'Andorra</div>
      <p>Explora l'estoc i la captació estimada de CO₂ de les pinedes de pi roig i pi negre,
      amb i sense una aclarida.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------- Configuració ----------
with st.sidebar:
    st.subheader("Rodal de partida")
    origen = st.radio(
        "Origen de les dades",
        ["Parcel·la de l'INF d'Andorra", "Dades del rodal", "Inventari propi (CSV)"],
        help="Una parcel·la real del 1r INF d'Andorra (2022), una densitat i àrea basimètrica, o el vostre inventari.",
    )
    rodal = None
    error = None
    codi_seleccionat = None
    if origen == "Parcel·la de l'INF d'Andorra":
        filtre = st.segmented_control("Espècie dominant", ["Totes", "Pi negre", "Pi roig"], default="Totes") or "Totes"
        sp_filtre = {"Pi negre": "Pu", "Pi roig": "Ps"}.get(filtre)
        opcions = pinars(resum, sp_filtre).sort_values("parcela")
        fila = st.selectbox(
            "Parcel·la", list(opcions.itertuples(index=False)), format_func=etiqueta_parcela_curta,
            help="Pinedes amb almenys un 70 % de l'àrea basimètrica de pi roig o pi negre.",
        )
        codi_seleccionat = fila.parcela
        rodal = rodal_de_parcela(resum, arbres_inf, codi_seleccionat)
    elif origen == "Dades del rodal":
        sp = NOMS_ESPECIE[st.selectbox("Espècie", list(NOMS_ESPECIE))]
        habitat = st.radio("Hàbitat", list(HABITATS), format_func=lambda h: HABITATS[h].capitalize(), horizontal=True)
        altitud = st.number_input("Altitud (m)", 900, 2500, 1900 if sp == "Pu" else 1500, step=50)
        c1, c2 = st.columns(2)
        n = c1.number_input("Densitat (peus/ha)", 50, 5000, 900, step=50)
        ab = c2.number_input("Àrea basimètrica (m²/ha)", 2.0, 90.0, 35.0, step=1.0)
        st.caption("La distribució de diàmetres es pren de la parcel·la andorrana més semblant.")
        rodal = rodal_analog(resum, arbres_inf, sp, habitat, float(altitud), float(n), float(ab))
    else:
        plantilla = pd.DataFrame({"especie": ["Pu", "Pu", "Ps"], "dbh": [18.5, 32.0, 24.2], "ht": [11.0, None, 14.5]})
        fitxer = st.file_uploader("Inventari (CSV)", type="csv",
                                  help="Columnes: especie (Ps, Pu o el nom), dbh en cm; opcionals: ht en m, peus_ha.")
        st.download_button("Descarregar la plantilla", plantilla.to_csv(index=False).encode(),
                           "plantilla_inventari.csv", "text/csv", type="tertiary")
        habitat = st.radio("Hàbitat", list(HABITATS), format_func=lambda h: HABITATS[h].capitalize(), horizontal=True)
        altitud = st.number_input("Altitud (m)", 900, 2500, 1900, step=50)
        superficie = st.number_input("Superfície de la parcel·la (m²)", 50.0, 5000.0, 452.0, step=10.0,
                                     help="Només si el fitxer no té la columna peus_ha. 452 m² = radi de 12 m.")
        if fitxer is not None:
            try:
                rodal = rodal_de_csv(pd.read_csv(fitxer, sep=None, engine="python"), habitat, float(altitud),
                                     superficie)
            except ValueError as e:
                error = str(e)

    st.divider()
    st.subheader("Intervenció")
    amb_intervencio = st.toggle("Comparar amb una aclarida", value=True)
    tipus = st.selectbox("Tipus d'aclarida", list(TIPUS), format_func=lambda t: NOMS_TIPUS_CURTS[t].capitalize(),
                         disabled=not amb_intervencio, help="  \n".join(TIPUS.values()))
    pes = st.slider("Àrea basimètrica extreta (%)", 5, 60, 30, step=1, disabled=not amb_intervencio)

    st.divider()
    st.subheader("Horitzó")
    anys = st.slider("Anys de projecció", 5, 80, 20, step=1, format="%d anys")

    st.divider()
    st.subheader("Pertorbacions")
    mode_pert = st.radio(
        "Escenari de pertorbacions", list(MODES_PERTORBACIO), index=1, format_func=MODES_PERTORBACIO.get,
        help="Foc, sequera, neu i vent, plagues… La mortalitat de fons (competència i causa desconeguda) "
             "s'aplica sempre.",
    )
    tipus_pert, any_pert, severitat_pert = "", 0, 0.0
    taxes = "; ".join(f"{NOMS_COMUNS[sp].lower()} {fmt(p['taxa_anual'] * 100, 2)} %"
                      for sp, p in parametres_pert.items())
    if mode_pert == "cap":
        st.caption("Només la mortalitat de fons (competència i causa desconeguda). Projecció optimista: "
                   "no hi ha cap foc, sequera, nevada, ventada ni plaga.")
    elif mode_pert == "risc_mitja":
        st.caption(f"A més de la mortalitat de fons, cada any mor la fracció de pins que es va perdre de mitjana per "
                   f"pertorbacions al Pirineu de Lleida entre 1989 i 2016 ({taxes} de l'AB). Al pi "
                   "negre, els pins gruixuts hi tenen més risc; al pi roig, totes les mides igual. És una mitjana: un rodal concret pot no patir-ne cap o patir-ne una de forta. "
                   "Els resultats mostren també l'interval entre l'escenari sense pertorbacions i un cas pessimista "
                   "(1 de cada 10 rodals).")
    else:
        tipus_pert = st.selectbox("Tipus d'episodi", list(EPISODIS), format_func=TIPUS_PERTORBACIO.get)
        any_pert = st.slider("Any de l'episodi", 1, anys, min(10, anys), step=1)
        sp_pert = especie_dominant(rodal.arbres) if rodal is not None else "Pu"
        dades_tipus = parametres_pert[sp_pert]["tipus"][tipus_pert]
        severitat_pert = st.slider(
            "Pins morts en l'episodi (% de l'AB)", 1, 100, max(1, round(dades_tipus["severitat"] * 100)), step=1,
            key=f"severitat_{tipus_pert}_{sp_pert}",
            help="Per defecte, la mitjana a les parcel·les de Lleida on aquesta causa va matar almenys el 10 % "
                 "de l'AB de pins.") / 100.0
        origen_sev = ("totes dues espècies" if dades_tipus["severitat_conjunta"]
                      else NOMS_COMUNS[sp_pert].lower())
        st.caption(f"L'any {any_pert} moren el {fmt(severitat_pert * 100)} % dels pins, de totes les mides. Valor per "
                   f"defecte: mitjana de {dades_tipus['parceles_afectades']} parcel·les de Lleida afectades "
                   f"({origen_sev}). La resta d'anys, només la mortalitat de fons.")

if rodal is None:
    if error:
        st.error(error)
    else:
        st.info("Puja un inventari en CSV per començar. Pots descarregar la plantilla a la barra lateral.")
    st.stop()

# ---------- Simulació ----------
parametres_simulacio = (mode_pert, tipus_pert, any_pert, severitat_pert)
sim = simular_cache(rodal.arbres, rodal.habitat, rodal.altitud, rodal.descripcio, anys, *parametres_simulacio)
serie = sim.serie
sim_amb, extrets = None, None
if amb_intervencio:
    rodal_aclarit, extrets = aclarir(rodal, tipus, pes / 100.0)
    sim_amb = simular_cache(rodal_aclarit.arbres, rodal_aclarit.habitat, rodal_aclarit.altitud,
                            rodal_aclarit.descripcio, anys, *parametres_simulacio)

# Banda d'incertesa del risc mitjà: entre «sense pertorbacions» i el cas pessimista (percentil 90).
banda, banda_amb = None, None
if mode_pert == "risc_mitja":
    def simular_banda(r: Rodal) -> tuple[pd.DataFrame, pd.DataFrame]:
        return tuple(simular_cache(r.arbres, r.habitat, r.altitud, r.descripcio, anys, m, "", 0, 0.0).serie
                     for m in ("cap", "pessimista"))

    banda = simular_banda(rodal)
    if sim_amb is not None:
        banda_amb = simular_banda(rodal_aclarit)


def interval_captacio(b: tuple[pd.DataFrame, pd.DataFrame] | None) -> str:
    if b is None:
        return ""
    alta, baixa = (s["co2"].iloc[-1] - s["co2"].iloc[0] for s in b)
    return f"<br>Interval: {fmt_signe(baixa)} a {fmt_signe(alta)} t CO₂/ha (cas pessimista – sense pertorbacions)"

pins = rodal.arbres[rodal.arbres["sp"].isin(ESPECIES)]
sp_dom = especie_dominant(rodal.arbres)
if mode_pert == "episodi":
    nom_pert = (f"episodi de {TIPUS_PERTORBACIO[tipus_pert].lower()} l'any {any_pert} "
                f"({fmt(severitat_pert * 100)} % de l'AB)")
else:
    nom_pert = MODES_PERTORBACIO[mode_pert].lower()
nom_aclarida = f"aclarida {NOMS_TIPUS_CURTS[tipus]} del {pes} % de l'AB"
escenaris = f"{marca(VERD)}{SENSE}" + (f" &nbsp;·&nbsp; {marca(BLAU)}{AMB.lower()}: {nom_aclarida}" if sim_amb else "")
st.markdown(
    f'<div class="context"><strong>{NOMS_COMUNS[sp_dom]}</strong> · {html.escape(rodal.descripcio, quote=False)} · '
    f"hàbitat {HABITATS[rodal.habitat]} · horitzó <strong>{anys} anys</strong> · pertorbacions: {nom_pert}"
    f"<br>{escenaris}</div>",
    unsafe_allow_html=True,
)

# ---------- Avisos ----------
cal = calibracio["especies"][ESPECIES[sp_dom]]["calibratge"]
ab_total = area_basimetrica(rodal.arbres)
altres = 1 - area_basimetrica(pins) / ab_total
avisos = []
if not cal["altitud_min"] <= rodal.altitud <= cal["altitud_max"]:
    avisos.append(f"L'altitud ({fmt(rodal.altitud)} m) queda fora del rang de calibratge del "
                  f"{NOMS_COMUNS[sp_dom].lower()} ({fmt(cal['altitud_min'])}–{fmt(cal['altitud_max'])} m).")
if ab_total > cal["ab_max"]:
    avisos.append(f"L'àrea basimètrica del rodal ({fmt(ab_total)} m²/ha) supera la màxima de les parcel·les de "
                  f"calibratge ({fmt(cal['ab_max'])} m²/ha): el creixement s'extrapola.")
if altres > 0.05:
    avisos.append(f"El {fmt(altres * 100)} % de l'àrea basimètrica és d'altres espècies: compten per a la "
                  "competència, però no creixen ni s'inclouen en el carboni.")
if avisos:
    titol = ("Projecció parcialment extrapolada" if any("extrapola" in a or "fora del rang" in a for a in avisos)
             else "Cal tenir en compte")
    st.markdown(f'<div class="avis"><strong>{titol}.</strong> Afecta tots dos escenaris; '
                f'{len(avisos)} {"avís" if len(avisos) == 1 else "avisos"}.</div>', unsafe_allow_html=True)
    with st.expander("Veure el detall"):
        for avis in avisos:
            st.markdown(f"- {avis}")

# ---------- Resultats ----------
co2_ara = serie["co2"].iloc[0]
captacio = serie["co2"].iloc[-1] - co2_ara
fila1 = st.columns(2, gap="medium")
fila1[0].markdown(targeta(
    "CO₂ emmagatzemat avui", fmt(co2_ara), "t CO₂/ha",
    f"Pins vius, part aèria i arrels · carboni aeri {fmt(serie['cat'].iloc[0])} t C/ha", GRIS_REFERENCIA,
), unsafe_allow_html=True)
fila1[1].markdown(targeta(
    f"Captació estimada en {anys} anys · sense aclarida", fmt_signe(captacio), "t CO₂/ha",
    f"{marca(VERD)}{fmt(captacio / anys, 1)} t CO₂/ha·any de mitjana · estoc final {fmt(serie['co2'].iloc[-1])} t CO₂/ha"
    + interval_captacio(banda), VERD,
), unsafe_allow_html=True)

if sim_amb is not None:
    co2_extret = co2_ara - sim_amb.serie["co2"].iloc[0]
    captacio_amb = sim_amb.serie["co2"].iloc[-1] - sim_amb.serie["co2"].iloc[0]
    diferencia = captacio_amb - captacio
    balanc = sim_amb.serie["co2"].iloc[-1] - serie["co2"].iloc[-1]
    fila2 = st.columns(2, gap="medium")
    fila2[0].markdown(targeta(
        f"Captació estimada en {anys} anys · amb aclarida", fmt_signe(captacio_amb), "t CO₂/ha",
        f"{marca(BLAU)}{fmt_signe(diferencia)} t CO₂/ha respecte de l'escenari sense aclarida"
        + interval_captacio(banda_amb), BLAU,
    ), unsafe_allow_html=True)
    fila2[1].markdown(targeta(
        "CO₂ retirat en l'aclarida", fmt(co2_extret), "t CO₂/ha",
        "Carboni dels arbres tallats. Que torni o no a l'atmosfera depèn de l'ús de la fusta.", TERRA,
    ), unsafe_allow_html=True)
    sentit = "més" if diferencia > 0 else "menys"
    st.markdown(
        f'<div class="lectura">En {anys} anys, el rodal aclarit capta <strong>{fmt(abs(diferencia))} t CO₂/ha '
        f"{sentit}</strong> que el rodal sense aclarir. Si se suma el CO₂ retirat, l'estoc en peu final del rodal "
        f"aclarit és <strong>{fmt(abs(balanc))} t CO₂/ha {'inferior' if balanc < 0 else 'superior'}</strong> al "
        "del rodal sense aclarir.</div>",
        unsafe_allow_html=True,
    )
else:
    fila2 = st.columns(2, gap="medium")
    fila2[0].markdown(targeta(
        "Àrea basimètrica avui", fmt(serie["ab"].iloc[0], 1), "m²/ha",
        f"{fmt(serie['densitat'].iloc[0])} peus/ha · diàmetre quadràtic mitjà {fmt(serie['dq'].iloc[0], 1)} cm",
        GRIS_REFERENCIA,
    ), unsafe_allow_html=True)
    fila2[1].markdown(targeta(
        f"Àrea basimètrica d'aquí a {anys} anys", fmt(serie["ab"].iloc[-1], 1), "m²/ha",
        f"{marca(VERD)}{fmt(serie['densitat'].iloc[-1])} peus/ha · mortalitat"
        f"{'' if mode_pert == 'cap' else ', pertorbacions'} i autoaclarida incloses", VERD,
    ), unsafe_allow_html=True)

if codi_seleccionat is not None and isinstance(fotos.get(codi_seleccionat), str):
    with st.expander(f"Fotografies de camp de la parcel·la {codi_seleccionat}"):
        st.image(fotos[codi_seleccionat], width="stretch")
        st.caption("Fotografies cap al nord, l'est, el sud i l'oest des del centre de la parcel·la. "
                   "Font: 1r Inventari Nacional Forestal d'Andorra, Andorra Recerca + Innovació (2022).")


# ---------- Gràfics ----------
def estil(fig: go.Figure, titol: str, subtitol: str, titol_y: str, alcada: int = 500) -> go.Figure:
    fig.update_layout(
        template="simple_white", height=alcada, separators=", ",
        font=dict(family=FONT, size=14, color=TINTA),
        title=dict(text=titol, subtitle=dict(text=subtitol, font=dict(size=14, color=GRIS)),
                   x=0, xanchor="left", y=0.97, yanchor="top", font=dict(size=18, weight=600)),
        legend=dict(orientation="h", x=0, xanchor="left", y=-0.2, yanchor="top", font=dict(size=14),
                    entrywidth=250, entrywidthmode="pixels"),
        margin=dict(l=8, r=16, t=84, b=8), yaxis_title=titol_y,
        paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF", hoverlabel=dict(font=dict(family=FONT, size=14)),
    )
    fig.update_xaxes(showgrid=False, linecolor=GRIS_REFERENCIA, tickfont=dict(size=13))
    fig.update_yaxes(showgrid=True, gridcolor="#EEF1EC", linecolor=GRIS_REFERENCIA, tickfont=dict(size=13),
                     rangemode="tozero")
    return fig


CONFIG_GRAFIC = {"displaylogo": False, "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"],
                 "locale": "ca"}

tabs = st.tabs(["Evolució", "Diàmetres", "Mapa", "Taula", "Model i validació", "Com funciona"])

with tabs[0]:
    variable = st.segmented_control("Variable", list(VARIABLES), default="co2",
                                    format_func=lambda v: VARIABLES[v][0], label_visibility="collapsed") or "co2"
    _, nom, unitats = VARIABLES[variable]
    inici = 1 if variable == "increment_brut" else 0
    fig = go.Figure()
    traces = [(serie, SENSE, VERD, "circle", banda)]
    if sim_amb:
        traces.append((sim_amb.serie, AMB, BLAU, "square", banda_amb))
    for s, etiqueta, color, simbol, b in traces:
        if b is not None:
            for k, (sb, extrem) in enumerate(zip(b, ("sense pertorbacions", "cas pessimista"))):
                fig.add_trace(go.Scatter(
                    x=sb["any"][inici:], y=sb[variable][inici:], mode="lines", line=dict(width=0, color=color),
                    fill="tonexty" if k else None, fillcolor=transparent(color, 0.14), legendgroup=etiqueta,
                    name=f"{etiqueta}: interval", showlegend=bool(k),
                    hovertemplate=f"{etiqueta}, {extrem}: %{{y:,.1f}} {unitats}<extra></extra>",
                ))
        fig.add_trace(go.Scatter(
            x=s["any"][inici:], y=s[variable][inici:], name=etiqueta, mode="lines+markers", legendgroup=etiqueta,
            line=dict(color=color, width=3), marker=dict(symbol=simbol, size=6, color=color),
            hovertemplate=f"{etiqueta}: %{{y:,.1f}} {unitats}<extra></extra>",
        ))
    fig.update_xaxes(title="Anys des d'avui")
    fig.update_layout(hovermode="x unified")
    escenaris_text = SENSE + (f" i {AMB.lower()} ({nom_aclarida})" if sim_amb else "")
    with st.container(border=True):
        st.plotly_chart(estil(fig, f"Evolució estimada · {nom}",
                              f"{escenaris_text} · pertorbacions: {nom_pert} · horitzó {anys} anys",
                              unitats), width="stretch", config=CONFIG_GRAFIC)
        if banda is not None:
            st.caption(NOTA_BANDA)

with tabs[1]:
    classes = np.arange(7.5, max(80.0, sim.arbres_finals["dbh"].max() + 5), 5.0)
    centres = (classes[:-1] + classes[1:]) / 2

    def histograma(a: pd.DataFrame) -> np.ndarray:
        a = a[a["sp"].isin(ESPECIES)]
        return np.histogram(a["dbh"], bins=classes, weights=a["f"])[0]

    fig = go.Figure()
    fig.add_bar(x=centres, y=histograma(sim.arbres_inicials), name="Avui", marker_color=GRIS_REFERENCIA)
    fig.add_bar(x=centres, y=histograma(sim.arbres_finals), name=f"{SENSE}, d'aquí a {anys} anys", marker_color=VERD)
    if sim_amb is not None:
        fig.add_bar(x=centres, y=histograma(extrets), name="Peus tallats", marker_color=TERRA)
        fig.add_bar(x=centres, y=histograma(sim_amb.arbres_finals), name=f"{AMB}, d'aquí a {anys} anys",
                    marker_color=BLAU)
    fig.update_layout(barmode="group", hovermode="x unified")
    fig.update_traces(hovertemplate="%{y:,.0f} peus/ha")
    fig.update_xaxes(title="Classe diamètrica (cm)", dtick=5)
    with st.container(border=True):
        st.plotly_chart(estil(fig, "Distribució de diàmetres", "Pins vius per classes de 5 cm", "peus/ha"),
                        width="stretch", config=CONFIG_GRAFIC)

with tabs[2]:
    mapa = estoc_actual_pinars()
    fig = go.Figure(go.Scattermap(
        lat=mapa["lat"], lon=mapa["lon"], mode="markers",
        marker=dict(size=10, color=mapa["co2"], colorscale=[[0, "#E3EEE6"], [0.5, "#5E9B7E"], [1, BOSC]], cmin=0,
                    showscale=True, colorbar=dict(title=dict(text="t CO₂/ha", font=dict(size=13)), thickness=12,
                                                      len=0.75, y=0.42)),
        text=[etiqueta_parcela(f) for f in mapa.itertuples(index=False)],
        customdata=mapa["co2"], hovertemplate="%{text}<br>%{customdata:,.0f} t CO₂/ha<extra></extra>",
    ))
    seleccionada = mapa[mapa["parcela"] == codi_seleccionat]
    if not seleccionada.empty:
        fig.add_trace(go.Scattermap(lat=seleccionada["lat"], lon=seleccionada["lon"], mode="markers",
                                    marker=dict(size=20, color=BLAU), hoverinfo="skip"))
    fig.update_layout(map=dict(style="carto-positron", center=dict(lat=42.545, lon=1.575), zoom=10.7), height=560,
                      margin=dict(l=0, r=0, t=0, b=0), showlegend=False, separators=", ",
                      font=dict(family=FONT))
    with st.container(border=True):
        st.plotly_chart(fig, width="stretch", config={"displaylogo": False, "locale": "ca"})
    st.caption(f"{len(mapa)} pinedes del 1r INF d'Andorra (2022), acolorits pel CO₂ emmagatzemat"
               + (". En blau, la parcel·la seleccionada." if not seleccionada.empty else "."))

with tabs[3]:
    columnes = {
        "any": "Any", "densitat": "Densitat (peus/ha)", "ab": "AB (m²/ha)", "dq": "Dq (cm)", "ho": "Ho (m)",
        "bat": "Biomassa aèria (t/ha)", "bst": "Biomassa arrels (t/ha)", "cat": "C aeri (t/ha)",
        "cst": "C arrels (t/ha)", "co2": "CO₂ (t/ha)", "increment_brut": "Creixement (t C/ha·any)",
        "morts_peus": "Morts (peus/ha·any)", "morts_pertorbacio": "Morts per pertorbacions (peus/ha·any)",
    }
    taula = serie[list(columnes)].assign(escenari=SENSE)
    if sim_amb is not None:
        taula = pd.concat([taula, sim_amb.serie[list(columnes)].assign(escenari=AMB)], ignore_index=True)
    taula = taula[["escenari"] + list(columnes)].rename(columns={"escenari": "Escenari"} | columnes)
    formats = {c: "{:,.2f}" for c in taula.columns if c not in ("Escenari", "Any", "Densitat (peus/ha)")}
    formats["Densitat (peus/ha)"] = "{:,.0f}"
    formats["Morts (peus/ha·any)"] = formats["Morts per pertorbacions (peus/ha·any)"] = "{:,.1f}"
    estil_taula = taula.style.format(formats, decimal=",", thousands=" ").hide(axis="index")
    with st.container(border=True):
        st.dataframe(estil_taula, hide_index=True, width="stretch", height=420)
    c1, c2 = st.columns([3, 1])
    with c1.expander("Glossari"):
        st.markdown(
            "- **AB**: àrea basimètrica, suma de les seccions dels troncs a 1,30 m, per hectàrea.\n"
            "- **Dq**: diàmetre quadràtic mitjà, el de l'arbre d'àrea basimètrica mitjana.\n"
            "- **Ho**: alçada dominant, mitjana del 20 % de pins més gruixuts.\n"
            "- **C**: carboni, el 50 % de la biomassa. **CO₂** = C × 44/12 (aeri i arrels).\n"
            "- **Creixement**: carboni que fixa el creixement dels arbres en un any, abans de la mortalitat.\n"
            "- **Morts**: peus morts en l'any per mortalitat de fons, pertorbacions i autoaclarida.\n"
            "- **Pertorbacions**: foc, sequera, neu i vent, plagues i altres causes externes de mort."
        )
    c2.download_button("Descarregar CSV", taula.to_csv(index=False).encode("utf-8"), "projeccio_carboni.csv",
                       "text/csv", type="secondary", width="stretch")

with tabs[4]:
    pins_calibratge = sum(p["calibratge"]["arbres"] for p in calibracio["especies"].values())
    st.markdown(
        '<div class="lectura-metode">Els models de creixement i mortalitat es van calibrar amb <strong>'
        f"{pins_calibratge:,} registres de pins".replace(",", " ") + "</strong> remesurats a les comarques "
        "pirinenques de Lleida (Pallars Sobirà, Alt Urgell, Alta Ribagorça i Cerdanya) entre l'IFN2 (1989–90) i "
        "l'IFN3 (2000–01) i entre l'IFN3 i l'IFN4 (2014–16), amb un terme de període; les prediccions fan servir "
        "el període més recent. La mortalitat del model és la <strong>de fons</strong> (competència i causa "
        "desconeguda): els pins morts per foc, sequera, neu, vent, plagues o altres pertorbacions, segons l'agent "
        "causant que anota l'IFN, es tracten a part com a escenaris. S'han validat amb <strong>validació creuada "
        "per parcel·les</strong>: cada parcel·la es prediu amb models ajustats sense ella."
        "</div>",
        unsafe_allow_html=True,
    )
    files = []
    for especie, p in calibracio["especies"].items():
        v = p["validacio"]
        files.append({
            "Espècie": especie,
            "Arbres": p["calibratge"]["arbres"],
            "Parcel·les": p["calibratge"]["parceles"],
            "R² de l'AB final": v["parceles"]["r2"],
            "Biaix de l'AB (m²/ha)": v["parceles"]["biaix"],
            "R² de l'increment d'AB": v["increment_ab_parceles"]["r2"],
            "R² sense pertorbacions": v["sense_pertorbacions"]["increment_ab_parceles"]["r2"],
            "AUC mortalitat de fons": v["mortalitat"]["auc"],
            "Morts de fons predits": v["mortalitat"]["morts_predits"],
            "Morts de fons observats": v["mortalitat"]["morts_observats"],
        })
    resum_validacio = pd.DataFrame(files).style.format(
        {"Arbres": "{:,.0f}", "Parcel·les": "{:,.0f}", "R² de l'AB final": "{:.2f}",
         "Biaix de l'AB (m²/ha)": "{:.2f}", "R² de l'increment d'AB": "{:.2f}", "R² sense pertorbacions": "{:.2f}",
         "AUC mortalitat de fons": "{:.2f}", "Morts de fons predits": "{:,.0f}",
         "Morts de fons observats": "{:,.0f}"}, decimal=",", thousands=" ").hide(axis="index")
    st.dataframe(resum_validacio, hide_index=True, width="stretch")
    st.caption("AB i increment d'AB: parcel·les de pineda, amb la mortalitat de fons i la pèrdua mitjana per "
               "pertorbacions (escenari de risc mitjà). «Sense pertorbacions»: només les parcel·les que hi van "
               "perdre menys de 2 m²/ha, amb la mortalitat de fons.")

    files = []
    for especie, p in calibracio["especies"].items():
        for tipus, c in p["pertorbacions"]["causes"].items():
            files.append({
                "Espècie": especie, "Pertorbació": TIPUS_PERTORBACIO[tipus],
                "Pèrdua anual d'AB (%)": c["taxa_perdua_ab"] * 100,
                "Parcel·les afectades": c["parceles_afectades"],
                "Probabilitat anual d'episodi (%)": c["p_anual"] * 100,
                "AB morta en l'episodi (%, mitjana)": (c["fraccio_mitjana"] or np.nan) * 100,
            })
    st.markdown("**Pertorbacions al Pirineu de Lleida (IFN2–IFN4)**")
    st.dataframe(pd.DataFrame(files).style.format(
        {"Pèrdua anual d'AB (%)": "{:.3f}", "Probabilitat anual d'episodi (%)": "{:.3f}",
         "AB morta en l'episodi (%, mitjana)": "{:.0f}"}, decimal=",", thousands=" ", na_rep="–"
    ).hide(axis="index"), hide_index=True, width="stretch")
    st.caption("Parcel·la afectada: la causa hi va matar almenys el 10 % de l'AB inicial de pins de l'espècie. "
               "Totes les parcel·les remesurades amb l'espècie, tots dos períodes junts.")

    validacio = pd.read_csv(DATA / "validacio_parceles.csv")
    fig = go.Figure()
    for especie, color, simbol in zip(ESPECIES.values(), [VERD, BLAU], ["circle", "square"]):
        v = validacio[validacio["especie"] == especie]
        fig.add_trace(go.Scatter(x=v["ab_2001_observada"], y=v["ab_2001_predita"], mode="markers", name=especie,
                                 marker=dict(color=color, size=7, opacity=0.55, symbol=simbol),
                                 hovertemplate="Observada %{x:.1f} · predita %{y:.1f} m²/ha<extra></extra>"))
    maxim = validacio[["ab_2001_observada", "ab_2001_predita"]].max().max()
    fig.add_trace(go.Scatter(x=[0, maxim], y=[0, maxim], mode="lines", name="Predicció perfecta",
                             line=dict(color=GRIS_REFERENCIA, dash="dash")))
    fig.update_xaxes(title="AB observada a la remesura (m²/ha)")
    with st.container(border=True):
        st.plotly_chart(estil(fig, "Àrea basimètrica predita i observada", "Parcel·les de pineda del Pirineu de Lleida, IFN3 i IFN4",
                              "AB predita (m²/ha)"), width="stretch", config=CONFIG_GRAFIC)

    fig = go.Figure()
    for especie, color, simbol in zip(ESPECIES.values(), [VERD, BLAU], ["circle", "square"]):
        comp = pd.DataFrame(calibracio["especies"][especie]["validacio"]["comparacio_informe"])
        etiquetes = [f"{a:.0f}–{b:.0f}" if b < 100 else f"més de {a:.0f}" for a, b in zip(comp["ab_min"], comp["ab_max"])]
        fig.add_trace(go.Scatter(x=etiquetes, y=comp["model_sobre_observat"], name=f"Aquest model · {especie}",
                                 mode="lines+markers", line=dict(color=color, width=3), marker=dict(symbol=simbol)))
        fig.add_trace(go.Scatter(x=etiquetes, y=comp["informe_sobre_observat"], name=f"Equacions IAVC · {especie}",
                                 mode="lines+markers", line=dict(color=color, width=2, dash="dot"),
                                 marker=dict(symbol=simbol, size=7)))
    fig.add_hline(y=1, line=dict(color=GRIS_REFERENCIA, dash="dash"))
    fig.update_xaxes(title="Àrea basimètrica de la parcel·la (m²/ha)")
    fig.update_traces(hovertemplate="%{y:.2f}<extra>%{fullData.name}</extra>")
    with st.container(border=True):
        st.plotly_chart(estil(fig, "Creixement en volum: predit / observat",
                              "Valor 1 = predicció exacta · parcel·les de Lleida per classes de densitat", "Ràtio"),
                        width="stretch", config=CONFIG_GRAFIC)
    st.markdown(
        '<div class="lectura-metode">Les equacions IAVC de l\'IFN2, que fa servir l\'informe de l\'INF1-AND per '
        "estimar l'embornal, no tenen en compte la competència i sobreestimen el creixement de les masses denses. "
        "Aquest model el reprodueix en totes les classes de densitat.</div>",
        unsafe_allow_html=True,
    )
    with st.expander("Coeficients dels models"):
        for especie, p in calibracio["especies"].items():
            st.markdown(f"**{especie}** · correcció de biaix {fmt(p['correccio_biaix'], 3)} · "
                        f"SDI màxim {fmt(p['sdi_max'])}")
            st.json({"creixement": p["creixement"], "mortalitat": p["mortalitat"]}, expanded=False)

with tabs[5]:
    c1, _ = st.columns([5, 2])
    with c1:
        st.markdown((ARREL / "docs" / "metode.md").read_text(encoding="utf-8"))

st.markdown(
    f"""
    <div class="peu">
      pyrenCarbon · versió {VERSIO} · <a href="https://github.com/jordiordonez/pyrenCarbon">codi font</a> · Estimacions orientatives; vegeu les limitacions a <em>Com funciona</em>.<br>
      Dades: 1r Inventari Nacional Forestal d'Andorra, BD-INF1-AND (Andorra Recerca + Innovació i CTFC, CC BY-SA 4.0);
      Segon, Tercer i Quart Inventari Forestal Nacional, província de Lleida (MITECO).
    </div>
    """,
    unsafe_allow_html=True,
)
