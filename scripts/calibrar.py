"""Calibra els models de creixement diametral i mortalitat i els valida.

Ús:
    python scripts/calibrar.py

Llegeix data/ifn_lleida_arbres.csv (scripts/preparar_ifn.py) i data/inf1and/ (INF1-AND) i escriu:
    data/calibracio.json          coeficients, límits d'autoaclarida i mètriques de validació
    data/validacio_arbres.csv     creixement observat i predit (validació creuada) de cada arbre
    data/validacio_parceles.csv   AB observada i predita a la mesura final de cada parcel·la de pinar i període

Els models s'ajusten amb les dues remesures juntes (IFN2-IFN3 i IFN3-IFN4) i un terme de període, p34
(1 = IFN3-IFN4). Per a les prediccions a Andorra es pren com a referència el període més recent: els
coeficients de calibracio.json porten l'efecte de p34 sumat a la constant, de manera que el simulador no
necessita saber el període. El coeficient de p34, el seu error estàndard i el període de referència es
desen a part («periode»), i els coeficients originals a «ajust».

Efectes aleatoris del creixement. El model de creixement té dues constants aleatòries niades: una de parcel·la
(u_p, comuna als dos períodes: qualitat de lloc) i una de parcel·la×període (v_pt: el que una parcel·la creix de
més o de menys en un període concret per raons pròpies, com el clima local, danys o gestió). El període continua
sent un efecte fix (p34). Les prediccions per a Andorra fan servir u = v = 0. Per a una parcel·la nova en un
període nou, la variància entre parcel·les és σu² + σv². La prova de raó de versemblances contra el model amb
només u_p es fa amb els dos models ajustats per màxima versemblança (ML); els coeficients i les variàncies que es
desen són els de l'ajust REML del model niat.

Mortalitat de fons i pertorbacions. L'IFN anota l'agent causant de cada arbre mort (scripts/preparar_ifn.py,
CAUSES). El model de mortalitat s'ajusta només amb les morts de fons (competència, causa desconeguda o sense
agent). Els arbres morts per una pertorbació (foc, sequera, neu i vent, causes biòtiques, altres) es treuen de
l'ajust: es tracten com a censurats, és a dir, se sap que no han mort de causa de fons durant el temps que
han viscut, però no quant de temps ha estat. És la simplificació habitual de riscos competitius quan no se
sap el moment de la mort; com que són pocs (un 5 % dels pins), sobreestimar-ne o subestimar-ne l'exposició
gairebé no canvia els coeficients.

Les pertorbacions s'estimen a part, a escala de parcel·la (estimar_pertorbacions), i es desen al bloc
«pertorbacions» de cada espècie: l'app les aplica com a escenaris.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

ARREL = Path(__file__).resolve().parent.parent
DATA = ARREL / "data"
sys.path.insert(0, str(ARREL))
from carboni.creixement import (  # noqa: E402
    EXPONENT_REINEKE,
    FORMULA_CREIXEMENT,
    FORMULA_MORTALITAT,
    PERIODE_REFERENCIA,
    TERME_PERIODE,
    TIPUS_PERTORBACIO,
    ModelEspecie,
    sdi,
    variables_arbre,
)
from carboni.equacions import alcada  # noqa: E402
from carboni.inventari import carregar_inf1and, resum_parceles  # noqa: E402

PLECS = 5
QUANTIL_SDI_MAX = 0.99
CODIS = {"Pinus sylvestris": "Ps", "Pinus uncinata": "Pu"}
PERIODES = ["IFN2-IFN3", "IFN3-IFN4"]
FORMULA_CREIXEMENT_PERIODE = f"{FORMULA_CREIXEMENT} + {TERME_PERIODE}"
FORMULA_MORTALITAT_PERIODE = f"{FORMULA_MORTALITAT} + {TERME_PERIODE}"
# Constant aleatòria de parcel·la×període, niada dins la de parcel·la (un nivell per període).
VC_PERIODE = {"periode": "0 + C(periode)"}
# Equacions de l'informe de carboni de l'INF1-AND (taules 8 i 9): IFN3-Lleida i IFN2-Lleida.
VOLUM_INFORME = {"Pinus sylvestris": (30.83, 0.000321), "Pinus uncinata": (77.94, 0.0002743)}
IAVC_INFORME = {"Pinus sylvestris": (3.9, 0.0162104, -0.0000019), "Pinus uncinata": (3.32, 0.0105025, -0.0000019)}
CLASSES_AB = [0, 20, 30, 40, 100]
# Una parcel·la és «afectada» per un tipus de pertorbació si aquesta causa hi mata almenys el 10 % de l'AB
# inicial de pins de l'espècie: separa els episodis de les morts aïllades.
LLINDAR_AFECTADA = 0.10
# Validació de parcel·les «sense pertorbacions»: menys de 2 m²/ha d'AB inicial morta per pertorbacions.
LLINDAR_VALIDACIO = 2.0
CLAU_PARCELA = ["periode", "parcela", "especie"]


def model_creixement(dades: pd.DataFrame, niat: bool = True):
    """Model mixt de creixement amb terme de període. niat: constant aleatòria de parcel·la (u_p) i de
    parcel·la×període (v_pt); si no, només u_p (model anterior, per a la prova de raó de versemblances)."""
    aleatoris = {"re_formula": "1", "vc_formula": VC_PERIODE} if niat else {}
    return smf.mixedlm(FORMULA_CREIXEMENT_PERIODE, dades, groups=dades["parcela"], **aleatoris)


def ajustar_creixement(dades: pd.DataFrame):
    """Model mixt niat (REML) amb terme de període. Retorna el model, els efectes fixos i la correcció de biaix de
    la retransformació logarítmica (estimador de raó de Snowdon, amb les prediccions dels efectes fixos)."""
    model = model_creixement(dades).fit(reml=True)
    fixos = {k: float(v) for k, v in model.fe_params.items()}
    pred_log = model.model.exog @ model.fe_params.to_numpy()
    correccio = float((dades["dd10"] + 1).mean() / np.exp(pred_log).mean())
    return model, fixos, correccio


def variancies_creixement(model) -> dict:
    """Components de variància (REML) del model niat i informació de convergència."""
    var_u, var_v = float(model.cov_re.iloc[0, 0]), float(model.vcomp[0])
    return {"sigma_u": float(np.sqrt(var_u)), "sigma_v": float(np.sqrt(var_v)), "sigma_e": float(np.sqrt(model.scale)),
            "var_u": var_u, "var_v": var_v, "var_e": float(model.scale),
            # Parcel·la nova en un període nou: u_p i v_pt sense informació.
            "sigma_entre_parceles": float(np.sqrt(var_u + var_v)),
            "convergit": bool(model.converged), "frontera_v": bool(var_v < 1e-8),
            "logvers_reml": float(model.llf)}


def prova_efecte_niat(dades: pd.DataFrame) -> dict:
    """Prova de raó de versemblances (ML) del model niat (u_p + v_pt) contra el de només u_p. La variància de v_pt
    és 0 sota la hipòtesi nul·la (a la frontera): el p de referència és el de la barreja 50:50 de χ²(0) i χ²(1),
    que és la meitat del de la χ²(1)."""
    niat, simple = model_creixement(dades).fit(reml=False), model_creixement(dades, niat=False).fit(reml=False)
    lr = float(2.0 * (niat.llf - simple.llf))
    p_chi2 = float(stats.chi2.sf(lr, 1))
    return {"metode": "ML", "logvers_niat": float(niat.llf), "logvers_parcela": float(simple.llf),
            "aic_niat": float(niat.aic), "aic_parcela": float(simple.aic), "lr": lr, "gl": 1,
            "p_chi2_1": p_chi2, "p_frontera": 0.5 * p_chi2,
            "convergit": {"niat": bool(niat.converged), "parcela": bool(simple.converged)}}


def ajustar_mortalitat(dades: pd.DataFrame):
    """GLM binomial cloglog de la mortalitat de fons. Els morts per pertorbacions es treuen (censurats)."""
    dades = fons(dades)
    familia = sm.families.Binomial(link=sm.families.links.CLogLog())
    return smf.glm(FORMULA_MORTALITAT_PERIODE, dades, family=familia, offset=np.log(dades["anys"])).fit()


def plegar_periode(coeficients: dict) -> dict:
    """Coeficients per al període de referència (IFN3-IFN4, p34 = 1): l'efecte del període passa a la constant."""
    c = dict(coeficients)
    c["Intercept"] += c.pop(TERME_PERIODE) * (1.0 if PERIODE_REFERENCIA == "IFN3-IFN4" else 0.0)
    return c


def resum_ajust(params: pd.Series, ee: pd.Series, p: pd.Series) -> dict:
    return {"coeficients": {k: float(v) for k, v in params.items()},
            "ee": {k: float(v) for k, v in ee.items()},
            "p": {k: float(v) for k, v in p.items()}}


def preparar(dades: pd.DataFrame) -> pd.DataFrame:
    dades = variables_arbre(dades)
    dades[TERME_PERIODE] = (dades["periode"] == "IFN3-IFN4").astype(int)
    dades["dd10"] = dades["id_dbh"] * 10.0
    dades["y"] = np.log(dades["dd10"] + 1.0)
    # Morts per pertorbacions: fora de l'ajust de mortalitat. m = mort de fons.
    dades["pertorbat"] = dades["mort"] & (dades["causa"].fillna("fons") != "fons")
    dades["m"] = (dades["mort"] & ~dades["pertorbat"]).astype(int)
    dades["g"] = np.pi * (dades["dbh"] / 200.0) ** 2 * dades["f"]
    return dades


def fons(dades: pd.DataFrame) -> pd.DataFrame:
    """Arbres que compten per a la mortalitat de fons: tots menys els morts per pertorbacions."""
    return dades[~dades["pertorbat"]]


def unitats_parcela(dades: pd.DataFrame) -> pd.DataFrame:
    """Per parcel·la, període i espècie: AB inicial de pins (g0), durada de l'interval i AB inicial dels pins
    morts per cada causa (m²/ha)."""
    u = dades.groupby(CLAU_PARCELA).agg(g0=("g", "sum"), anys=("anys", "first"))
    morts = dades[dades["mort"]].pivot_table(index=CLAU_PARCELA, columns="causa", values="g", aggfunc="sum")
    return u.join(morts.reindex(columns=["fons", *TIPUS_PERTORBACIO])).fillna(0.0)


def estimar_pertorbacions(dades: pd.DataFrame) -> dict:
    """Freqüència i severitat de cada tipus de pertorbació a les parcel·les (totes les parcel·les remesurades on
    hi ha l'espècie, no només les pinedes: les pertorbacions són rares i calen totes les parcel·les; la fracció es
    mesura sobre l'AB de pins de l'espècie, de manera que no depèn de la composició).

    - taxa_perdua_ab: AB inicial morta per la causa / suma de (AB inicial · anys). És la fracció de l'AB de
      pins que es perd cada any de mitjana, comptant totes les morts d'aquesta causa (també les aïllades).
    - parceles_afectades: parcel·les on la causa mata almenys LLINDAR_AFECTADA de l'AB inicial.
    - taxa_episodis: parcel·les afectades / anys-parcel·la observats (taxa de Poisson); p_anual = 1 − e^(−taxa)
      és la probabilitat anual que una parcel·la sigui afectada.
    - fraccio_mitjana i fraccio_mediana: fracció de l'AB inicial morta a les parcel·les afectades.
    """
    u = unitats_parcela(dades)
    exposicio = float((u["g0"] * u["anys"]).sum())
    anys_parcela = float(u["anys"].sum())
    causes = {}
    for causa in TIPUS_PERTORBACIO:
        fraccio = u[causa] / u["g0"]
        afectada = fraccio >= LLINDAR_AFECTADA
        taxa = float(afectada.sum()) / anys_parcela
        causes[causa] = {
            "taxa_perdua_ab": float(u[causa].sum()) / exposicio,
            "parceles_afectades": int(afectada.sum()),
            "parceles_amb_morts": int((u[causa] > 0).sum()),
            "taxa_episodis": taxa,
            "p_anual": float(1.0 - np.exp(-taxa)),
            "fraccio_mitjana": float(fraccio[afectada].mean()) if afectada.any() else None,
            "fraccio_mediana": float(fraccio[afectada].median()) if afectada.any() else None,
        }
    return {
        "parceles": int(len(u)),
        "anys_parcela": anys_parcela,
        "causes": causes,
        "taxa_perdua_ab_total": float(sum(c["taxa_perdua_ab"] for c in causes.values())),
        "taxa_perdua_ab_fons_observada": float(u["fons"].sum()) / exposicio,
    }


def resum_pertorbacions(dades: pd.DataFrame) -> dict:
    """Estimació de pertorbacions amb tots dos períodes junts, i per període com a referència."""
    return estimar_pertorbacions(dades) | {
        "per_periode": {p: {k: v for k, v in estimar_pertorbacions(d).items() if k != "causes"}
                        for p, d in dades.groupby("periode")}}


def morts_per_causa(dades: pd.DataFrame) -> dict:
    morts = dades[dades["mort"]]
    return {"per_causa": {c: int(n) for c, n in morts["causa"].value_counts().items()},
            "sense_agent": int(morts["agent"].isna().sum())}


def calibratge_mortalitat(observat: np.ndarray, predit: np.ndarray, grups: int = 10) -> list[dict]:
    """Morts observades i predites per decils de probabilitat predita."""
    t = pd.DataFrame({"o": observat, "p": predit})
    t["decil"] = pd.qcut(t["p"].rank(method="first"), grups, labels=False)
    r = t.groupby("decil").agg(n=("o", "size"), p_mitjana=("p", "mean"), observats=("o", "sum"), predits=("p", "sum"))
    return [{"n": int(f.n), "p_mitjana": float(f.p_mitjana), "observats": int(f.observats),
             "predits": float(f.predits)} for f in r.itertuples()]


def vius(dades: pd.DataFrame) -> pd.DataFrame:
    return dades[~dades["mort"] & (dades["dd10"] > -1)]


def metriques(observat: np.ndarray, predit: np.ndarray) -> dict:
    residu = predit - observat
    ss_tot = float(np.sum((observat - observat.mean()) ** 2))
    return {
        "n": int(len(observat)),
        "biaix": float(residu.mean()),
        "rmse": float(np.sqrt(np.mean(residu**2))),
        "r2": float(1.0 - np.sum(residu**2) / ss_tot),
        "mitjana_observada": float(observat.mean()),
    }


def validar_especie(dades: pd.DataFrame) -> pd.DataFrame:
    """Validació creuada per parcel·les: cada parcel·la (amb els seus dos períodes, si en té) es prediu amb
    models que no l'han vist (model niat; la parcel·la de prova té u = v = 0). Les prediccions fan servir el
    coeficient del període de cada arbre. La mortalitat
    de fons i la taxa de pèrdua per pertorbacions també s'estimen sense la parcel·la."""
    pred_creix = pd.Series(np.nan, index=dades.index)
    pred_mort = pd.Series(np.nan, index=dades.index)
    pred_pert = pd.Series(np.nan, index=dades.index)
    plecs = GroupKFold(n_splits=PLECS)
    for entrenament, prova in plecs.split(dades, groups=dades["parcela"]):
        entren = dades.iloc[entrenament]
        prova_df = dades.iloc[prova]
        _, fixos, correccio = ajustar_creixement(vius(entren))
        mortalitat = {k: float(v) for k, v in ajustar_mortalitat(entren).params.items()}
        model = ModelEspecie(especie="", creixement=fixos, correccio_biaix=correccio, mortalitat=mortalitat,
                             sdi_max=np.inf)
        pred_creix.loc[prova_df.index] = model.increment_diametral(prova_df)
        p_anual = model.mortalitat_anual(prova_df)
        pred_mort.loc[prova_df.index] = 1.0 - (1.0 - p_anual) ** prova_df["anys"]
        taxa = estimar_pertorbacions(entren)["taxa_perdua_ab_total"]
        pred_pert.loc[prova_df.index] = 1.0 - (1.0 - taxa) ** prova_df["anys"]
    return dades.assign(id_predit=pred_creix, p_mort_periode=pred_mort, p_pert_periode=pred_pert)


def metriques_parceles(parceles: pd.DataFrame, predita: str) -> dict:
    """AB final i increment d'AB de les parcel·les, observats vs predits (columna `predita`)."""
    return {
        "parceles": metriques(parceles["ab_2001_observada"].to_numpy(), parceles[predita].to_numpy()),
        "increment_ab_parceles": metriques((parceles["ab_2001_observada"] - parceles["ab_1990"]).to_numpy(),
                                           (parceles[predita] - parceles["ab_1990"]).to_numpy()),
    }


def metriques_validacio(validat: pd.DataFrame, parceles: pd.DataFrame) -> dict:
    """Mortalitat: només la de fons (els morts per pertorbacions no hi compten). Parcel·les: (1) totes, amb el
    model complet (mortalitat de fons i pèrdua esperada per pertorbacions); (2) totes, només amb la mortalitat de
    fons; (3) les que perden menys de LLINDAR_VALIDACIO m²/ha per pertorbacions, només amb la de fons."""
    v = vius(validat)
    f = fons(validat)
    mort_obs = f["m"].to_numpy()
    mort_pred = f["p_mort_periode"].to_numpy()
    estables = parceles[parceles["perdua_pertorbacions"] <= LLINDAR_VALIDACIO]
    return {
        "creixement": metriques(v["id_dbh"].to_numpy(), v["id_predit"].to_numpy()),
        "mortalitat": {
            "tipus": "fons",
            "n": int(len(f)),
            "censurats_pertorbacio": int(validat["pertorbat"].sum()),
            "auc": float(roc_auc_score(mort_obs, mort_pred)),
            "morts_observats": int(mort_obs.sum()),
            "morts_predits": float(mort_pred.sum()),
            "calibratge_decils": calibratge_mortalitat(mort_obs, mort_pred),
        },
        **metriques_parceles(parceles, "ab_2001_predita"),
        "nomes_fons": metriques_parceles(parceles, "ab_2001_predita_fons"),
        "sense_pertorbacions": metriques_parceles(estables, "ab_2001_predita_fons")
        | {"llindar_m2ha": LLINDAR_VALIDACIO, "parceles_excloses": int(len(parceles) - len(estables))},
    }


def validar_parceles(dades: pd.DataFrame) -> pd.DataFrame:
    """AB a la mesura final dels pins presents a la inicial: observada vs predita (creixement i mortalitat
    esperats). Per compatibilitat, les columnes es diuen ab_1990 (mesura inicial) i ab_2001_* (final).
    ab_2001_predita inclou la pèrdua esperada per pertorbacions; ab_2001_predita_fons, només la mortalitat de
    fons. perdua_pertorbacions és l'AB inicial dels pins morts per pertorbacions (observada)."""
    d = dades.copy()
    d["g_1990"] = np.pi * (d["dbh"] / 200.0) ** 2 * d["f"]
    d["g_2001_obs"] = np.where(d["mort"], 0.0, np.pi * (d["dbh_ifn3"] / 200.0) ** 2 * d["f"])
    dbh_pred = d["dbh"] + d["id_predit"] * d["anys"]
    d["g_2001_fons"] = np.pi * (dbh_pred / 200.0) ** 2 * d["f"] * (1.0 - d["p_mort_periode"])
    d["g_2001_pred"] = d["g_2001_fons"] * (1.0 - d["p_pert_periode"])
    d["g_pert"] = np.where(d["pertorbat"], d["g_1990"], 0.0)
    parceles = d.groupby(["periode", "parcela", "especie"]).agg(
        ab_1990=("g_1990", "sum"),
        ab_2001_observada=("g_2001_obs", "sum"),
        ab_2001_predita=("g_2001_pred", "sum"),
        ab_2001_predita_fons=("g_2001_fons", "sum"),
        perdua_pertorbacions=("g_pert", "sum"),
        ab_parcela=("ab", "first"),
        altitud=("altitud", "first"),
        anys=("anys", "first"),
    ).reset_index()
    # Només pinars: els pins de l'espècie fan almenys el 70 % de l'AB de la parcel·la.
    return parceles[parceles["ab_1990"] >= 0.7 * parceles["ab_parcela"]].round(4)


def comparar_amb_informe(especie: str, validat: pd.DataFrame) -> list[dict]:
    """Creixement en volum observat vs predit pel model i per les equacions IAVC de l'informe, per densitat."""
    v = validat[~validat["mort"]].copy()
    a, b = VOLUM_INFORME[especie]
    ia, ib, ic = IAVC_INFORME[especie]
    sp = pd.Series(CODIS[especie], index=v.index)

    def volum(dbh: pd.Series) -> np.ndarray:
        return a + b * (dbh.to_numpy() * 10.0) ** 2 * alcada(sp, dbh, "m")

    v0 = volum(v["dbh"])
    v["observat"] = (volum(v["dbh_ifn3"]) - v0) / v["anys"] * v["f"]
    v["model"] = (volum(v["dbh"] + v["id_predit"] * v["anys"]) - v0) / v["anys"] * v["f"]
    v["informe"] = (ia + ib * v0 + ic * v0**2) * v["f"]
    v["classe"] = pd.cut(v["ab"], CLASSES_AB)
    taula = v.groupby("classe", observed=True).agg(
        arbres=("observat", "size"), observat=("observat", "sum"), model=("model", "sum"), informe=("informe", "sum")
    )
    return [
        {"ab_min": float(c.left), "ab_max": float(c.right), "arbres": int(f.arbres),
         "model_sobre_observat": float(f.model / f.observat), "informe_sobre_observat": float(f.informe / f.observat)}
        for c, f in taula.iterrows()
    ]


def sdi_maxim(especie: str, lleida: pd.DataFrame, andorra: pd.DataFrame) -> dict:
    """Percentil 99 de l'SDI de les pinedes d'Andorra i de l'estat inicial de les parcel·les de Lleida (cada
    parcel·la hi compta un cop per període: l'IFN2 i l'IFN3 són dos estats diferents del rodal)."""
    and_ = andorra[(andorra["especie_dominant"] == CODIS[especie]) & (andorra["fraccio_dominant"] >= 0.7)]
    sdi_and = sdi(and_["n"], and_["dq"])
    parc = lleida[lleida["especie"] == especie].groupby(["periode", "parcela"]).agg(n=("n", "first"), ab=("ab", "first"))
    dq = np.sqrt(4 * parc["ab"] / (np.pi * parc["n"])) * 100
    sdi_lle = sdi(parc["n"], dq)
    sdi_and, sdi_lle = pd.Series(sdi_and), pd.Series(sdi_lle)
    tots = pd.concat([sdi_and, sdi_lle])
    return {
        "sdi_max": float(tots.quantile(QUANTIL_SDI_MAX)),
        "sdi_p50_andorra": float(sdi_and.median()),
        "sdi_p50_lleida": float(sdi_lle.median()),
        "n_parceles": int(len(tots)),
    }


def mostra(dades: pd.DataFrame) -> dict:
    return {"arbres": int(len(dades)), "parceles": int(dades["parcela"].nunique()), "morts": int(dades["mort"].sum()),
            "morts_fons": int(dades["m"].sum()), "morts_pertorbacio": int(dades["pertorbat"].sum()),
            "arbres_creixement": int(len(vius(dades))), "anys_min": int(dades["anys"].min()),
            "anys_max": int(dades["anys"].max())}


def main() -> None:
    lleida = preparar(pd.read_csv(DATA / "ifn_lleida_arbres.csv"))
    andorra = resum_parceles(*carregar_inf1and(DATA / "inf1and"))
    comarques = sorted(lleida["comarca"].unique())

    sortida = {
        "exponent_reineke": EXPONENT_REINEKE,
        "dades_calibratge": {
            "font": "IFN2, IFN3 i IFN4 de Lleida (MITECO), parcel·les remesurades sense tallades",
            "comarques": comarques,
            "periodes": PERIODES,
            "terme_periode": TERME_PERIODE,
            "periode_referencia": PERIODE_REFERENCIA,
        },
        "pertorbacions": {
            "descripcio": "Mortalitat per pertorbacions (agent causant de l'IFN), estimada a part de la de fons. "
                          "Fracció de l'AB inicial de pins de l'espècie; tots dos períodes junts; totes les "
                          "parcel·les remesurades on hi ha l'espècie.",
            "tipus": list(TIPUS_PERTORBACIO),
            "llindar_afectada": LLINDAR_AFECTADA,
            "totes_especies": resum_pertorbacions(lleida),
        },
        "especies": {},
    }
    validacions_arbres, validacions_parceles = [], []
    for especie, dades in lleida.groupby("especie"):
        model_c, fixos, correccio = ajustar_creixement(vius(dades))
        model_m = ajustar_mortalitat(dades)
        validat = validar_especie(dades)
        parceles = validar_parceles(validat)
        metr = metriques_validacio(validat, parceles)
        metr["comparacio_informe"] = comparar_amb_informe(especie, validat)
        metr["per_periode"] = {
            periode: metriques_validacio(validat[validat["periode"] == periode], parceles[parceles["periode"] == periode])
            | {"comparacio_informe": comparar_amb_informe(especie, validat[validat["periode"] == periode])}
            for periode in PERIODES
        }
        ee_c, ee_m = model_c.bse_fe[TERME_PERIODE], model_m.bse[TERME_PERIODE]
        sortida["especies"][especie] = {
            "creixement": plegar_periode(fixos),
            "correccio_biaix": correccio,
            "mortalitat": plegar_periode({k: float(v) for k, v in model_m.params.items()}),
            **sdi_maxim(especie, lleida, andorra),
            "pertorbacions": resum_pertorbacions(dades),
            "periode": {
                "terme": TERME_PERIODE,
                "codificacio": "p34 = 1 a IFN3-IFN4, 0 a IFN2-IFN3",
                "referencia": PERIODE_REFERENCIA,
                "creixement": {"coef": fixos[TERME_PERIODE], "ee": float(ee_c),
                               "p": float(model_c.pvalues[TERME_PERIODE]),
                               "efecte_relatiu": float(np.exp(fixos[TERME_PERIODE]) - 1.0)},
                "mortalitat": {"coef": float(model_m.params[TERME_PERIODE]), "ee": float(ee_m),
                               "p": float(model_m.pvalues[TERME_PERIODE]),
                               "rao_de_riscos": float(np.exp(model_m.params[TERME_PERIODE]))},
            },
            "ajust": {
                "creixement": resum_ajust(model_c.fe_params, model_c.bse_fe, model_c.pvalues[model_c.fe_params.index])
                | {"efectes_aleatoris": "u_p (parcel·la) + v_pt (parcel·la×període), REML",
                   **variancies_creixement(model_c), "prova_lr_parcela_periode": prova_efecte_niat(vius(dades))},
                "mortalitat": resum_ajust(model_m.params, model_m.bse, model_m.pvalues)
                | {"convergit": bool(model_m.converged)},
            },
            "calibratge": {
                **mostra(dades),
                "altitud_min": float(dades["altitud"].min()),
                "altitud_max": float(dades["altitud"].max()),
                "ab_max": float(dades["ab"].max()),
                "parceles_dos_periodes": int(dades.groupby("parcela")["periode"].nunique().eq(2).sum()),
                "comarques": comarques,
                "per_periode": {p: mostra(d) | {"morts_causa": morts_per_causa(d)}
                                for p, d in dades.groupby("periode")},
            },
            "validacio": metr,
        }
        validacions_arbres.append(validat.loc[~validat["mort"], ["periode", "especie", "parcela", "dbh", "bal", "ab",
                                                                 "altitud", "id_dbh", "id_predit"]])
        validacions_parceles.append(parceles)
        print(especie, json.dumps({k: metr[k] for k in ("creixement", "mortalitat", "parceles", "sense_pertorbacions")},
                                  indent=1, ensure_ascii=False))
        print("  pertorbacions:",
              json.dumps(sortida["especies"][especie]["pertorbacions"]["causes"], ensure_ascii=False))
        print("  període:", json.dumps(sortida["especies"][especie]["periode"], ensure_ascii=False))
        print("  efectes aleatoris:", json.dumps(sortida["especies"][especie]["ajust"]["creixement"]["sigma_u"]),
              json.dumps(sortida["especies"][especie]["ajust"]["creixement"]["sigma_v"]),
              json.dumps(sortida["especies"][especie]["ajust"]["creixement"]["prova_lr_parcela_periode"]["lr"]))

    (DATA / "calibracio.json").write_text(json.dumps(sortida, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.concat(validacions_arbres).round(4).to_csv(DATA / "validacio_arbres.csv", index=False)
    pd.concat(validacions_parceles).to_csv(DATA / "validacio_parceles.csv", index=False)


if __name__ == "__main__":
    main()
