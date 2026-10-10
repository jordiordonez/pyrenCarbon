"""Calibra els models de creixement diametral i mortalitat i els valida.

Ús:
    python scripts/calibrar.py

Llegeix data/ifn_lleida_arbres.csv (scripts/preparar_ifn.py) i data/inf1and/ (INF1-AND) i escriu:
    data/calibracio.json          coeficients, límits d'autoaclarida i mètriques de validació
    data/validacio_arbres.csv     creixement observat i predit (validació creuada) de cada arbre
    data/validacio_parceles.csv   AB observada i predita el 2001 a cada parcel·la de pinar
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

ARREL = Path(__file__).resolve().parent.parent
DATA = ARREL / "data"
sys.path.insert(0, str(ARREL))
from carboni.creixement import (  # noqa: E402
    EXPONENT_REINEKE,
    FORMULA_CREIXEMENT,
    FORMULA_MORTALITAT,
    ModelEspecie,
    sdi,
    variables_arbre,
)
from carboni.equacions import alcada  # noqa: E402
from carboni.inventari import carregar_inf1and, resum_parceles  # noqa: E402

PLECS = 5
QUANTIL_SDI_MAX = 0.99
CODIS = {"Pinus sylvestris": "Ps", "Pinus uncinata": "Pu"}
# Equacions de l'informe de carboni de l'INF1-AND (taules 8 i 9): IFN3-Lleida i IFN2-Lleida.
VOLUM_INFORME = {"Pinus sylvestris": (30.83, 0.000321), "Pinus uncinata": (77.94, 0.0002743)}
IAVC_INFORME = {"Pinus sylvestris": (3.9, 0.0162104, -0.0000019), "Pinus uncinata": (3.32, 0.0105025, -0.0000019)}
CLASSES_AB = [0, 20, 30, 40, 100]


def ajustar_creixement(dades: pd.DataFrame) -> tuple[dict, float]:
    model = smf.mixedlm(FORMULA_CREIXEMENT, dades, groups=dades["parcela"]).fit(reml=True)
    fixos = {k: float(v) for k, v in model.fe_params.items()}
    # Correcció de biaix de la retransformació logarítmica (estimador de raó de Snowdon).
    pred_log = model.model.exog @ model.fe_params.to_numpy()
    correccio = float((dades["dd10"] + 1).mean() / np.exp(pred_log).mean())
    return fixos, correccio


def ajustar_mortalitat(dades: pd.DataFrame) -> dict:
    familia = sm.families.Binomial(link=sm.families.links.CLogLog())
    model = smf.glm(FORMULA_MORTALITAT, dades, family=familia, offset=np.log(dades["anys"])).fit()
    return {k: float(v) for k, v in model.params.items()}


def preparar(dades: pd.DataFrame) -> pd.DataFrame:
    dades = variables_arbre(dades)
    dades["dd10"] = dades["id_dbh"] * 10.0
    dades["y"] = np.log(dades["dd10"] + 1.0)
    dades["m"] = dades["mort"].astype(int)
    return dades


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


def validar_especie(dades: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Validació creuada per parcel·les: cada parcel·la es prediu amb models que no l'han vist."""
    vius = dades[~dades["mort"] & (dades["dd10"] > -1)]
    pred_creix = pd.Series(np.nan, index=dades.index)
    pred_mort = pd.Series(np.nan, index=dades.index)
    plecs = GroupKFold(n_splits=PLECS)
    for entrenament, prova in plecs.split(dades, groups=dades["parcela"]):
        entren = dades.iloc[entrenament]
        prova_df = dades.iloc[prova]
        fixos, correccio = ajustar_creixement(entren[~entren["mort"] & (entren["dd10"] > -1)])
        mortalitat = ajustar_mortalitat(entren)
        model = ModelEspecie(especie="", creixement=fixos, correccio_biaix=correccio, mortalitat=mortalitat,
                             sdi_max=np.inf)
        pred_creix.loc[prova_df.index] = model.increment_diametral(prova_df)
        p_anual = model.mortalitat_anual(prova_df)
        pred_mort.loc[prova_df.index] = 1.0 - (1.0 - p_anual) ** prova_df["anys"]
    dades = dades.assign(id_predit=pred_creix, p_mort_periode=pred_mort)

    v = dades.loc[vius.index]
    mort_obs = dades["m"].to_numpy()
    mort_pred = dades["p_mort_periode"].to_numpy()
    resultat = {
        "creixement": metriques(v["id_dbh"].to_numpy(), v["id_predit"].to_numpy()),
        "mortalitat": {
            "n": int(len(dades)),
            "auc": float(roc_auc_score(mort_obs, mort_pred)),
            "morts_observats": int(mort_obs.sum()),
            "morts_predits": float(mort_pred.sum()),
        },
    }
    return dades, resultat


def validar_parceles(dades: pd.DataFrame) -> pd.DataFrame:
    """AB del 2001 dels pins presents el 1990: observada vs predita (creixement i mortalitat esperats)."""
    d = dades.copy()
    d["g_1990"] = np.pi * (d["dbh"] / 200.0) ** 2 * d["f"]
    d["g_2001_obs"] = np.where(d["mort"], 0.0, np.pi * (d["dbh_ifn3"] / 200.0) ** 2 * d["f"])
    dbh_pred = d["dbh"] + d["id_predit"] * d["anys"]
    d["g_2001_pred"] = np.pi * (dbh_pred / 200.0) ** 2 * d["f"] * (1.0 - d["p_mort_periode"])
    parceles = d.groupby(["parcela", "especie"]).agg(
        ab_1990=("g_1990", "sum"),
        ab_2001_observada=("g_2001_obs", "sum"),
        ab_2001_predita=("g_2001_pred", "sum"),
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
    codi = {"Pinus sylvestris": "Ps", "Pinus uncinata": "Pu"}[especie]
    and_ = andorra[(andorra["especie_dominant"] == codi) & (andorra["fraccio_dominant"] >= 0.7)]
    sdi_and = sdi(and_["n"], and_["dq"])
    parc = lleida[lleida["especie"] == especie].groupby("parcela").agg(n=("n", "first"), ab=("ab", "first"))
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


def main() -> None:
    lleida = preparar(pd.read_csv(DATA / "ifn_lleida_arbres.csv"))
    andorra = resum_parceles(*carregar_inf1and(DATA / "inf1and"))

    sortida = {"exponent_reineke": EXPONENT_REINEKE, "especies": {}}
    validacions_arbres, validacions_parceles = [], []
    for especie, dades in lleida.groupby("especie"):
        vius = dades[~dades["mort"] & (dades["dd10"] > -1)]
        fixos, correccio = ajustar_creixement(vius)
        mortalitat = ajustar_mortalitat(dades)
        validat, metr = validar_especie(dades)
        parceles = validar_parceles(validat)
        metr["parceles"] = metriques(parceles["ab_2001_observada"].to_numpy(), parceles["ab_2001_predita"].to_numpy())
        metr["comparacio_informe"] = comparar_amb_informe(especie, validat)
        increment_obs = parceles["ab_2001_observada"] - parceles["ab_1990"]
        increment_pred = parceles["ab_2001_predita"] - parceles["ab_1990"]
        metr["increment_ab_parceles"] = metriques(increment_obs.to_numpy(), increment_pred.to_numpy())
        sortida["especies"][especie] = {
            "creixement": fixos,
            "correccio_biaix": correccio,
            "mortalitat": mortalitat,
            **sdi_maxim(especie, lleida, andorra),
            "calibratge": {
                "arbres": int(len(dades)),
                "parceles": int(dades["parcela"].nunique()),
                "morts": int(dades["mort"].sum()),
                "altitud_min": float(dades["altitud"].min()),
                "altitud_max": float(dades["altitud"].max()),
                "ab_max": float(dades["ab"].max()),
            },
            "validacio": metr,
        }
        validacions_arbres.append(
            validat.loc[~validat["mort"], ["especie", "parcela", "dbh", "bal", "ab", "altitud", "id_dbh", "id_predit"]]
        )
        validacions_parceles.append(parceles)
        print(especie, json.dumps(metr, indent=1, ensure_ascii=False))

    (DATA / "calibracio.json").write_text(json.dumps(sortida, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.concat(validacions_arbres).round(4).to_csv(DATA / "validacio_arbres.csv", index=False)
    pd.concat(validacions_parceles).to_csv(DATA / "validacio_parceles.csv", index=False)


if __name__ == "__main__":
    main()
