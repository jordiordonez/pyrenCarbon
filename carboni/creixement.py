"""Simulador de creixement arbre a arbre per a boscos de pi roig i de pi negre.

Cada any:
1. increment diametral de cada pi, en funció de la mida, la competència (BAL i AB) i l'altitud
   (model lineal mixt calibrat amb les remesures IFN2–IFN3 i IFN3–IFN4 de les comarques pirinenques de
   Lleida, amb un terme de període);
2. mortalitat de fons individual, en funció de la mida i la competència (GLM binomial cloglog, mateixes dades,
   ajustat només amb els pins morts per competència o causa desconeguda);
3. pertorbacions (foc, sequera, neu i vent, plagues…), segons l'escenari triat (Pertorbacio):
   - sense pertorbacions: només la mortalitat de fons;
   - risc mitjà: cada any es perd la fracció de pins que es va perdre de mitjana per pertorbacions a les
     parcel·les de calibratge (taxa_perdua_ab_total de calibracio.json). Al pi roig és igual per a totes les
     mides; al pi negre, el risc anual de cada pi creix amb el diàmetre (GLM cloglog, bloc «mida» de
     calibracio.json), amb la constant calibrada perquè la pèrdua d'AB a les parcel·les de calibratge sigui la
     mateixa taxa;
   - episodi puntual: l'any indicat mor una fracció (severitat) dels pins de totes les mides;
4. autoaclarida: si l'índex de densitat de Reineke supera el màxim observat, moren els pins més petits.

Els arbres d'altres espècies compten per a la competència però no creixen ni moren.

Els coeficients de data/calibracio.json corresponen al període de referència (IFN3–IFN4, el més recent):
l'efecte del terme de període ja és a la constant, de manera que el simulador no necessita cap entrada de
període.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd

from carboni.equacions import (
    CO2_PER_C,
    ESPECIES,
    FRACCIO_CARBONI,
    alcada,
    biomassa_aeria_kg,
    biomassa_subterrania_kg,
)
from carboni.inventari import Rodal, factor_alcada_local

FORMULA_CREIXEMENT = "y ~ lnd + dbh + bal_l + lnab + alt"
FORMULA_MORTALITAT = "m ~ inv_d + dbh + bal"
# Risc anual de mort per pertorbacions segons la mida (espècies on s'ha adoptat; scripts/calibrar.py).
FORMULA_PERTORBACIO_MIDA = "y ~ dbh"
# Terme de període del calibratge (1 = IFN3–IFN4) i període que fan servir les prediccions.
TERME_PERIODE = "p34"
PERIODE_REFERENCIA = "IFN3-IFN4"
EXPONENT_REINEKE = 1.605
DQ_REFERENCIA_REINEKE = 25.0
# Tipus de pertorbació (grups de causa de scripts/preparar_ifn.py, CAUSES, llevat de «fons»).
TIPUS_PERTORBACIO = {
    "foc": "Foc",
    "sequera": "Sequera",
    "neu_vent": "Neu o vent",
    "biotic": "Plaga o causa biòtica",
    "altres": "Altres (aprofitaments, despreniments, erosió)",
}
# Tipus que es poden simular com a episodi puntual.
EPISODIS = ("foc", "sequera", "neu_vent", "biotic")
MODES_PERTORBACIO = {"cap": "Sense pertorbacions", "risc_mitja": "Risc mitjà", "episodi": "Episodi puntual"}
# Amb menys parcel·les afectades, la severitat per defecte d'un episodi es pren de totes dues espècies juntes.
MIN_PARCELES_SEVERITAT = 3
# Cas pessimista de la banda d'incertesa: percentil de la pèrdua acumulada per pertorbacions en l'horitzó,
# estimat amb una simulació Monte Carlo (nombre de rodals simulats i llavor fixa, perquè sigui reproduïble).
QUANTIL_PESSIMISTA = 0.9
SIMULACIONS_MC = 20_000
LLAVOR_MC = 20261016


def variables_arbre(arbres: pd.DataFrame) -> pd.DataFrame:
    """Afegeix les variables explicatives a partir de dbh, bal, ab i altitud."""
    a = arbres.copy()
    a["lnd"] = np.log(a["dbh"])
    a["bal_l"] = a["bal"] / np.log(a["dbh"] + 1.0)
    a["lnab"] = np.log(a["ab"])
    a["alt"] = a["altitud"] / 1000.0
    a["inv_d"] = 1.0 / a["dbh"]
    return a


def sdi(n, dq):
    return np.asarray(n) * np.power(np.asarray(dq) / DQ_REFERENCIA_REINEKE, EXPONENT_REINEKE)


def _predictor_lineal(coeficients: dict[str, float], variables: pd.DataFrame) -> np.ndarray:
    resultat = np.full(len(variables), coeficients.get("Intercept", 0.0))
    for nom, valor in coeficients.items():
        if nom != "Intercept":
            resultat = resultat + valor * variables[nom].to_numpy(dtype=float)
    return resultat


@dataclass(frozen=True)
class ModelEspecie:
    especie: str
    creixement: dict[str, float]
    correccio_biaix: float
    mortalitat: dict[str, float]
    sdi_max: float

    def increment_diametral(self, arbres: pd.DataFrame) -> np.ndarray:
        """Increment diametral anual (cm/any)."""
        v = variables_arbre(arbres)
        dd10 = np.exp(_predictor_lineal(self.creixement, v)) * self.correccio_biaix - 1.0
        return np.maximum(dd10, 0.0) / 10.0

    def mortalitat_anual(self, arbres: pd.DataFrame) -> np.ndarray:
        v = variables_arbre(arbres)
        return 1.0 - np.exp(-np.exp(_predictor_lineal(self.mortalitat, v)))


def carregar_models(ruta: Path) -> dict[str, ModelEspecie]:
    dades = json.loads(Path(ruta).read_text(encoding="utf-8"))
    codi_per_nom = {nom: codi for codi, nom in ESPECIES.items()}
    return {
        codi_per_nom[nom]: ModelEspecie(
            especie=nom,
            creixement=p["creixement"],
            correccio_biaix=p["correccio_biaix"],
            mortalitat=p["mortalitat"],
            sdi_max=p["sdi_max"],
        )
        for nom, p in dades["especies"].items()
    }


@dataclass(frozen=True)
class Pertorbacio:
    """Escenari de pertorbacions. `taxa_anual`: fracció dels pins de cada espècie (codi) que mor cada any en el
    mode «risc_mitja»; `risc_mida`: per a les espècies que en tenen, coeficients del risc anual segons la mida
    (cloglog), que substitueixen la taxa uniforme; `escala`: multiplicador del risc anual de cada espècie (1 = risc
    mitjà; més gran en el cas pessimista). `any_episodi` i `severitat` (fracció dels pins que mor): mode
    «episodi»."""

    mode: str = "cap"
    taxa_anual: dict[str, float] = field(default_factory=dict)
    risc_mida: dict[str, dict[str, float]] = field(default_factory=dict)
    escala: dict[str, float] = field(default_factory=dict)
    tipus: str = ""
    any_episodi: int = 0
    severitat: float = 0.0

    def __post_init__(self):
        if self.mode not in MODES_PERTORBACIO:
            raise ValueError(f"Mode de pertorbació desconegut: {self.mode}")
        if not 0.0 <= self.severitat <= 1.0:
            raise ValueError("La severitat ha d'estar entre 0 i 1.")

    def supervivencia(self, any_: int, arbres: pd.DataFrame, es_pi: np.ndarray) -> np.ndarray:
        """Fracció de cada registre que sobreviu a les pertorbacions l'any `any_`. `arbres`: estat a l'inici de
        l'any, amb la competència (_competencia)."""
        sp = arbres["sp"].to_numpy()
        factor = np.ones(len(sp))
        if self.mode == "risc_mitja":
            for codi, taxa in self.taxa_anual.items():
                idx = (sp == codi) & es_pi
                coeficients = self.risc_mida.get(codi)
                escala = self.escala.get(codi, 1.0)
                if coeficients and idx.any():
                    risc = np.exp(_predictor_lineal(coeficients, variables_arbre(arbres.loc[idx])))
                    factor[idx] = np.exp(-escala * risc)
                else:
                    factor[idx] = (1.0 - taxa) ** escala
        elif self.mode == "episodi" and any_ == self.any_episodi:
            factor[es_pi] = 1.0 - self.severitat
        return factor


SENSE_PERTORBACIONS = Pertorbacio()


def carregar_pertorbacions(ruta: Path) -> dict[str, dict]:
    """Paràmetres de pertorbacions per espècie (codi): taxa anual total de pèrdua d'AB, coeficients del risc
    segons la mida (None si l'espècie fa servir la taxa uniforme) i, per tipus, probabilitat anual d'episodi i
    severitat per defecte (fracció mitjana d'AB morta a les parcel·les afectades; si l'espècie en té menys de
    MIN_PARCELES_SEVERITAT, la de totes dues espècies juntes)."""
    dades = json.loads(Path(ruta).read_text(encoding="utf-8"))
    conjunt = dades["pertorbacions"]["totes_especies"]["causes"]
    codi_per_nom = {nom: codi for codi, nom in ESPECIES.items()}
    resultat = {}
    for nom, p in dades["especies"].items():
        pert = p["pertorbacions"]
        tipus = {}
        for t, c in pert["causes"].items():
            propi = c["parceles_afectades"] >= MIN_PARCELES_SEVERITAT
            font = c if propi else conjunt[t]
            tipus[t] = {"p_anual": c["p_anual"], "taxa_episodis": c["taxa_episodis"],
                        "severitat": font["fraccio_mitjana"] or 0.0, "fraccions": font.get("fraccions", []),
                        "parceles_afectades": font["parceles_afectades"], "severitat_conjunta": not propi}
        mida = pert.get("mida") or {}
        taxa = pert["taxa_perdua_ab_total"]
        episodica = sum(t["taxa_episodis"] * np.mean(t["fraccions"]) for t in tipus.values() if t["fraccions"])
        resultat[codi_per_nom[nom]] = {"taxa_anual": taxa, "tipus": tipus,
                                       "taxa_cronica": max(0.0, taxa - episodica),
                                       "risc_mida": mida.get("coeficients") if mida.get("adoptat") else None}
    return resultat


def pertorbacio_risc_mitja(parametres: dict[str, dict]) -> Pertorbacio:
    return Pertorbacio("risc_mitja", taxa_anual={sp: p["taxa_anual"] for sp, p in parametres.items()},
                       risc_mida={sp: p["risc_mida"] for sp, p in parametres.items() if p.get("risc_mida")})


def perdua_acumulada_mc(parametres_especie: dict, anys: int, n: int = SIMULACIONS_MC,
                        llavor: int = LLAVOR_MC) -> np.ndarray:
    """Fracció de l'AB de pins perduda per pertorbacions en `anys` anys a n rodals simulats (Monte Carlo).

    Model: els episodis de cada tipus arriben com un procés de Poisson (taxa_episodis, episodis per rodal i any) i
    cada episodi mata una fracció de l'AB presa a l'atzar de les observades a les parcel·les afectades de Lleida
    (distribució empírica). A més, cada any es perd la part «crònica» de la taxa (morts aïllades, per sota del
    llindar d'episodi), de manera que la pèrdua mitjana anual és la del risc mitjà."""
    rng = np.random.default_rng(llavor)
    log_supervivencia = np.full(n, anys * np.log1p(-parametres_especie["taxa_cronica"]))
    for tipus in parametres_especie["tipus"].values():
        fraccions = np.asarray(tipus["fraccions"], dtype=float)
        if not len(fraccions) or anys <= 0:
            continue
        episodis = rng.poisson(tipus["taxa_episodis"] * anys, n)
        severitat = np.minimum(rng.choice(fraccions, episodis.sum()), 1.0 - 1e-9)
        log_supervivencia += np.bincount(np.repeat(np.arange(n), episodis), weights=np.log1p(-severitat),
                                         minlength=n)
    return 1.0 - np.exp(log_supervivencia)


def escala_pessimista(parametres_especie: dict, anys: int, quantil: float = QUANTIL_PESSIMISTA) -> float:
    """Multiplicador del risc anual mitjà perquè la pèrdua acumulada en `anys` anys sigui el percentil `quantil`
    del Monte Carlo, repartida igual cada any. Mai per sota d'1: si en l'horitzó menys d'1 de cada 10 rodals pateix
    un episodi, el percentil 90 és més baix que la mitjana i el cas pessimista es queda en el risc mitjà."""
    if anys <= 0:
        return 1.0
    perdua = float(np.quantile(perdua_acumulada_mc(parametres_especie, anys), quantil))
    risc_mitja = -np.log1p(-parametres_especie["taxa_anual"])
    return max(1.0, -np.log1p(-perdua) / anys / risc_mitja)


def pertorbacio_pessimista(parametres: dict[str, dict], anys: int) -> Pertorbacio:
    """Risc mitjà amb el risc anual multiplicat perquè la pèrdua en l'horitzó sigui la del percentil 90."""
    return replace(pertorbacio_risc_mitja(parametres),
                   escala={sp: escala_pessimista(p, anys) for sp, p in parametres.items()})


def pertorbacio_episodi(tipus: str, any_episodi: int, severitat: float) -> Pertorbacio:
    if tipus not in EPISODIS:
        raise ValueError(f"Tipus d'episodi desconegut: {tipus}")
    return Pertorbacio("episodi", tipus=tipus, any_episodi=int(any_episodi), severitat=float(severitat))


def _competencia(arbres: pd.DataFrame) -> pd.DataFrame:
    a = arbres.copy()
    a["g"] = np.pi * (a["dbh"] / 200.0) ** 2 * a["f"]
    ordre = np.argsort(-a["dbh"].to_numpy(), kind="stable")
    g_ordenat = a["g"].to_numpy()[ordre]
    bal = np.empty(len(a))
    bal[ordre] = np.cumsum(g_ordenat) - g_ordenat
    a["bal"] = bal
    a["ab"] = a["g"].sum()
    return a


def _sdi_rodal(dbh: np.ndarray, f: np.ndarray) -> float:
    n = f.sum()
    if n <= 0:
        return 0.0
    dq = np.sqrt(np.sum(f * dbh**2) / n)
    return float(sdi(n, dq))


def _autoaclarida(arbres: pd.DataFrame, es_pi: np.ndarray, sdi_max: float) -> tuple[pd.DataFrame, float]:
    """Elimina densitat dels pins més petits fins que el rodal no supera sdi_max. Retorna peus/ha morts."""
    dbh = arbres["dbh"].to_numpy()
    f = arbres["f"].to_numpy().copy()
    if _sdi_rodal(dbh, f) <= sdi_max:
        return arbres, 0.0
    candidats = np.where(es_pi)[0]
    candidats = candidats[np.argsort(dbh[candidats], kind="stable")]
    disponible = f[candidats].sum()

    def f_amb_retirada(retirada: float) -> np.ndarray:
        nou = f.copy()
        acumulat = np.cumsum(f[candidats])
        anterior = acumulat - f[candidats]
        retirat = np.clip(retirada - anterior, 0.0, f[candidats])
        nou[candidats] = f[candidats] - retirat
        return nou

    baix, alt = 0.0, disponible
    for _ in range(40):
        mig = (baix + alt) / 2.0
        if _sdi_rodal(dbh, f_amb_retirada(mig)) > sdi_max:
            baix = mig
        else:
            alt = mig
    resultat = arbres.copy()
    resultat["f"] = f_amb_retirada(alt)
    return resultat, alt


def estoc(arbres: pd.DataFrame, habitat: str, factor_alcada: float) -> pd.DataFrame:
    """Afegeix alçada, biomassa i carboni per arbre (només pins)."""
    a = arbres[arbres["sp"].isin(ESPECIES)].copy()
    a["ht"] = alcada(a["sp"], a["dbh"], habitat) * factor_alcada
    a["bat"] = biomassa_aeria_kg(a["sp"], a["dbh"], a["ht"]) * a["f"] / 1000.0
    a["bst"] = biomassa_subterrania_kg(a["sp"], a["dbh"]) * a["f"] / 1000.0
    a["cat"] = a["bat"] * FRACCIO_CARBONI
    a["cst"] = a["bst"] * FRACCIO_CARBONI
    return a


def _resum_any(any_: int, arbres: pd.DataFrame, habitat: str, factor_alcada: float, morts_n: float,
               morts_c: float) -> dict:
    pins = estoc(arbres, habitat, factor_alcada)
    n = pins["f"].sum()
    ab = (np.pi * (pins["dbh"] / 200.0) ** 2 * pins["f"]).sum()
    dq = np.sqrt(4.0 * ab / (np.pi * n)) * 100.0 if n > 0 else 0.0
    c_total = pins["cat"].sum() + pins["cst"].sum()
    return {
        "any": any_,
        "densitat": n,
        "ab": ab,
        "dq": dq,
        "ho": float(pins.nlargest(max(1, int(round(len(pins) * 0.2))), "dbh")["ht"].mean()) if len(pins) else 0.0,
        "sdi": _sdi_rodal(arbres["dbh"].to_numpy(), arbres["f"].to_numpy()),
        "bat": pins["bat"].sum(),
        "bst": pins["bst"].sum(),
        "cat": pins["cat"].sum(),
        "cst": pins["cst"].sum(),
        "carboni": c_total,
        "co2": c_total * CO2_PER_C,
        "morts_peus": morts_n,
        "morts_carboni": morts_c,
    }  # morts_peus: mortalitat de fons, pertorbacions i autoaclarida


@dataclass(frozen=True)
class Simulacio:
    serie: pd.DataFrame
    arbres_inicials: pd.DataFrame
    arbres_finals: pd.DataFrame
    factor_alcada: float


def simular(rodal: Rodal, models: dict[str, ModelEspecie], anys: int,
            pertorbacio: Pertorbacio = SENSE_PERTORBACIONS) -> Simulacio:
    if not np.isfinite(rodal.altitud):
        raise ValueError("Cal l'altitud del rodal per simular el creixement.")
    factor = factor_alcada_local(rodal)
    arbres = rodal.arbres[["sp", "dbh", "f"]].copy().reset_index(drop=True)
    arbres["altitud"] = rodal.altitud
    es_pi = arbres["sp"].isin(models).to_numpy()
    pins_inicials = arbres.loc[es_pi]
    sp_dominant = (
        pins_inicials.assign(g=pins_inicials["dbh"] ** 2 * pins_inicials["f"]).groupby("sp")["g"].sum().idxmax()
    )
    # Un rodal que ja parteix per sobre del màxim no s'aclareix de cop: el seu límit és la densitat inicial.
    sdi_max = max(models[sp_dominant].sdi_max, _sdi_rodal(arbres["dbh"].to_numpy(), arbres["f"].to_numpy()))

    files = [_resum_any(0, arbres, rodal.habitat, factor, 0.0, 0.0) | {"increment_brut": 0.0,
                                                                       "morts_pertorbacio": 0.0}]
    inicials = estoc(arbres, rodal.habitat, factor)
    for any_ in range(1, anys + 1):
        comp = _competencia(arbres)
        increment = np.zeros(len(arbres))
        p_mort = np.zeros(len(arbres))
        for sp, model in models.items():
            idx = (arbres["sp"] == sp).to_numpy()
            if idx.any():
                increment[idx] = model.increment_diametral(comp.loc[idx])
                p_mort[idx] = model.mortalitat_anual(comp.loc[idx])
        abans = estoc(arbres, rodal.habitat, factor)
        crescuts = estoc(arbres.assign(dbh=arbres["dbh"] + increment), rodal.habitat, factor)
        increment_brut = float((crescuts["cat"] + crescuts["cst"]).sum() - (abans["cat"] + abans["cst"]).sum())
        f = arbres["f"].to_numpy()
        supervivents_fons = f * (1.0 - p_mort)
        morts_pert = supervivents_fons * (1.0 - pertorbacio.supervivencia(any_, comp, es_pi))
        morts = f - supervivents_fons + morts_pert
        arbres["f"] = f - morts
        arbres["dbh"] = arbres["dbh"] + increment
        arbres, retirats = _autoaclarida(arbres, es_pi, sdi_max)
        f_abans = abans["f"].to_numpy()
        f_despres = arbres.loc[abans.index, "f"].to_numpy()
        fraccio_morta = np.divide(f_abans - f_despres, f_abans, out=np.zeros_like(f_abans), where=f_abans > 0)
        c_mort = float(np.sum((abans["cat"] + abans["cst"]).to_numpy() * fraccio_morta))
        resum = _resum_any(any_, arbres, rodal.habitat, factor, float(morts[es_pi].sum() + retirats), c_mort)
        files.append(resum | {"increment_brut": increment_brut, "morts_pertorbacio": float(morts_pert.sum())})
    return Simulacio(
        serie=pd.DataFrame(files),
        arbres_inicials=inicials,
        arbres_finals=estoc(arbres, rodal.habitat, factor),
        factor_alcada=factor,
    )
