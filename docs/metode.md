**En resum.** L'app simula cada pi de la parcel·la any a any: creix, pot morir (per competència o per una pertorbació, segons l'escenari) i, si el rodal arriba a la densitat màxima, s'autoaclareix. Amb els diàmetres i les alçades de cada any calcula la biomassa i el carboni. Tot el model es basa en dades públiques.

### Fonts

| Peça | Font |
|---|---|
| Rodals de partida | 1r Inventari Nacional Forestal d'Andorra, BD-INF1-AND (ARI i CTFC, 2022), 194 parcel·les |
| Creixement i mortalitat | Inventari Forestal Nacional, IFN2 (1989–90), IFN3 (2000–01) i IFN4 (2014–16) del Pallars Sobirà, l'Alt Urgell, l'Alta Ribagorça i la Cerdanya (Lleida): dues remesures, 19.603 registres de pins en 784 parcel·les |
| Alçada–diàmetre | Equacions pròpies de l'INF-AND per a pi roig i pi negre, en hàbitat xeròfil i mesòfil |
| Biomassa aèria | IEFC Catalunya, BAT = a·Dn^b·Ht^c |
| Biomassa subterrània | INIA Espanya, BST = a·Dn^b |
| Carboni | 0,5 t C per tona de biomassa; CO₂ = C · 44/12 |

### Creixement diametral

Model lineal mixt (efectes aleatoris niats de parcel·la i de parcel·la×període) ajustat per espècie amb els pins vius de les parcel·les de l'IFN de les quatre comarques pirinenques de Lleida remesurades sense tallades. Les dues remesures (IFN2→IFN3 i IFN3→IFN4) s'ajusten juntes amb un terme de període:

$$\ln(\Delta D_{10} + 1) = \beta_0 + \beta_1 \ln D + \beta_2 D + \beta_3 \frac{BAL}{\ln(D+1)} + \beta_4 \ln AB + \beta_5\,\text{altitud} + \beta_6\,P + u_p + v_{pt} + e$$

on $\Delta D_{10}$ és el creixement diametral en 10 anys (cm), $D$ el diàmetre normal (cm), $BAL$ l'àrea basimètrica dels arbres més gruixuts (m²/ha) i $AB$ la de la parcel·la; $P$ val 1 a la remesura IFN3→IFN4 i 0 a la IFN2→IFN3. $u_p$ és la constant aleatòria de la parcel·la, comuna a les dues remesures (qualitat de lloc), i $v_{pt}$ la de la parcel·la en un període concret: recull que una parcel·la no creix igual a les dues remesures per raons pròpies (clima local, danys, gestió). Sense $v_{pt}$, les observacions dels dos períodes d'una parcel·la es tractarien com si fossin més independents del que són i l'error estàndard de $\beta_6$ quedaria subestimat; el model niat s'ajusta clarament millor (raó de versemblances de 105 per al pi roig i 49 per al pi negre). Les desviacions típiques (REML) són $\sigma_u$ = 0,25 i $\sigma_v$ = 0,11 (pi roig) i 0,16 i 0,13 (pi negre). Les prediccions fan servir els efectes fixos del període més recent ($P = 1$, $u = v = 0$) i una correcció de biaix de la retransformació logarítmica (estimador de raó). Entre les dues remesures, el creixement del pi roig va baixar un 4,1 % (p = 0,01); el del pi negre, un 3,6 %, però la diferència no és significativa (p = 0,09).

### Mortalitat de fons

L'IFN anota, per a cada arbre mort, l'**agent causant** (camp *Agente* de l'IFN3 i l'IFN4; documentador de l'IFN4, annex 8). Els agents s'agrupen així:

| Grup | Agents (codi IFN) |
|---|---|
| **Fons** | competència o dominància (316), causa desconeguda (200), sense dany (100) o sense agent anotat |
| Foc | incendi (421), llamp (413) |
| Sequera | sequera (412), gelada (414) |
| Neu o vent | neu (410), vent (411), calamarsa (415) |
| Biòtic | biòtic genèric (300), fongs (310), insectes (311), vesc (312), epífits (313), fauna (314), bestiar (315) |
| Altres | maquinària (320), aprofitament (321), persones (322), abiòtic genèric (400), despreniments (422), erosió (423) |

Només un 38 % dels pins morts a les parcel·les de calibratge són de **fons** (569 de 1.484); la resta són pertorbacions, sobretot biòtiques i de neu o vent. El model de mortalitat s'ajusta **només amb les morts de fons**: els pins morts per pertorbacions es treuen de l'ajust, és a dir, es tracten com a censurats (no se sap quan van morir, però sí que fins aleshores no havien mort de fons). És la simplificació habitual dels riscos competitius; com que són pocs pins (un 5 %), el tractament de la seva exposició gairebé no canvia els coeficients.

Model binomial amb enllaç log-log complementari i la durada de l'interval com a *offset*, que dona directament una probabilitat anual:

$$p_{anual} = 1 - \exp\!\big(-\exp(\gamma_0 + \gamma_1/D + \gamma_2 D + \gamma_3 BAL + \gamma_4 P)\big)$$

Recull la mortalitat dels arbres petits i dominats i la dels grans i vells. També es fa servir el període més recent: el risc de mort de fons hi és 1,29 vegades el de la remesura anterior per al pi roig i 1,01 vegades per al pi negre (no significatiu). Sense les pertorbacions, la capacitat de discriminació del model millora (AUC de 0,60 a 0,67 per al pi roig i de 0,59 a 0,62 per al pi negre).

### Pertorbacions

Les pertorbacions s'estimen a escala de parcel·la, amb tots dos períodes junts i totes les parcel·les remesurades on hi ha l'espècie (són poc freqüents i calen totes; la fracció es mesura sobre l'àrea basimètrica de pins de l'espècie, de manera que no depèn de la composició). Per a cada tipus:

- **Pèrdua anual d'AB**: AB inicial dels pins morts per aquesta causa dividida per la suma d'AB inicial × anys de totes les parcel·les. Compta totes les morts, també les aïllades. En total, un 0,34 % anual de l'AB de pi roig i un 0,42 % de la de pi negre (la mortalitat de fons observada és d'un 0,21 % i un 0,26 %).
- **Parcel·les afectades**: on la causa mata almenys el 10 % de l'AB inicial de pins. La probabilitat anual d'episodi és $1 - e^{-\lambda}$, amb $\lambda$ = parcel·les afectades / anys-parcel·la observats.
- **Severitat**: fracció mitjana (i mediana) de l'AB morta a les parcel·les afectades.

| | Pi roig: pèrdua anual | parcel·les afectades | severitat mitjana | Pi negre: pèrdua anual | parcel·les afectades | severitat mitjana |
|---|---|---|---|---|---|---|
| Foc | 0,030 % | 3 | 57 % | 0,005 % | 1 | (16 %) |
| Sequera | 0,023 % | 5 | 68 % | 0,001 % | 1 | (26 %) |
| Neu o vent | 0,096 % | 20 | 18 % | 0,158 % | 23 | 22 % |
| Biòtic | 0,164 % | 36 | 27 % | 0,197 % | 33 | 22 % |
| Altres | 0,029 % | 11 | 39 % | 0,062 % | 13 | 35 % |

Els valors es desen al bloc `pertorbacions` de `data/calibracio.json`. L'app ofereix tres escenaris, que s'apliquen igual al rodal amb i sense aclarida:

- **Sense pertorbacions**: només la mortalitat de fons. És una projecció optimista.
- **Risc mitjà** (per defecte): cada any es perd, a més, la fracció d'AB que es va perdre de mitjana per pertorbacions (0,34 % al pi roig i 0,42 % al pi negre). Al pi roig, igual a totes les mides; al pi negre, amb un risc que creix amb el diàmetre (vegeu més avall). És determinista: representa l'esperança de moltes parcel·les, no una parcel·la concreta. Amb aquest escenari la mortalitat total del model és comparable a l'observada.
- **Episodi puntual**: l'any triat mor una fracció dels pins de totes les mides (foc, sequera, neu o vent, plaga). Per defecte, la severitat mitjana del tipus a les parcel·les afectades de l'espècie; si l'espècie en té menys de tres, la de totes dues espècies juntes (foc 47 %, sequera 61 %). La resta d'anys, només la mortalitat de fons.

#### Risc segons la mida

A les dades, els pins negres gruixuts moren més per pertorbacions que els prims. Per a cada espècie s'ha ajustat un model binomial cloglog de la probabilitat de morir per una pertorbació, amb la durada de l'interval com a *offset*, i s'ha comparat amb el model amb només constant, que equival al risc uniforme. Els pins morts de fons són el risc competidor i es tracten com a censurats, igual que es fa amb els morts per pertorbacions a la mortalitat de fons (comptar-los com a supervivents dona el mateix resultat). Com que les morts per pertorbacions van per episodis, agrupades per parcel·les, a més de la raó de versemblances (que suposa arbres independents) s'ha fet una prova de Wald amb errors robustos per parcel·la i s'ha mirat l'efecte del diàmetre dins de cada parcel·la afectada.

| | Pi negre | Pi roig |
|---|---|---|
| Raó de versemblances del diàmetre (1 gl) | 23,4 (p < 0,001) | 4,9 (p = 0,03) |
| ΔAIC respecte del risc uniforme | −21,4 | −2,9 |
| Wald robust per parcel·la | p < 0,001 | p = 0,39 |
| Efecte dins de la parcel·la, per cm | +0,012 (p = 0,02) | −0,028 (p < 0,001) |
| Per causa (robust) | biòtic +0,026 (p < 0,001); neu o vent +0,008 (p = 0,49) | cap de significatiu |

Al pi negre el risc creix un 1,9 % per cada centímetre de diàmetre (×2,1 entre 10 i 50 cm), sobretot per les causes biòtiques; amb la neu i el vent la tendència va en el mateix sentit, però no és significativa. Afegir-hi el BAL, ln D o D² no millora clarament el model (ΔAIC ≤ 2 o pitjor). Al pi roig no hi ha un patró clar: el diàmetre sol no és significatiu amb errors robustos, dins de les parcel·les afectades moren més els pins prims i la relació té forma de U (més risc als prims i, pel foc, als gruixuts). Els models amb D² i BAL milloren l'AIC (ΔAIC −10 a −30), però empitjoren la validació creuada de les parcel·les (R² de l'increment d'AB de 0,177 a 0,168) i no milloren el calibratge dels pins gruixuts. Al pi roig es manté el risc uniforme.

Mortalitat anual per pertorbacions per classes diamètriques (% d'arbres per any, validació creuada per parcel·les):

| Pi negre, D (cm) | 7,5–15 | 15–25 | 25–35 | 35–45 | > 45 |
|---|---|---|---|---|---|
| Observada | 0,32 | 0,35 | 0,40 | 0,56 | 0,70 |
| Risc uniforme | 0,41 | 0,41 | 0,41 | 0,41 | 0,41 |
| Risc segons la mida | 0,29 | 0,34 | 0,40 | 0,48 | 0,64 |
| Total (fons + pertorbacions) observada | 0,67 | 0,57 | 0,60 | 0,84 | 1,11 |
| Total amb el risc uniforme | 0,74 | 0,67 | 0,65 | 0,64 | 0,81 |
| Total amb el risc segons la mida | 0,62 | 0,59 | 0,64 | 0,71 | 1,04 |

La validació creuada de les parcel·les de pineda de pi negre (307) millora una mica: AB final R² 0,934 → 0,935, RMSE 3,57 → 3,53 m²/ha, biaix −0,16 → −0,12 m²/ha; increment d'AB R² 0,32 → 0,33. Per això el risc segons la mida s'aplica al pi negre i no al pi roig. La constant es recalibra perquè, a les parcel·les de calibratge, la pèrdua d'AB esperada sigui la mateixa taxa observada (0,42 % anual): el model canvia el repartiment entre mides, no el total. En un rodal amb més pins gruixuts que la mitjana de calibratge la pèrdua és més alta, i al revés. Els coeficients són al bloc `pertorbacions.mida` de `calibracio.json`.

L'episodi puntual es manté igual per a totes les mides: la severitat és la fracció d'AB que tria l'usuari, i dins de les parcel·les afectades l'efecte de la mida és petit.

### Autoaclarida

L'índex de densitat de Reineke, $SDI = N\,(D_q/25)^{1,605}$, no pot superar el percentil 99 observat a les parcel·les de pineda d'Andorra i de Lleida (a Lleida, l'estat inicial de cada remesura). Si el creixement el faria superar, moren els pins més prims fins a tornar-hi. Un rodal que ja parteix per sobre del màxim manté la seva densitat inicial com a límit.

### Aclarides

L'aclarida treu la fracció indicada de l'àrea basimètrica de pins: per baix (primer els més prims), per dalt (primer els més gruixuts) o sistemàtica (la mateixa proporció de cada mida). Després el rodal aclarit se simula amb els mateixos models; com que té menys competència, cada arbre hi creix més.

### Validació

- **Validació creuada per parcel·les** a Lleida (parcel·les de pineda): amb la mortalitat de fons i la pèrdua mitjana per pertorbacions, l'àrea basimètrica predita a la remesura (11–15 anys) té R² = 0,92–0,93 i un biaix d'entre −0,07 i −0,16 m²/ha. Només amb la mortalitat de fons, el biaix és de +0,9 a +1,1 m²/ha (el model no hi posa les pertorbacions). A les parcel·les que van perdre menys de 2 m²/ha per pertorbacions (394 de 444 de pi roig i 247 de 307 de pi negre), la mortalitat de fons sola dona R² = 0,94–0,96 i un biaix de +0,3 a +0,4 m²/ha. El nombre de morts de fons predit coincideix amb l'observat.
- **Estoc a Andorra**: el carboni aeri mitjà de cada hàbitat de pineda queda a menys d'un 3 % del publicat a l'informe de l'INF1-AND (taula 26, per hàbitat).
- **Embornal**: l'informe de l'INF1-AND estima el creixement amb les equacions IAVC de l'IFN2, que no tenen en compte la competència. A les parcel·les de Lleida amb més de 40 m²/ha aquestes equacions sobreestimen el creixement observat en volum entre 1,9 i 2,4 vegades, mentre que aquest model s'hi ajusta amb un error de menys del ±8 % en totes les classes de densitat. Per això l'embornal que dona l'app per a les pinedes denses d'Andorra és més baix que el de l'informe.

### Limitacions

- Les **pertorbacions** són una mitjana del Pirineu de Lleida entre 1989 i 2016, amb poques parcel·les cremades o afectades per sequera; no inclouen l'augment esperat amb el canvi climàtic. El risc mitjà és determinista. Al pi negre el risc creix amb el diàmetre; al pi roig, i en els episodis puntuals, es reparteix igual entre totes les mides.
- No hi ha **incorporació** de peus nous (regeneració), de manera que a llarg termini la captació de rodals oberts s'infravalora.
- El creixement està calibrat al **Pirineu de Lleida (1989–2016)**: el terme de període només recull la diferència mitjana entre les dues remesures; no inclou l'efecte del canvi climàtic posterior ni diferències de qualitat d'estació pròpies d'Andorra. Quan hi hagi la remesura de l'INF-AND es podrà calibrar amb dades andorranes.
- Les masses andorranes són més denses que les de calibratge; per sobre de l'àrea basimètrica màxima de Lleida el creixement s'extrapola i l'app ho avisa.
- Les **altres espècies** del rodal compten per a la competència, però no creixen ni s'inclouen en el carboni.
- El carboni de la **fusta extreta** es mostra a part: que torni o no a l'atmosfera depèn de l'ús que se'n faci.
- La **captació neta** és la variació del carboni dels pins vius: els arbres que moren en surten, tot i que en realitat queden com a fusta morta i es descomponen lentament. Com que els pins negres més gruixuts tenen més mortalitat (de fons, un 0,46 % anual a partir de 45 cm; amb les pertorbacions, un 1,1 % a les dades i un 1,0 % al model), una aclarida per dalt n'evita les pèrdues i surt afavorida en captació neta. L'estoc en peu al final, l'estructura del rodal i la biodiversitat lligada als arbres grossos també compten a l'hora de decidir.
