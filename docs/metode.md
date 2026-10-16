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

on $\Delta D_{10}$ és el creixement diametral en 10 anys (cm), $D$ el diàmetre normal (cm), $BAL$ l'àrea basimètrica dels arbres més gruixuts (m²/ha) i $AB$ la de la parcel·la; $P$ val 1 a la remesura IFN3→IFN4 i 0 a la IFN2→IFN3. $u_p$ és la constant aleatòria de la parcel·la, comuna a les dues remesures (qualitat de lloc), i $v_{pt}$ la de la parcel·la en un període concret: recull que una parcel·la no creix igual a les dues remesures per raons pròpies (clima local, danys, gestió). Sense $v_{pt}$, les observacions dels dos períodes d'una parcel·la es tractarien com si fossin més independents del que són i l'error estàndard de $\beta_6$ quedaria subestimat; el model niat s'ajusta clarament millor (raó de versemblances de 105 per al pi roig i 49 per al pi negre). Les desviacions típiques (REML) són $\sigma_u$ = 0,25 i $\sigma_v$ = 0,11 (pi roig) i 0,16 i 0,13 (pi negre). Les prediccions fan servir els efectes fixos del període més recent ($P = 1$, $u = v = 0$) i una correcció de biaix de la retransformació logarítmica (estimador de raó). Entre les dues remesures, el creixement va baixar un 4,1 % (pi roig, p = 0,01) i un 3,6 % (pi negre, p = 0,09).

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
- **Risc mitjà** (per defecte): cada any mor, a més, la fracció de pins que es va perdre de mitjana per pertorbacions (0,34 % o 0,42 % de l'AB segons l'espècie), igual a totes les mides. És determinista: representa l'esperança de moltes parcel·les, no una parcel·la concreta. Amb aquest escenari la mortalitat total del model és comparable a l'observada.
- **Episodi puntual**: l'any triat mor una fracció dels pins de totes les mides (foc, sequera, neu o vent, plaga). Per defecte, la severitat mitjana del tipus a les parcel·les afectades de l'espècie; si l'espècie en té menys de tres, la de totes dues espècies juntes (foc 47 %, sequera 61 %). La resta d'anys, només la mortalitat de fons.

### Autoaclarida

L'índex de densitat de Reineke, $SDI = N\,(D_q/25)^{1,605}$, no pot superar el percentil 99 observat a les parcel·les de pineda d'Andorra i de Lleida (a Lleida, l'estat inicial de cada remesura). Si el creixement el faria superar, moren els pins més prims fins a tornar-hi. Un rodal que ja parteix per sobre del màxim manté la seva densitat inicial com a límit.

### Aclarides

L'aclarida treu la fracció indicada de l'àrea basimètrica de pins: per baix (primer els més prims), per dalt (primer els més gruixuts) o sistemàtica (la mateixa proporció de cada mida). Després el rodal aclarit se simula amb els mateixos models; com que té menys competència, cada arbre hi creix més.

### Validació

- **Validació creuada per parcel·les** a Lleida (parcel·les de pineda): amb la mortalitat de fons i la pèrdua mitjana per pertorbacions, l'àrea basimètrica predita a la remesura (11–15 anys) té R² = 0,92–0,93 i un biaix d'entre −0,07 i −0,16 m²/ha. Només amb la mortalitat de fons, el biaix és de +0,9 a +1,1 m²/ha (el model no hi posa les pertorbacions). A les parcel·les que van perdre menys de 2 m²/ha per pertorbacions (394 de 444 de pi roig i 247 de 307 de pi negre), la mortalitat de fons sola dona R² = 0,94–0,96 i un biaix de +0,3 a +0,4 m²/ha. El nombre de morts de fons predit coincideix amb l'observat.
- **Estoc a Andorra**: el carboni aeri mitjà de cada hàbitat de pineda queda a menys d'un 3 % del publicat a l'informe de l'INF1-AND (taula 26, per hàbitat).
- **Embornal**: l'informe de l'INF1-AND estima el creixement amb les equacions IAVC de l'IFN2, que no tenen en compte la competència. A les parcel·les de Lleida amb més de 40 m²/ha aquestes equacions sobreestimen el creixement observat en volum entre 1,9 i 2,4 vegades, mentre que aquest model s'hi ajusta amb un error de menys del ±8 % en totes les classes de densitat. Per això l'embornal que dona l'app per a les pinedes denses d'Andorra és més baix que el de l'informe.

### Limitacions

- Les **pertorbacions** són una mitjana del Pirineu de Lleida entre 1989 i 2016, amb poques parcel·les cremades o afectades per sequera; no inclouen l'augment esperat amb el canvi climàtic. El risc mitjà és determinista i reparteix la mortalitat entre totes les mides; un episodi també mata la mateixa fracció de cada mida. A les dades, en canvi, la neu i les plagues maten més pins negres grossos: el model en subestima una mica la mortalitat.
- No hi ha **incorporació** de peus nous (regeneració), de manera que a llarg termini la captació de rodals oberts s'infravalora.
- El creixement està calibrat al **Pirineu de Lleida (1989–2016)**: el terme de període només recull la diferència mitjana entre les dues remesures; no inclou l'efecte del canvi climàtic posterior ni diferències de qualitat d'estació pròpies d'Andorra. Quan hi hagi la remesura de l'INF-AND es podrà calibrar amb dades andorranes.
- Les masses andorranes són més denses que les de calibratge; per sobre de l'àrea basimètrica màxima de Lleida el creixement s'extrapola i l'app ho avisa.
- Les **altres espècies** del rodal compten per a la competència, però no creixen ni s'inclouen en el carboni.
- El carboni de la **fusta extreta** es mostra a part: que torni o no a l'atmosfera depèn de l'ús que se'n faci.
- La **captació neta** és la variació del carboni dels pins vius: els arbres que moren en surten, tot i que en realitat queden com a fusta morta i es descomponen lentament. Com que els pins negres més gruixuts tenen més mortalitat (de fons, un 0,46 % anual a partir de 45 cm; amb les pertorbacions, un 1,2 % a les dades), una aclarida per dalt n'evita les pèrdues i surt afavorida en captació neta. L'estoc en peu al final, l'estructura del rodal i la biodiversitat lligada als arbres grossos també compten a l'hora de decidir.
