**En resum.** L'app simula cada pi de la parcel·la any a any: creix, pot morir i, si el rodal arriba a la densitat màxima, s'autoaclareix. Amb els diàmetres i les alçades de cada any calcula la biomassa i el carboni. Tot el model es basa en dades públiques.

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

Model lineal mixt (efecte aleatori de parcel·la) ajustat per espècie amb els pins vius de les parcel·les de l'IFN de les quatre comarques pirinenques de Lleida remesurades sense tallades. Les dues remesures (IFN2→IFN3 i IFN3→IFN4) s'ajusten juntes amb un terme de període:

$$\ln(\Delta D_{10} + 1) = \beta_0 + \beta_1 \ln D + \beta_2 D + \beta_3 \frac{BAL}{\ln(D+1)} + \beta_4 \ln AB + \beta_5\,\text{altitud} + \beta_6\,P$$

on $\Delta D_{10}$ és el creixement diametral en 10 anys (cm), $D$ el diàmetre normal (cm), $BAL$ l'àrea basimètrica dels arbres més gruixuts (m²/ha) i $AB$ la de la parcel·la; $P$ val 1 a la remesura IFN3→IFN4 i 0 a la IFN2→IFN3. Les prediccions fan servir els efectes fixos del període més recent ($P = 1$) i una correcció de biaix de la retransformació logarítmica (estimador de raó). Entre les dues remesures, el creixement va baixar un 2,7 % (pi roig) i un 3,9 % (pi negre).

### Mortalitat

Model binomial amb enllaç log-log complementari i la durada de l'interval com a *offset*, que dona directament una probabilitat anual:

$$p_{anual} = 1 - \exp\!\big(-\exp(\gamma_0 + \gamma_1/D + \gamma_2 D + \gamma_3 BAL + \gamma_4 P)\big)$$

Recull la mortalitat dels arbres petits i dominats i la dels grans i vells. També es fa servir el període més recent: el risc de mort hi és 1,10 vegades el de la remesura anterior per al pi roig (no significatiu) i 0,84 vegades per al pi negre.

### Autoaclarida

L'índex de densitat de Reineke, $SDI = N\,(D_q/25)^{1,605}$, no pot superar el percentil 99 observat a les parcel·les de pineda d'Andorra i de Lleida (a Lleida, l'estat inicial de cada remesura). Si el creixement el faria superar, moren els pins més prims fins a tornar-hi. Un rodal que ja parteix per sobre del màxim manté la seva densitat inicial com a límit.

### Aclarides

L'aclarida treu la fracció indicada de l'àrea basimètrica de pins: per baix (primer els més prims), per dalt (primer els més gruixuts) o sistemàtica (la mateixa proporció de cada mida). Després el rodal aclarit se simula amb els mateixos models; com que té menys competència, cada arbre hi creix més.

### Validació

- **Validació creuada per parcel·les** a Lleida: l'àrea basimètrica predita a la remesura (11–15 anys) té R² = 0,91–0,94 i un biaix d'entre −0,11 i −0,16 m²/ha. El nombre total de morts predit coincideix amb l'observat.
- **Estoc a Andorra**: el carboni aeri mitjà de cada hàbitat de pineda queda a menys d'un 3 % del publicat a l'informe de l'INF1-AND (taula 26, per hàbitat).
- **Embornal**: l'informe de l'INF1-AND estima el creixement amb les equacions IAVC de l'IFN2, que no tenen en compte la competència. A les parcel·les de Lleida amb més de 40 m²/ha aquestes equacions sobreestimen el creixement observat en volum entre 1,9 i 2,4 vegades, mentre que aquest model s'hi ajusta amb un error de menys del ±8 % en totes les classes de densitat. Per això l'embornal que dona l'app per a les pinedes denses d'Andorra és més baix que el de l'informe.

### Limitacions

- No hi ha **incorporació** de peus nous (regeneració), de manera que a llarg termini la captació de rodals oberts s'infravalora.
- El creixement està calibrat al **Pirineu de Lleida (1989–2016)**: el terme de període només recull la diferència mitjana entre les dues remesures; no inclou l'efecte del canvi climàtic posterior ni diferències de qualitat d'estació pròpies d'Andorra. Quan hi hagi la remesura de l'INF-AND es podrà calibrar amb dades andorranes.
- Les masses andorranes són més denses que les de calibratge; per sobre de l'àrea basimètrica màxima de Lleida el creixement s'extrapola i l'app ho avisa.
- Les **altres espècies** del rodal compten per a la competència, però no creixen ni s'inclouen en el carboni.
- El carboni de la **fusta extreta** es mostra a part: que torni o no a l'atmosfera depèn de l'ús que se'n faci.
- La **captació neta** és la variació del carboni dels pins vius: els arbres que moren en surten, tot i que en realitat queden com a fusta morta i es descomponen lentament. Com que a les dades els pins més gruixuts tenen més mortalitat (sobretot el pi negre, d'un 1,2 % anual a partir de 45 cm), una aclarida per dalt n'evita les pèrdues i surt afavorida en captació neta. L'estoc en peu al final, l'estructura del rodal i la biodiversitat lligada als arbres grossos també compten a l'hora de decidir.
