### Esquema

L'app simula cada pi de la parcel·la any a any: creix, pot morir i, si el rodal arriba a la densitat màxima, s'autoaclareix. Amb els diàmetres i les alçades de cada any calcula la biomassa i el carboni. Tot el model es basa en dades públiques.

| Peça | Font |
|---|---|
| Rodals de partida | 1r Inventari Nacional Forestal d'Andorra, BD-INF1-AND (ARI i CTFC, 2022), 194 parcel·les |
| Creixement i mortalitat | Inventari Forestal Nacional, IFN2 (1989–90) i IFN3 (2000–01) de Lleida, 15.536 pins remesurats |
| Alçada–diàmetre | Equacions pròpies de l'INF-AND per a pi roig i pi negre, en hàbitat xeròfil i mesòfil |
| Biomassa aèria | IEFC Catalunya, BAT = a·Dn^b·Ht^c |
| Biomassa subterrània | INIA Espanya, BST = a·Dn^b |
| Carboni | 0,5 t C per tona de biomassa; CO₂ = C · 44/12 |

### Creixement diametral

Model lineal mixt (efecte aleatori de parcel·la) ajustat per espècie amb els pins vius de les parcel·les de l'IFN de Lleida remesurades sense tallades:

$$\ln(\Delta D_{10} + 1) = \beta_0 + \beta_1 \ln D + \beta_2 D + \beta_3 \frac{BAL}{\ln(D+1)} + \beta_4 \ln AB + \beta_5\,\text{altitud}$$

on $\Delta D_{10}$ és el creixement diametral en 10 anys (cm), $D$ el diàmetre normal (cm), $BAL$ l'àrea basimètrica dels arbres més gruixuts (m²/ha) i $AB$ la de la parcel·la. Les prediccions fan servir els efectes fixos i una correcció de biaix de la retransformació logarítmica (estimador de raó).

### Mortalitat

Model binomial amb enllaç log-log complementari i la durada de l'interval com a *offset*, que dona directament una probabilitat anual:

$$p_{anual} = 1 - \exp\!\big(-\exp(\gamma_0 + \gamma_1/D + \gamma_2 D + \gamma_3 BAL)\big)$$

Recull la mortalitat dels arbres petits i dominats i la dels grans i vells.

### Autoaclarida

L'índex de densitat de Reineke, $SDI = N\,(D_q/25)^{1,605}$, no pot superar el percentil 99 observat a les parcel·les de pinar d'Andorra i de Lleida. Si el creixement el faria superar, moren els pins més prims fins a tornar-hi. Un rodal que ja parteix per sobre del màxim manté la seva densitat inicial com a límit.

### Aclarides

L'aclarida treu la fracció indicada de l'àrea basimètrica de pins: per baix (primer els més prims), per dalt (primer els més gruixuts) o sistemàtica (la mateixa proporció de cada mida). Després el rodal aclarit se simula amb els mateixos models; com que té menys competència, cada arbre hi creix més.

### Validació

- **Validació creuada per parcel·les** a Lleida: l'àrea basimètrica predita al cap d'11 anys té R² = 0,94–0,95 i un biaix d'entre −0,15 i −0,22 m²/ha. El nombre total de morts predit coincideix amb l'observat.
- **Estoc a Andorra**: el carboni aeri mitjà de cada estrat de pinar queda a menys d'un 10 % del publicat a l'informe de l'INF1-AND.
- **Embornal**: l'informe de l'INF1-AND estima el creixement amb les equacions IAVC de l'IFN2, que no tenen en compte la competència. A les parcel·les de Lleida amb més de 40 m²/ha aquestes equacions sobreestimen el creixement observat en volum entre 2 i 2,4 vegades, mentre que aquest model s'hi ajusta amb un error de menys del ±8 % en totes les classes de densitat. Per això l'embornal que dona l'app per a les pinedes denses d'Andorra és més baix que el de l'informe.

### Limitacions

- No hi ha **incorporació** de peus nous (regeneració), de manera que a llarg termini la captació de rodals oberts s'infravalora.
- El creixement està calibrat a **Lleida (1990–2001)**: no inclou l'efecte del canvi climàtic posterior ni diferències de qualitat d'estació pròpies d'Andorra. Quan hi hagi la remesura de l'INF-AND es podrà calibrar amb dades andorranes.
- Les masses andorranes són més denses que les de calibratge; per sobre de l'àrea basimètrica màxima de Lleida el creixement s'extrapola i l'app ho avisa.
- Les **altres espècies** del rodal compten per a la competència, però no creixen ni s'inclouen en el carboni.
- El carboni de la **fusta extreta** es mostra a part: que torni o no a l'atmosfera depèn de l'ús que se'n faci.
- La **captació neta** és la variació del carboni dels pins vius: els arbres que moren en surten, tot i que en realitat queden com a fusta morta i es descomponen lentament. Com que a les dades els pins més gruixuts tenen més mortalitat (sobretot el pi negre, d'un 1,3 % anual a partir de 45 cm), una aclarida per dalt n'evita les pèrdues i surt afavorida en captació neta. L'estoc en peu al final, l'estructura del rodal i la biodiversitat lligada als arbres grossos també compten a l'hora de decidir.
