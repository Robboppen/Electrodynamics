<h1 align="center">Resumen del curso — Electromagnetismo Intermedio (LFIS322)</h1>

<p align="center"><em>Formulario comentado: las ecuaciones centrales de cada unidad, de dónde salen, cuándo valen y qué significan.</em></p>

<p align="center"><a href="README.md">← volver al README</a></p>

Unidades SI y notación de Griffiths (4.ª ed.) y Jackson (3.ª ed.). Vector separación: $\boldsymbol{\mathfrak R}\equiv\mathbf r-\mathbf r'$, con $\mathfrak R=|\boldsymbol{\mathfrak R}|$.

---

## Índice

| # | Unidad | Clases (programa) | Simulaciones |
|---|---|---|---|
| 1 | [Herramientas matemáticas](#s1) | 1–3 | — |
| 2 | [Electrostática](#s2) | 3–4 | — |
| 3 | [Energía electrostática, conductores y capacidad](#s3) | 5 | — |
| 4 | [Problemas de contorno: unicidad, Green, imágenes, separación de variables](#s4) | 6–9 | — |
| 5 | [Expansión multipolar](#s5) | 12–13 | — |
| 6 | [Magnetostática](#s6) | 13–15 | — |
| 7 | [Campos en la materia: dieléctricos y materiales magnéticos](#s7) | 16–18 | — |
| 8 | [Inducción y ecuaciones de Maxwell](#s8) | 18–20 | — |
| 9 | [Leyes de conservación: energía, momento y momento angular](#s9) | 20, 23–24 | — |
| 10 | [Potenciales, gauge y potenciales retardados](#s10) | 25 | [Guía 6 · P13–P14](EI_Simulaciones/Guia_06_Radiacion_Covarianza/03_maxwell_accion_gauge.ipynb) |
| 11 | [Liénard–Wiechert: campos de una carga en movimiento](#s11) | 26 | [Guía 5 · P1–P3](EI_Simulaciones/Guia_05_Radiacion_Lienard_Wiechert/01_tiempo_retardado_lienard_wiechert.ipynb) |
| 12 | [Radiación](#s12) | 27–28 | [Guía 5](EI_Simulaciones/Guia_05_Radiacion_Lienard_Wiechert/) · [Guía 6 · P1–P4](EI_Simulaciones/Guia_06_Radiacion_Covarianza/01_radiacion.ipynb) |
| 13 | [Relatividad especial y formulación covariante](#s13) | Guía 6 (2026) | [Guía 6 · P5–P15](EI_Simulaciones/Guia_06_Radiacion_Covarianza/) |

---

<a id="s1"></a>
## 1. Herramientas matemáticas

**Notación indicial.** Índices repetidos se suman (convención de Einstein). Los dos símbolos fundamentales son la delta de Kronecker $\delta_{ij}$ y el símbolo de Levi-Civita $\epsilon_{ijk}$ (totalmente antisimétrico, $\epsilon_{123}=1$):

```math
\mathbf A\cdot\mathbf B=A_iB_i,\qquad (\mathbf A\times\mathbf B)_i=\epsilon_{ijk}A_jB_k,\qquad
\boxed{\epsilon_{ijk}\,\epsilon_{imn}=\delta_{jm}\delta_{kn}-\delta_{jn}\delta_{km}}
```

Con esta identidad salen todas las identidades vectoriales; por ejemplo, $\nabla\times(\nabla\times\mathbf A)=\nabla(\nabla\cdot\mathbf A)-\nabla^2\mathbf A$. Dos identidades estructurales del curso:

```math
\nabla\cdot(\nabla\times\mathbf A)=0,\qquad \nabla\times(\nabla f)=0 .
```

Son la razón por la que $\mathbf B=\nabla\times\mathbf A$ y $\mathbf E=-\nabla V$ resuelven automáticamente dos de las ecuaciones de Maxwell.

**Teoremas integrales.**

```math
\int_{\mathcal V}(\nabla\cdot\mathbf F)\,d\tau=\oint_{\mathcal S}\mathbf F\cdot d\mathbf a\quad\text{(Gauss)},\qquad
\int_{\mathcal S}(\nabla\times\mathbf F)\cdot d\mathbf a=\oint_{\mathcal C}\mathbf F\cdot d\mathbf l\quad\text{(Stokes)}.
```

**Delta de Dirac.** Representa cargas puntuales como distribuciones, $\rho(\mathbf r)=q\,\delta^3(\mathbf r-\mathbf r_0)$, con $\int f(\mathbf r)\delta^3(\mathbf r-\mathbf a)d\tau=f(\mathbf a)$. El resultado clave, que resuelve la paradoja de que $\nabla\cdot(\hat{\mathbf r}/r^2)$ "es cero" pero su flujo no:

```math
\boxed{\nabla\cdot\frac{\hat{\mathbf r}}{r^2}=4\pi\,\delta^3(\mathbf r)\qquad\Longleftrightarrow\qquad \nabla^2\frac1r=-4\pi\,\delta^3(\mathbf r)}
```

**Teorema de Helmholtz.** Un campo que se anula en el infinito queda determinado por su divergencia y su rotor. Por eso la electrostática y la magnetostática se plantean como "fijar $\nabla\cdot$ y $\nabla\times$".

**Coordenadas curvilíneas ortogonales** $(u_1,u_2,u_3)$ con factores de escala $h_i$ (esféricas: $h_r=1$, $h_\theta=r$, $h_\phi=r\sin\theta$):

```math
\nabla f=\sum_i\frac{1}{h_i}\frac{\partial f}{\partial u_i}\hat{\mathbf e}_i,\qquad
\nabla\cdot\mathbf F=\frac{1}{h_1h_2h_3}\sum_i\frac{\partial}{\partial u_i}\!\left(\frac{h_1h_2h_3}{h_i}F_i\right),\qquad
\nabla^2 f=\frac{1}{h_1h_2h_3}\sum_i\frac{\partial}{\partial u_i}\!\left(\frac{h_1h_2h_3}{h_i^2}\frac{\partial f}{\partial u_i}\right).
```

*Material:* `Material electro, guías, pruebas/Electro 2026/Notacion_indicial_y_curvilinear_coords.pdf`.

---

<a id="s2"></a>
## 2. Electrostática

**Ley de Coulomb y superposición.** Todo parte de aquí:

```math
\mathbf E(\mathbf r)=\frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r')\,\hat{\boldsymbol{\mathfrak R}}}{\mathfrak R^2}\,d\tau' .
```

Aplicando $\nabla\cdot$ y $\nabla\times$ (y usando la delta de la sección 1) se obtienen las **ecuaciones de la electrostática**:

```math
\boxed{\nabla\cdot\mathbf E=\frac{\rho}{\varepsilon_0},\qquad \nabla\times\mathbf E=0}
\qquad\Longrightarrow\qquad
\oint\mathbf E\cdot d\mathbf a=\frac{Q_{\rm enc}}{\varepsilon_0}\quad\text{(ley de Gauss)} .
```

**Potencial.** Como $\nabla\times\mathbf E=0$, existe $V$ con $\mathbf E=-\nabla V$:

```math
V(\mathbf r)=\frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r')}{\mathfrak R}\,d\tau',\qquad
\boxed{\nabla^2V=-\frac{\rho}{\varepsilon_0}}\ \ \text{(Poisson)},\qquad \nabla^2V=0\ \ \text{(Laplace, sin cargas)} .
```

**Condiciones de borde** sobre una superficie con densidad $\sigma$ (normal $\hat{\mathbf n}$ de "abajo" hacia "arriba"):

```math
E^\perp_{\rm arriba}-E^\perp_{\rm abajo}=\frac{\sigma}{\varepsilon_0},\qquad \mathbf E^\parallel_{\rm arriba}=\mathbf E^\parallel_{\rm abajo},\qquad V\ \text{continuo},\qquad
\frac{\partial V_{\rm arriba}}{\partial n}-\frac{\partial V_{\rm abajo}}{\partial n}=-\frac{\sigma}{\varepsilon_0}.
```

**Conductores (equilibrio).** $\mathbf E=0$ adentro; $\rho=0$ adentro (toda la carga en la superficie); $V$ constante en todo el conductor; justo afuera, $\mathbf E=(\sigma/\varepsilon_0)\hat{\mathbf n}$, perpendicular a la superficie.

---

<a id="s3"></a>
## 3. Energía electrostática, conductores y capacidad

El trabajo para ensamblar una distribución de carga es

```math
\boxed{W=\frac12\int\rho V\,d\tau=\frac{\varepsilon_0}{2}\int_{\text{todo el espacio}}E^2\,d\tau},\qquad u_E=\frac{\varepsilon_0}{2}E^2 .
```

**Comentarios.**
- La versión discreta $W=\frac12\sum_iq_iV(\mathbf r_i)$ **excluye** la autoenergía. La versión con $E^2$ la incluye, y para cargas puntuales diverge: es el primer indicio del problema de la autoenergía clásica, que reaparece en la reacción de radiación (sección 12).
- La energía **no** es lineal en los campos: $W_{\rm total}\neq W_1+W_2$, porque aparece el término cruzado $\varepsilon_0\int\mathbf E_1\cdot\mathbf E_2\,d\tau$.

**Fuerza sobre un conductor** (presión electrostática hacia afuera) y fuerza por **desplazamiento virtual**:

```math
f=\frac{\sigma^2}{2\varepsilon_0},\qquad F_x=-\left(\frac{\partial W}{\partial x}\right)_{Q}=+\left(\frac{\partial W}{\partial x}\right)_{V}.
```

**Capacidad.** Para dos conductores con $\pm Q$:

```math
C=\frac{Q}{\Delta V},\qquad W=\frac12C(\Delta V)^2=\frac{Q^2}{2C};\qquad
C_{\rm placas}=\frac{\varepsilon_0A}{d},\qquad C_{\rm esferas}=4\pi\varepsilon_0\frac{ab}{b-a}.
```

Con varios conductores se generaliza a una matriz de coeficientes: $Q_i=\sum_jC_{ij}V_j$.

---

<a id="s4"></a>
## 4. Problemas de contorno: unicidad, Green, imágenes y separación de variables

**Teoremas de unicidad.** La solución de Poisson en un volumen $\mathcal V$ es única si se da $V$ en la frontera (**Dirichlet**), o $\partial V/\partial n$ (**Neumann**, salvo una constante). Para conductores, basta dar la carga total de cada uno. *Consecuencia práctica:* cualquier método que produzca una solución que cumpla las condiciones de borde (imágenes, adivinanza, separación) da **la** solución.

**Teorema de Green y funciones de Green.** Con $\nabla'^2G(\mathbf r,\mathbf r')=-4\pi\delta^3(\mathbf r-\mathbf r')$, $G=1/\mathfrak R+F$ y $\nabla'^2F=0$ (convención de Jackson):

```math
\boxed{V(\mathbf r)=\frac{1}{4\pi\varepsilon_0}\int_{\mathcal V}\rho(\mathbf r')\,G_D(\mathbf r,\mathbf r')\,d\tau'-\frac{1}{4\pi}\oint_{\mathcal S}V(\mathbf r')\frac{\partial G_D}{\partial n'}\,da'}\qquad\big(G_D=0\ \text{en}\ \mathcal S\big)
```

```math
V(\mathbf r)=\langle V\rangle_{\mathcal S}+\frac{1}{4\pi\varepsilon_0}\int_{\mathcal V}\rho\,G_N\,d\tau'+\frac{1}{4\pi}\oint_{\mathcal S}\frac{\partial V}{\partial n'}\,G_N\,da'\qquad\Big(\frac{\partial G_N}{\partial n'}=-\frac{4\pi}{S}\ \text{en}\ \mathcal S\Big)
```

$G$ es el potencial de una carga puntual unitaria *en presencia de la frontera*. Encontrar $G$ una vez resuelve el problema para **cualquier** $\rho$ y cualquier dato de borde. **Teorema del valor medio:** en una región sin carga, $V$ en un punto es el promedio de $V$ sobre cualquier esfera centrada en él. Por eso $V$ no tiene máximos ni mínimos locales (no hay equilibrio electrostático estable, teorema de Earnshaw).

**Método de imágenes.** Se reemplaza la frontera por cargas ficticias fuera de la región de interés que reproducen las condiciones de borde:
- Plano conductor a tierra: imagen $-q$ en el punto reflejado.
- Esfera a tierra de radio $R$ con $q$ a distancia $a>R$ del centro:

```math
q'=-\frac{R}{a}\,q,\qquad b=\frac{R^2}{a}\quad(\text{distancia de la imagen al centro}).
```

Las funciones de Green de Dirichlet del plano y de la esfera se obtienen exactamente así.

**Separación de variables.** En cartesianas aparecen exponenciales y senos (series de Fourier). En esféricas con simetría azimutal:

```math
\boxed{V(r,\theta)=\sum_{\ell=0}^{\infty}\left(A_\ell\,r^\ell+\frac{B_\ell}{r^{\ell+1}}\right)P_\ell(\cos\theta)},\qquad
\int_{-1}^{1}P_\ell(x)P_{\ell'}(x)\,dx=\frac{2}{2\ell+1}\delta_{\ell\ell'} .
```

Los polinomios de Legendre tienen la **función generatriz** $\dfrac{1}{\sqrt{1-2xt+t^2}}=\sum_\ell t^\ell P_\ell(x)$, que es exactamente el origen de la expansión multipolar. Ejemplo canónico: esfera conductora en un campo uniforme, $V=-E_0\left(r-\dfrac{R^3}{r^2}\right)\cos\theta$, $\sigma=3\varepsilon_0E_0\cos\theta$.

*Material:* `Electro 2026/Funciones_de_green.pdf`, `Guia_Quiz5_Green_Radiacion.pdf`.

---

<a id="s5"></a>
## 5. Expansión multipolar

Lejos de una distribución localizada, la función generatriz de la sección 4 da

```math
\frac{1}{\mathfrak R}=\sum_{\ell=0}^\infty\frac{r'^\ell}{r^{\ell+1}}P_\ell(\cos\alpha)
\;\Longrightarrow\;
V(\mathbf r)=\frac{1}{4\pi\varepsilon_0}\left[\frac{Q}{r}+\frac{\mathbf p\cdot\hat{\mathbf r}}{r^2}+\frac12\sum_{ij}Q_{ij}\frac{\hat r_i\hat r_j}{r^3}+\cdots\right]
```

con

```math
Q=\int\rho\,d\tau,\qquad \mathbf p=\int\mathbf r'\rho\,d\tau',\qquad Q_{ij}=\int\rho\,(3x_i'x_j'-r'^2\delta_{ij})\,d\tau' .
```

Cada término decae una potencia más rápido. El primer momento no nulo es independiente del origen; los siguientes no (por ejemplo, $\mathbf p$ depende del origen si $Q\neq0$).

**Dipolo eléctrico.**

```math
\boxed{\mathbf E_{\rm dip}=\frac{1}{4\pi\varepsilon_0}\frac{3(\mathbf p\cdot\hat{\mathbf r})\hat{\mathbf r}-\mathbf p}{r^3}}
\qquad\text{en campo externo:}\quad \mathbf N=\mathbf p\times\mathbf E,\quad U=-\mathbf p\cdot\mathbf E,\quad \mathbf F=\nabla(\mathbf p\cdot\mathbf E).
```

La expansión es la base de la radiación multipolar (sección 12): el término dipolar domina cuando la fuente es pequeña frente a $\lambda$.

---

<a id="s6"></a>
## 6. Magnetostática

**Corrientes y conservación de la carga.**

```math
\boxed{\frac{\partial\rho}{\partial t}+\nabla\cdot\mathbf J=0}\qquad(\text{magnetostática: }\nabla\cdot\mathbf J=0),\qquad \mathbf F=q(\mathbf E+\mathbf v\times\mathbf B).
```

La fuerza magnética no hace trabajo: $\mathbf v\cdot(\mathbf v\times\mathbf B)=0$.

**Biot–Savart** (la "ley del inverso del cuadrado" magnética) y las ecuaciones que implica:

```math
\mathbf B(\mathbf r)=\frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r')\times\hat{\boldsymbol{\mathfrak R}}}{\mathfrak R^2}\,d\tau'
\quad\Longrightarrow\quad
\boxed{\nabla\cdot\mathbf B=0,\qquad\nabla\times\mathbf B=\mu_0\mathbf J},\qquad \oint\mathbf B\cdot d\mathbf l=\mu_0I_{\rm enc}\ \ \text{(Ampère)} .
```

**Potencial vector.** $\nabla\cdot\mathbf B=0\Rightarrow\mathbf B=\nabla\times\mathbf A$. $\mathbf A$ está definido salvo un gradiente; en el **gauge de Coulomb** $\nabla\cdot\mathbf A=0$:

```math
\nabla^2\mathbf A=-\mu_0\mathbf J,\qquad \mathbf A(\mathbf r)=\frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r')}{\mathfrak R}\,d\tau' .
```

Cada componente cumple una ecuación de Poisson: todo lo aprendido en electrostática se reutiliza.

**Condiciones de borde** con corriente superficial $\mathbf K$:

```math
B^\perp\ \text{continuo},\qquad \mathbf B_{\rm arriba}-\mathbf B_{\rm abajo}=\mu_0\,\mathbf K\times\hat{\mathbf n},\qquad \mathbf A\ \text{continuo}.
```

**Expansión multipolar magnética.** No hay monopolo. El término dominante es el **dipolo magnético**:

```math
\mathbf m=I\int d\mathbf a=\frac12\int\mathbf r\times\mathbf J\,d\tau,\qquad
\mathbf A_{\rm dip}=\frac{\mu_0}{4\pi}\frac{\mathbf m\times\hat{\mathbf r}}{r^2},\qquad
\boxed{\mathbf B_{\rm dip}=\frac{\mu_0}{4\pi}\frac{3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r}-\mathbf m}{r^3}}
```

En un campo externo: $\mathbf N=\mathbf m\times\mathbf B$, $U=-\mathbf m\cdot\mathbf B$ y $\mathbf F=\nabla(\mathbf m\cdot\mathbf B)$. Lejos de la fuente, el dipolo magnético tiene la misma forma que el eléctrico.

---

<a id="s7"></a>
## 7. Campos en la materia: dieléctricos y materiales magnéticos

**Dieléctricos.** Un medio polarizado (momento dipolar por unidad de volumen $\mathbf P$) equivale a cargas ligadas:

```math
\rho_b=-\nabla\cdot\mathbf P,\qquad \sigma_b=\mathbf P\cdot\hat{\mathbf n};\qquad
\mathbf D\equiv\varepsilon_0\mathbf E+\mathbf P,\qquad \boxed{\nabla\cdot\mathbf D=\rho_f} .
```

En un medio lineal, $\mathbf P=\varepsilon_0\chi_e\mathbf E$ y $\mathbf D=\varepsilon\mathbf E$, con $\varepsilon=\varepsilon_0(1+\chi_e)=\varepsilon_0\varepsilon_r$. Condiciones de borde: $D^\perp_{\rm arriba}-D^\perp_{\rm abajo}=\sigma_f$ y $\mathbf E^\parallel$ continuo.

Resultados canónicos:

```math
\text{esfera uniformemente polarizada:}\ \ \mathbf E_{\rm in}=-\frac{\mathbf P}{3\varepsilon_0};\qquad
\text{esfera dieléctrica en }\mathbf E_0\text{ uniforme:}\ \ \mathbf E_{\rm in}=\frac{3}{\varepsilon_r+2}\,\mathbf E_0 .
```

**Materiales magnéticos.** Una magnetización $\mathbf M$ equivale a corrientes ligadas:

```math
\mathbf J_b=\nabla\times\mathbf M,\qquad \mathbf K_b=\mathbf M\times\hat{\mathbf n};\qquad
\mathbf H\equiv\frac{\mathbf B}{\mu_0}-\mathbf M,\qquad \boxed{\nabla\times\mathbf H=\mathbf J_f} .
```

- **Diamagnetismo:** $\chi_m<0$, inducido (Lenz a escala atómica).
- **Paramagnetismo:** $\chi_m>0$, alineación de momentos permanentes contra la agitación térmica.
- **Ferromagnetismo:** alineación espontánea en dominios e histéresis, por lo que no es lineal.

En medios lineales, $\mathbf M=\chi_m\mathbf H$ y $\mathbf B=\mu\mathbf H$. Condiciones de borde: $B^\perp$ continuo, $\mathbf H^\parallel_{\rm arriba}-\mathbf H^\parallel_{\rm abajo}=\mathbf K_f\times\hat{\mathbf n}$.

**Potencial escalar magnético.** Donde $\mathbf J_f=0$, $\mathbf H=-\nabla\Phi_M$ y de $\nabla\cdot\mathbf B=0$:

```math
\nabla^2\Phi_M=\nabla\cdot\mathbf M\qquad\Big(\rho_M\equiv-\nabla\cdot\mathbf M\ \text{actúa como "carga magnética"}\Big);\qquad
\text{esfera uniformemente magnetizada:}\ \ \mathbf B_{\rm in}=\frac23\mu_0\mathbf M .
```

---

<a id="s8"></a>
## 8. Inducción y ecuaciones de Maxwell

**Ley de Faraday.** Para circuitos fijos o en movimiento, la fem es la integral de la fuerza por unidad de carga, y se cumple la *regla del flujo*:

```math
\mathcal E=\oint(\mathbf E+\mathbf v\times\mathbf B)\cdot d\mathbf l=-\frac{d\Phi_B}{dt}
\qquad\Longrightarrow\qquad
\boxed{\nabla\times\mathbf E=-\frac{\partial\mathbf B}{\partial t}} .
```

En un circuito en movimiento la fem viene de $\mathbf v\times\mathbf B$ (fuerza magnética); en uno fijo con $\mathbf B$ variable, del campo eléctrico inducido. Que ambas den $-d\Phi/dt$ fue una de las motivaciones de Einstein para la relatividad.

**Corriente de desplazamiento.** Ampère con $\nabla\times\mathbf B=\mu_0\mathbf J$ es inconsistente con la conservación de la carga, porque $\nabla\cdot(\nabla\times\mathbf B)=0$ exigiría $\nabla\cdot\mathbf J=0$. Maxwell lo corrige agregando $\mu_0\varepsilon_0\,\partial\mathbf E/\partial t$:

```math
\boxed{\begin{aligned}
\nabla\cdot\mathbf E&=\frac{\rho}{\varepsilon_0}, &\qquad \nabla\times\mathbf E&=-\frac{\partial\mathbf B}{\partial t},\\
\nabla\cdot\mathbf B&=0, &\qquad \nabla\times\mathbf B&=\mu_0\mathbf J+\mu_0\varepsilon_0\frac{\partial\mathbf E}{\partial t}.
\end{aligned}}
```

En la materia: $\nabla\cdot\mathbf D=\rho_f$, $\nabla\times\mathbf H=\mathbf J_f+\partial\mathbf D/\partial t$, y las ecuaciones homogéneas no cambian. **Consecuencia:** en el vacío, $\nabla^2\mathbf E=\mu_0\varepsilon_0\,\partial_t^2\mathbf E$. Son ondas que viajan a $c=1/\sqrt{\mu_0\varepsilon_0}$: la luz es electromagnética.

---

<a id="s9"></a>
## 9. Leyes de conservación: energía, momento y momento angular

**Teorema de Poynting.** Es la conservación local de la energía:

```math
\boxed{\frac{\partial u}{\partial t}+\nabla\cdot\mathbf S=-\mathbf J\cdot\mathbf E},\qquad
u=\frac{\varepsilon_0}{2}E^2+\frac{1}{2\mu_0}B^2,\qquad \mathbf S=\frac{1}{\mu_0}\mathbf E\times\mathbf B .
```

$-\mathbf J\cdot\mathbf E$ es la potencia que el campo entrega a las cargas por unidad de volumen, y $\mathbf S$ es el flujo de energía. Para una onda plana, $\mathbf S=cu\,\hat{\mathbf k}$ y $\langle S\rangle=\frac12\varepsilon_0cE_0^2$.

**Momento del campo y tensor de esfuerzos de Maxwell.**

```math
T_{ij}=\varepsilon_0\Big(E_iE_j-\tfrac12\delta_{ij}E^2\Big)+\frac{1}{\mu_0}\Big(B_iB_j-\tfrac12\delta_{ij}B^2\Big),\qquad
\mathbf g=\mu_0\varepsilon_0\,\mathbf S=\frac{\mathbf S}{c^2},
```

```math
\boxed{\mathbf F_{\rm sobre\ cargas\ en\ \mathcal V}=\oint_{\mathcal S}\overleftrightarrow{\mathbf T}\cdot d\mathbf a-\mu_0\varepsilon_0\frac{d}{dt}\int_{\mathcal V}\mathbf S\,d\tau}
\qquad\Longleftrightarrow\qquad
\frac{\partial}{\partial t}\big(\mathbf g_{\rm mec}+\mathbf g\big)=\nabla\cdot\overleftrightarrow{\mathbf T} .
```

$-\overleftrightarrow{\mathbf T}$ es el flujo de momento. Ejemplos: presión de radiación $P=I/c$ (absorción) o $2I/c$ (reflexión), fuerza sobre conductores y solenoides.

**Momento angular del campo.** $\boldsymbol\ell=\mathbf r\times\mathbf g$. Ejemplo clásico: al apagar la corriente de un solenoide con un cilindro cargado alrededor, el cilindro empieza a girar. El momento angular "guardado" en los campos estáticos se transfiere a las cargas, de modo que el momento angular total se conserva.

---

<a id="s10"></a>
## 10. Potenciales, gauge y potenciales retardados

Las ecuaciones homogéneas se resuelven idénticamente con potenciales:

```math
\mathbf B=\nabla\times\mathbf A,\qquad \mathbf E=-\nabla V-\frac{\partial\mathbf A}{\partial t};\qquad
\boxed{\mathbf A\to\mathbf A+\nabla\lambda,\qquad V\to V-\frac{\partial\lambda}{\partial t}}\ \ \text{(invariancia de gauge)} .
```

- **Gauge de Coulomb** ($\nabla\cdot\mathbf A=0$): $\nabla^2V=-\rho/\varepsilon_0$ es instantáneo (no viola causalidad: $\mathbf E$ y $\mathbf B$ sí son retardados), pero la ecuación para $\mathbf A$ queda acoplada.
- **Gauge de Lorenz** ($\nabla\cdot\mathbf A+\frac1{c^2}\partial_tV=0$) desacopla todo en ecuaciones de onda:

```math
\boxed{\Box V=\frac{\rho}{\varepsilon_0},\qquad \Box\mathbf A=\mu_0\mathbf J},\qquad \Box\equiv\frac{1}{c^2}\frac{\partial^2}{\partial t^2}-\nabla^2 .
```

**Función de Green retardada y potenciales retardados.**

```math
G_{\rm ret}(\mathbf r,t;\mathbf r',t')=\frac{\delta\!\left(t'-t+\mathfrak R/c\right)}{4\pi\mathfrak R}
\quad\Longrightarrow\quad
\boxed{V(\mathbf r,t)=\frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf r',t_r)}{\mathfrak R}\,d\tau',\qquad
\mathbf A(\mathbf r,t)=\frac{\mu_0}{4\pi}\int\frac{\mathbf J(\mathbf r',t_r)}{\mathfrak R}\,d\tau'},\qquad t_r=t-\frac{\mathfrak R}{c}.
```

La información viaja por el **cono de luz pasado**: el campo en $(\mathbf r,t)$ depende de la fuente en el tiempo retardado.

**Simulaciones:** dos gauges con el mismo $\mathbf B$; el vínculo de Gauss solo se mantiene si $\partial_\mu J^\mu=0$; cono de luz y gauge de Lorenz verificados numéricamente → [Guía 6, cuaderno 03](EI_Simulaciones/Guia_06_Radiacion_Covarianza/03_maxwell_accion_gauge.ipynb).

---

<a id="s11"></a>
## 11. Liénard–Wiechert: campos de una carga en movimiento

Para una carga puntual en $\mathbf w(t)$, con $\boldsymbol{\mathfrak R}=\mathbf r-\mathbf w(t_r)$, $\hat{\mathbf n}=\hat{\boldsymbol{\mathfrak R}}$, $\boldsymbol\beta=\mathbf v/c$ y $\kappa=1-\hat{\mathbf n}\cdot\boldsymbol\beta$, la integral retardada da un jacobiano $1/\kappa$, porque $dt/dt_r=\kappa$:

```math
\boxed{V=\frac{q}{4\pi\varepsilon_0}\frac{1}{\kappa\,\mathfrak R}\bigg|_{t_r},\qquad \mathbf A=\frac{\mathbf v}{c^2}V}
```

Derivando se obtienen los campos, que se separan en dos términos:

```math
\mathbf E=\underbrace{\frac{q}{4\pi\varepsilon_0}\frac{(\hat{\mathbf n}-\boldsymbol\beta)(1-\beta^2)}{\kappa^3\mathfrak R^2}}_{\text{velocidad}\ \propto\ \mathfrak R^{-2}}
+\underbrace{\frac{q}{4\pi\varepsilon_0c}\frac{\hat{\mathbf n}\times\big[(\hat{\mathbf n}-\boldsymbol\beta)\times\dot{\boldsymbol\beta}\big]}{\kappa^3\mathfrak R}}_{\text{radiación (aceleración)}\ \propto\ \mathfrak R^{-1}},
\qquad \mathbf B=\frac1c\,\hat{\mathbf n}\times\mathbf E .
```

- Solo el término de **aceleración** decae como $1/\mathfrak R$ y transporta energía al infinito, porque $\mathfrak R^2S$ es finito: **solo las cargas aceleradas radian**.
- **Carga con velocidad constante:** el campo apunta desde la posición *presente* y se comprime en un "panqueque" perpendicular al movimiento (Heaviside):

```math
\mathbf E=\frac{q}{4\pi\varepsilon_0}\frac{1-\beta^2}{\big(1-\beta^2\sin^2\theta\big)^{3/2}}\frac{\hat{\mathbf R}}{R^2},\qquad \mathbf B=\frac{1}{c^2}\mathbf v\times\mathbf E .
```

**Simulaciones:** compresión temporal $dt/dt_r=\kappa$; escalamiento $R^{-2}$ vs. $R^{-1}$; flujo por esferas; **construcción de Purcell** → [Guía 5, cuaderno 01](EI_Simulaciones/Guia_05_Radiacion_Lienard_Wiechert/01_tiempo_retardado_lienard_wiechert.ipynb). Campo de Heaviside obtenido boosteando Coulomb → [Guía 6, cuaderno 02](EI_Simulaciones/Guia_06_Radiacion_Covarianza/02_cuadrivectores_tensor.ipynb).

---

<a id="s12"></a>
## 12. Radiación

**Carga lenta: fórmula de Larmor.** Para $v\ll c$, el campo de radiación es $\mathbf E_{\rm rad}=\frac{q}{4\pi\varepsilon_0c^2\mathfrak R}\,\hat{\mathbf n}\times(\hat{\mathbf n}\times\mathbf a)$. Este campo es **transversal** y proporcional a la proyección de $\mathbf a$ perpendicular a la línea de visión:

```math
\frac{dP}{d\Omega}=\frac{q^2a^2}{16\pi^2\varepsilon_0c^3}\sin^2\theta,\qquad
\boxed{P=\frac{q^2a^2}{6\pi\varepsilon_0c^3}=\frac{\mu_0q^2a^2}{6\pi c}}\quad\text{(Larmor)} .
```

El patrón es una "dona" alrededor de $\mathbf a$: nulo a lo largo de la aceleración y máximo perpendicular a ella.

**Radiación dipolar eléctrica** (fuente de tamaño $d\ll\lambda$; se reemplaza $q\mathbf a\to\ddot{\mathbf p}$). Con $\mathbf p=p_0\cos\omega t\,\hat{\mathbf z}$:

```math
\mathbf E=-\frac{\mu_0p_0\omega^2}{4\pi}\frac{\sin\theta}{r}\cos[\omega(t-r/c)]\,\hat{\boldsymbol\theta},\qquad
\boxed{\langle P\rangle=\frac{\mu_0p_0^2\omega^4}{12\pi c}}\ \ (\propto\omega^4:\ \text{cielo azul});\qquad
P(t)=\frac{\mu_0\,|\ddot{\mathbf p}(t_r)|^2}{6\pi c}\ \ \text{(fuente arbitraria)}.
```

**Radiación dipolar magnética.** Con $\mathbf m=m_0\cos\omega t\,\hat{\mathbf z}$, $\langle P\rangle=\dfrac{\mu_0m_0^2\omega^4}{12\pi c^3}$, suprimida respecto de la eléctrica por un factor $(m_0/p_0c)^2\sim(v/c)^2$. Los campos tienen la forma de los del dipolo eléctrico con $\mathbf E\leftrightarrow c\mathbf B$.

**Carga relativista: Liénard y efecto faro.**

```math
\boxed{P=\frac{q^2\gamma^6}{6\pi\varepsilon_0c^3}\Big(a^2-\big|\boldsymbol\beta\times\mathbf a\big|^2\Big)},\qquad
\frac{dP(t_r)}{d\Omega}=\frac{q^2}{16\pi^2\varepsilon_0c^3}\frac{\big|\hat{\mathbf n}\times[(\hat{\mathbf n}-\boldsymbol\beta)\times\mathbf a]\big|^2}{\kappa^5}.
```

Para $\mathbf a\parallel\mathbf v$, $P\propto\gamma^6a^2$; para $\mathbf a\perp\mathbf v$ (sincrotrón), $P\propto\gamma^4a^2$. El factor $\kappa^{-5}$ concentra la emisión en un cono de ángulo $\sim1/\gamma$ hacia adelante. El detector mide $\kappa^{-6}$, porque cuenta en su propio tiempo $dt=\kappa\,dt_r$.

**Reacción de radiación (Abraham–Lorentz).** Una carga que radia pierde energía, así que debe sentir una fuerza:

```math
\mathbf F_{\rm rad}=\frac{\mu_0q^2}{6\pi c}\,\dot{\mathbf a}=m\tau\,\dot{\mathbf a},\qquad \tau=\frac{q^2}{6\pi\varepsilon_0mc^3}\ \ (\tau_e\approx6.3\times10^{-24}\,\mathrm s).
```

Esta ecuación tiene soluciones desbocadas y preaceleración: la electrodinámica clásica de partículas puntuales no es del todo consistente, por el mismo problema de autoenergía de la sección 3. Para un oscilador, el tiempo de decaimiento es $\tau_E=1/(\omega^2\tau)\sim10\,$ns en el visible, el orden de las vidas medias atómicas.

**Simulaciones:**
- Órbita elíptica, dipolo de Hertz (líneas que se desprenden), efecto faro → [Guía 6, cuaderno 01](EI_Simulaciones/Guia_06_Radiacion_Covarianza/01_radiacion.ipynb).
- Larmor con fotones Monte Carlo, dipolo rotante, amortiguamiento radiativo → [Guía 5, cuaderno 02](EI_Simulaciones/Guia_05_Radiacion_Lienard_Wiechert/02_larmor_dipolos.ipynb).
- Colapso del átomo clásico, enfriamiento sincrotrón, $\kappa^{-5}$ vs. $\kappa^{-6}$ → [Guía 5, cuaderno 03](EI_Simulaciones/Guia_05_Radiacion_Lienard_Wiechert/03_radiacion_relativista.ipynb).

---

<a id="s13"></a>
## 13. Relatividad especial y formulación covariante

**Transformaciones de Lorentz y cuadrivectores.** Con $x^\mu=(ct,\mathbf x)$ y $\eta_{\mu\nu}=\mathrm{diag}(1,-1,-1,-1)$, un boost en $x$ es

```math
ct'=\gamma(ct-\beta x),\quad x'=\gamma(x-\beta ct),\qquad s^2=x^\mu x_\mu=c^2t^2-|\mathbf x|^2\ \text{invariante};\qquad
U^\mu=\gamma(c,\mathbf v),\quad p^\mu=\Big(\frac{\mathcal E}{c},\mathbf p\Big),\quad p^\mu p_\mu=m^2c^2 .
```

**Objetos electromagnéticos.**

```math
J^\mu=(c\rho,\mathbf J),\qquad A^\mu=\Big(\frac{V}{c},\mathbf A\Big),\qquad \boxed{F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu},\qquad
F^{\mu\nu}=\begin{pmatrix}0&-E_x/c&-E_y/c&-E_z/c\\ E_x/c&0&-B_z&B_y\\ E_y/c&B_z&0&-B_x\\ E_z/c&-B_y&B_x&0\end{pmatrix}.
```

$\mathbf E$ y $\mathbf B$ no son vectores independientes: son las 6 componentes de un tensor antisimétrico. Bajo un boost con velocidad $\mathbf v$ se mezclan:

```math
\mathbf E'_\parallel=\mathbf E_\parallel,\quad \mathbf B'_\parallel=\mathbf B_\parallel,\qquad
\mathbf E'_\perp=\gamma(\mathbf E+\mathbf v\times\mathbf B)_\perp,\quad \mathbf B'_\perp=\gamma\Big(\mathbf B-\frac{\mathbf v\times\mathbf E}{c^2}\Big)_\perp .
```

**Invariantes de Lorentz** (lo único que todos los observadores acuerdan):

```math
\boxed{F_{\mu\nu}F^{\mu\nu}=2\Big(B^2-\frac{E^2}{c^2}\Big),\qquad \mathbf E\cdot\mathbf B}
```

Si $\mathbf E\cdot\mathbf B=0$: con $E>cB$ existe un sistema con $\mathbf B'=0$; con $E<cB$, uno con $\mathbf E'=0$ (la deriva $\mathbf E\times\mathbf B$). Si ambos invariantes son nulos, el campo es **nulo**, como la radiación.

**Maxwell covariante, Bianchi y fuerza de Lorentz.**

```math
\boxed{\partial_\mu F^{\mu\nu}=\mu_0J^\nu}\ \ (\text{Gauss + Ampère–Maxwell}),\qquad
\boxed{\partial_\lambda F_{\mu\nu}+\partial_\mu F_{\nu\lambda}+\partial_\nu F_{\lambda\mu}=0}\ \ (\nabla\cdot\mathbf B=0\ \text{+ Faraday}),
```

```math
\frac{dp^\mu}{d\tau}=qF^{\mu\nu}U_\nu\quad\Longrightarrow\quad \frac{d\mathbf p}{dt}=q(\mathbf E+\mathbf v\times\mathbf B),\qquad \frac{d\mathcal E}{dt}=q\,\mathbf E\cdot\mathbf v .
```

La conservación de la carga, $\partial_\mu J^\mu=0$, es *obligatoria*: se sigue de $\partial_\nu\partial_\mu F^{\mu\nu}\equiv0$.

**Principio de acción y gauge.**

```math
S=\int d^4x\left[-\frac{1}{4\mu_0}F_{\mu\nu}F^{\mu\nu}-J^\mu A_\mu\right]\ \xrightarrow{\ \delta A_\nu\ }\ \partial_\mu F^{\mu\nu}=\mu_0J^\nu;\qquad
A_\mu\to A_\mu-\partial_\mu\chi\ \Rightarrow\ \delta S=\int\chi\,\partial_\mu J^\mu\,d^4x .
```

La acción es invariante de gauge si y solo si la carga se conserva. En gauge de Lorenz ($\partial_\mu A^\mu=0$) queda $\Box A^\nu=\mu_0J^\nu$, cuya solución retardada para una carga puntual son los potenciales de Liénard–Wiechert, $A^\mu\propto U^\mu/(U\cdot\mathfrak R)$.

**Tensor de energía-momento.** $T^{\mu\nu}=\frac{1}{\mu_0}\Big[F^{\mu\alpha}F_\alpha{}^{\nu}+\tfrac14\eta^{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\Big]$, con $T^{00}=u$ y $T^{0i}=S_i/c$. Agrupa el teorema de Poynting y el tensor de Maxwell de la sección 9 en $\partial_\mu T^{\mu\nu}=-F^{\nu\lambda}J_\lambda$. Para una onda plana, $T^{\mu\nu}\propto k^\mu k^\nu$ con $k^\mu k_\mu=0$: un flujo de "partículas sin masa".

**Simulaciones:** cuadricorriente en el espacio-tiempo, boosts de campos, integración de la ecuación covariante (la deriva es un boost), acción discreta en una red, plano de invariantes → [Guía 6](EI_Simulaciones/Guia_06_Radiacion_Covarianza/).

---

<p align="center"><a href="README.md">← volver al README</a> · <a href="#índice">↑ índice</a></p>

*Bibliografía:* D. J. Griffiths, *Introduction to Electrodynamics* (4.ª ed.), caps. 1–12 · J. D. Jackson, *Classical Electrodynamics* (3.ª ed.), caps. 1–6, 9, 11, 14 · L. D. Landau y E. M. Lifshitz, *The Classical Theory of Fields*.