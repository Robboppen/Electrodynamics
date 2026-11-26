"""
em_utils.py — Herramientas comunes para las simulaciones de la Guía 6
(Electromagnetismo Intermedio, U. de Valparaíso, 2026).

Convenciones (idénticas a la guía):
    x^mu = (ct, x, y, z),  eta = diag(1, -1, -1, -1)
    A^mu = (phi/c, A),     J^mu = (c rho, J)
    F_{mu nu} = d_mu A_nu - d_nu A_mu
    =>  F_{0i} = E_i / c,   F_{ij} = -eps_{ijk} B_k
        F^{0i} = -E_i / c,  F^{ij} = -eps_{ijk} B_k

Unidades: salvo que se diga lo contrario, trabajamos en unidades
adimensionales con c = 1 y q/(4 pi eps0) = 1. Las fórmulas se escriben
con c explícito para que el paso a SI sea directo.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# ---------------------------------------------------------------------------
# Rutas y estilo
# ---------------------------------------------------------------------------
RAIZ = os.path.dirname(os.path.abspath(__file__))
DIR_FIG = os.path.join(RAIZ, "figuras")
os.makedirs(DIR_FIG, exist_ok=True)

ETA = np.diag([1.0, -1.0, -1.0, -1.0])

# Símbolo de Levi-Civita en 3D como arreglo (i, j, k)
LEVI = np.zeros((3, 3, 3))
for (i, j, k), s in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
                     (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
    LEVI[i, j, k] = s

COLORES = {"E": "#c0392b", "B": "#2471a3", "S": "#7d3c98",
           "a": "#d35400", "v": "#1e8449", "carga": "#17202a"}


def estilo():
    """Estilo de figuras sobrio y legible (sirve igual en GitHub y en papel)."""
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "animation.embed_limit": 50,
    })


def guardar(fig, nombre: str):
    """Guarda una figura en figuras/<nombre>.png y la devuelve."""
    ruta = os.path.join(DIR_FIG, f"{nombre}.png")
    fig.savefig(ruta, bbox_inches="tight")
    print(f"[guardado] figuras/{nombre}.png")
    return ruta


def guardar_gif(anim: animation.FuncAnimation, nombre: str, fps: int = 20, dpi: int = 80):
    """Guarda una animación como GIF en figuras/<nombre>.gif."""
    ruta = os.path.join(DIR_FIG, f"{nombre}.gif")
    anim.save(ruta, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    print(f"[guardado] figuras/{nombre}.gif  ({os.path.getsize(ruta)/1e6:.2f} MB)")
    return ruta


# ---------------------------------------------------------------------------
# Tensor electromagnético
# ---------------------------------------------------------------------------
def F_abajo(E, B, c: float = 1.0):
    """F_{mu nu} a partir de E y B (vectores de 3 componentes, admite arreglos (...,3))."""
    E = np.asarray(E, float)
    B = np.asarray(B, float)
    F = np.zeros(E.shape[:-1] + (4, 4))
    F[..., 0, 1:] = E / c
    F[..., 1:, 0] = -E / c
    F[..., 1:, 1:] = -np.einsum("ijk,...k->...ij", LEVI, B)
    return F


def subir(F_low):
    """F^{mu nu} = eta^{mu a} eta^{nu b} F_{ab}."""
    return np.einsum("ma,...ab,bn->...mn", ETA, F_low, ETA)


def campos_desde_F_arriba(F_up, c: float = 1.0):
    """(E, B) a partir de F^{mu nu}:  E_i = c F^{i0},  B_k = -1/2 eps_{ijk} F^{ij}."""
    E = c * F_up[..., 1:, 0]
    B = -0.5 * np.einsum("ijk,...ij->...k", LEVI, F_up[..., 1:, 1:])
    return E, B


def invariantes(E, B, c: float = 1.0):
    """Devuelve (E^2 - c^2 B^2, E·B)."""
    E = np.asarray(E, float)
    B = np.asarray(B, float)
    return (np.sum(E * E, -1) - c**2 * np.sum(B * B, -1), np.sum(E * B, -1))


def FF(F_low):
    """Contracción completa F_{mu nu} F^{mu nu}."""
    return np.einsum("...mn,...mn->...", F_low, subir(F_low))


# ---------------------------------------------------------------------------
# Transformaciones de Lorentz
# ---------------------------------------------------------------------------
def matriz_boost(beta_vec):
    """Lambda^mu_nu para un boost puro hacia un sistema S' con velocidad beta_vec*c
    respecto de S. Actúa sobre cuadrivectores contravariantes: x'^mu = Lambda^mu_nu x^nu."""
    b = np.asarray(beta_vec, float)
    b2 = b @ b
    if b2 >= 1:
        raise ValueError("|beta| debe ser < 1")
    L = np.eye(4)
    if b2 == 0:
        return L
    g = 1.0 / np.sqrt(1 - b2)
    L[0, 0] = g
    L[0, 1:] = -g * b
    L[1:, 0] = -g * b
    L[1:, 1:] += (g - 1) * np.outer(b, b) / b2
    return L


def boost_campos(E, B, beta_vec, c: float = 1.0):
    """Transforma (E, B) a S' usando F'^{mu nu} = Lambda^mu_a Lambda^nu_b F^{ab}."""
    L = matriz_boost(beta_vec)
    Fu = subir(F_abajo(E, B, c))
    Fu_p = np.einsum("ma,...ab,nb->...mn", L, Fu, L)
    return campos_desde_F_arriba(Fu_p, c)


# ---------------------------------------------------------------------------
# Radiación
# ---------------------------------------------------------------------------
def dPdOmega_lienard(n, beta, a, c: float = 1.0):
    """Distribución angular de la potencia emitida por unidad de tiempo RETARDADO,
        dP/dOmega = q^2/(16 pi^2 eps0 c^3) |n x ((n - beta) x a)|^2 / (1 - n·beta)^5,
    devuelta en unidades de q^2/(16 pi^2 eps0 c^3).  n: (...,3), beta y a: (3,)."""
    n = np.asarray(n, float)
    beta = np.asarray(beta, float)
    a = np.asarray(a, float)
    u = np.cross(n - beta, a)
    w = np.cross(n, u)
    return np.sum(w * w, -1) / (1 - n @ beta) ** 5


def esfera(nt: int = 181, nphi: int = 361):
    """Malla (theta, phi) y vectores unitarios n de forma (nt, nphi, 3)."""
    th = np.linspace(0, np.pi, nt)
    ph = np.linspace(0, 2 * np.pi, nphi)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    n = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    return TH, PH, n


def integrar_esfera(f, nmu: int = 800, nphi: int = 400):
    """Integra f(n) sobre la esfera con Gauss-Legendre en mu = cos(theta)
    y la regla del trapecio (exacta para funciones periódicas) en phi."""
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    MU, PH = np.meshgrid(mu, ph, indexing="ij")
    s = np.sqrt(1 - MU**2)
    n = np.stack([s * np.cos(PH), s * np.sin(PH), MU], -1)
    vals = f(n)
    return np.sum(vals * wmu[:, None]) * (2 * np.pi / nphi)


# ---------------------------------------------------------------------------
# Liénard–Wiechert
# ---------------------------------------------------------------------------
def tiempo_retardado(r, t, w_fun, c: float = 1.0, t_min=None, iters: int = 60):
    """Resuelve  t - t_r = |r - w(t_r)|/c  para cada punto r (...,3) por bisección.
    Para |v| < c la función g(t_r) = t - t_r - |r - w(t_r)|/c es estrictamente
    decreciente, así que la raíz es única (un solo punto en el cono de luz pasado)."""
    r = np.asarray(r, float)
    dist_max = np.linalg.norm(r, axis=-1) + 10.0
    lo = (t - dist_max / c) if t_min is None else np.full(r.shape[:-1], t_min)
    lo = np.broadcast_to(lo, r.shape[:-1]).astype(float).copy()
    hi = np.full(r.shape[:-1], float(t))
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        g = t - mid - np.linalg.norm(r - w_fun(mid), axis=-1) / c
        lo = np.where(g > 0, mid, lo)
        hi = np.where(g > 0, hi, mid)
    return 0.5 * (lo + hi)


def campos_LW(r, t, w_fun, v_fun, a_fun, c: float = 1.0):
    """Campos de Liénard–Wiechert de una carga puntual (unidades q/(4 pi eps0) = 1).
    Devuelve (E_vel, E_rad, B, t_r) con
        E_vel = (n - beta)(1 - beta^2) / [(1 - n·beta)^3 R^2]
        E_rad = n x ((n - beta) x a) / [c^2 (1 - n·beta)^3 R]
        B     = n x (E_vel + E_rad) / c
    """
    tr = tiempo_retardado(r, t, w_fun, c)
    w = w_fun(tr)
    beta = v_fun(tr) / c
    acc = a_fun(tr)
    Rv = r - w
    R = np.linalg.norm(Rv, axis=-1, keepdims=True)
    n = Rv / R
    k = 1 - np.sum(n * beta, -1, keepdims=True)
    b2 = np.sum(beta * beta, -1, keepdims=True)
    E_vel = (n - beta) * (1 - b2) / (k**3 * R**2)
    E_rad = np.cross(n, np.cross(n - beta, acc)) / (c**2 * k**3 * R)
    B = np.cross(n, E_vel + E_rad) / c
    return E_vel, E_rad, B, tr


def potenciales_LW(r, t, w_fun, v_fun, c: float = 1.0):
    """Potenciales de Liénard–Wiechert (q/(4 pi eps0) = 1):
        phi = 1 / [R (1 - n·beta)],   A = beta phi / c."""
    tr = tiempo_retardado(r, t, w_fun, c)
    Rv = r - w_fun(tr)
    R = np.linalg.norm(Rv, axis=-1)
    n = Rv / R[..., None]
    beta = v_fun(tr) / c
    phi = 1.0 / (R * (1 - np.sum(n * beta, -1)))
    A = beta * phi[..., None] / c
    return phi, A
