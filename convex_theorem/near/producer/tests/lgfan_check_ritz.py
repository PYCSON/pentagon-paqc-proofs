"""float Rayleigh-Ritz value of lambda_1 of a convex polygon with the one-piece trial space (product of the edge functions
times a tensor Legendre polynomial): an upper bound of lambda_1 up to rounding, independent of the fan trial space.
Used by lgfan_check.py / lgpoly_check.py."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import lgcert as LG


def ritz(V, d=8):
    V = np.asarray(V, float); c = (V.max(0) + V.min(0)) / 2; sc = (V.max(0) - V.min(0)).max() / 2; U = (V - c) / sc; nv = len(U)
    x, y, wq = LG._quad(U, d + nv + 16); B = np.ones_like(x); Bx = np.zeros_like(x); By = np.zeros_like(x)
    for k in range(nv):
        p, q = U[k], U[(k + 1) % nv]; nx, ny = -(q - p)[1], (q - p)[0]; l = nx * (x - p[0]) + ny * (y - p[1]); Bx, By = Bx * l + B * nx, By * l + B * ny; B = B * l
    LX, DX = LG._legvals(d, x); LY, DY = LG._legvals(d, y); P0 = np.einsum('ip,jp->ijp', LX, LY).reshape((d + 1) ** 2, -1)
    Px = np.einsum('ip,jp->ijp', DX, LY).reshape((d + 1) ** 2, -1); Py = np.einsum('ip,jp->ijp', LX, DY).reshape((d + 1) ** 2, -1)
    F = B * P0; Fx = Bx * P0 + B * Px; Fy = By * P0 + B * Py; M = (F * wq) @ F.T; K = (Fx * wq) @ Fx.T + (Fy * wq) @ Fy.T
    em, Qm = np.linalg.eigh(M); keep = em > 1e-11 * em.max(); T0 = Qm[:, keep] / np.sqrt(em[keep])
    return np.linalg.eigvalsh(T0.T @ K @ T0)[0] / sc ** 2, wq.sum() * sc ** 2
