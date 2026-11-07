"""
Minimal P1 finite element solver on a rectangle, used as an independent
reference for the Monte Carlo (Feynman-Kac) solutions.

It solves   L u = b . grad(u) + div(a grad u) = -g   in D = (x0, x1) x (y0, y1),
            u = phi on dD,                           (a constant, symmetric > 0)

with the weak form (the same one used in ``fem/BVPUSC.edp`` for a = I):

    int_D (a grad u) . grad v  -  int_D (b . grad u) v  =  int_D g v      for all v in H^1_0.

The mesh is the usual structured one: nx * ny squares, each split in two
triangles.  ``g`` is integrated with the 3-point edge-midpoint rule (exact for
quadratics).
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.tri as mtri


class FEMSolution:
    def __init__(self, nodes, triangles, u):
        self.nodes, self.triangles, self.u = nodes, triangles, u
        self._tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], triangles)
        self._interp = mtri.LinearTriInterpolator(self._tri, u)

    def __call__(self, points):
        """Evaluate the P1 function at the given points, shape (M, 2)."""
        points = np.atleast_2d(points)
        return np.asarray(self._interp(points[:, 0], points[:, 1]))


def solve_dirichlet_poisson(a, b, g, phi, lo, hi, nx, ny=None):
    """Return a ``FEMSolution`` of  L u = -g,  u = phi on the boundary."""
    ny = nx if ny is None else ny
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    # --- structured mesh ----------------------------------------------------
    xs = np.linspace(lo[0], hi[0], nx + 1)
    ys = np.linspace(lo[1], hi[1], ny + 1)
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    nodes = np.column_stack([XX.ravel(), YY.ravel()])
    idx = np.arange((nx + 1) * (ny + 1)).reshape(nx + 1, ny + 1)
    n00, n10 = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()
    n01, n11 = idx[:-1, 1:].ravel(), idx[1:, 1:].ravel()
    tris = np.vstack([np.column_stack([n00, n10, n11]),
                      np.column_stack([n00, n11, n01])])

    # --- element geometry ---------------------------------------------------
    P = nodes[tris]                                    # (T, 3, 2)
    x, y = P[:, :, 0], P[:, :, 1]
    # gradient of the barycentric basis: grad(phi_i) = (b_i, c_i) / (2 A)
    bi = np.column_stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]])
    ci = np.column_stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]])
    area = 0.5 * ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
                  - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    grad = np.stack([bi, ci], axis=2) / (2.0 * area)[:, None, None]   # (T, 3, 2)

    # --- element matrices -----------------------------------------------------
    # stiffness: K_ij = A * (a grad_j) . grad_i
    K = area[:, None, None] * np.einsum("tik,kl,tjl->tij", grad, a, grad)
    # convection: C_ij = int (b . grad phi_j) phi_i = A/3 * (b . grad_j)
    Cj = np.einsum("tjk,k->tj", grad, b)
    C = (area / 3.0)[:, None, None] * np.broadcast_to(Cj[:, None, :], (len(tris), 3, 3))
    Aloc = K - C

    # load: int g phi_i with the edge-midpoint rule
    mids = 0.5 * (P[:, [0, 1, 2], :] + P[:, [1, 2, 0], :])      # (T, 3, 2): edges 01,12,20
    gm = g(mids.reshape(-1, 2)).reshape(-1, 3)
    # phi_i at the midpoint of edge k is 1/2 if vertex i belongs to edge k, else 0
    W = np.array([[0.5, 0.0, 0.5],     # vertex 0 lies on edges 01 and 20
                  [0.5, 0.5, 0.0],     # vertex 1 lies on edges 01 and 12
                  [0.0, 0.5, 0.5]])    # vertex 2 lies on edges 12 and 20
    F = (area / 3.0)[:, None] * (gm @ W.T)

    rows = np.repeat(tris, 3, axis=1).ravel()
    cols = np.tile(tris, (1, 3)).ravel()
    Amat = sp.coo_matrix((Aloc.ravel(), (rows, cols)), shape=(len(nodes),) * 2).tocsr()
    rhs = np.bincount(tris.ravel(), weights=F.ravel(), minlength=len(nodes))

    # --- Dirichlet conditions -------------------------------------------------
    on_bd = (np.isclose(nodes[:, 0], lo[0]) | np.isclose(nodes[:, 0], hi[0])
             | np.isclose(nodes[:, 1], lo[1]) | np.isclose(nodes[:, 1], hi[1]))
    u = np.zeros(len(nodes))
    u[on_bd] = phi(nodes[on_bd])
    free = ~on_bd
    rhs_free = rhs[free] - Amat[free][:, on_bd] @ u[on_bd]
    u[free] = spla.spsolve(Amat[free][:, free].tocsc(), rhs_free)
    return FEMSolution(nodes, tris, u)
