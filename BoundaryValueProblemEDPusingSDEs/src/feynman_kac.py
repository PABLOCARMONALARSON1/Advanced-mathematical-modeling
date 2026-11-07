"""
Monte Carlo solver for the Dirichlet-Poisson problem via its stochastic
representation (Dynkin / Feynman-Kac formula).

    L u = -g   in D,        u = phi   on dD,

    L = sum_i b_i(x) d/dx_i + sum_{ij} a_ij(x) d^2/(dx_i dx_j),   a = (1/2) sigma sigma^T.

The solution is  u(x) = E^x[ phi(X_tau) + int_0^tau g(X_s) ds ]  where X is the
Ito diffusion  dX = b(X) dt + sigma(X) dW  and  tau  is its first exit time
from D.

The diffusion is discretised with Euler-Maruyama.  Checking the exit only at
the grid times misses the paths that leave and come back inside between two
steps, which biases the exit time upwards and the estimate by O(sqrt(dt))
(Gobet 2000).  Shifting the boundary inward by  BETA * sigma_n * sqrt(dt)
(Broadie-Glasserman-Kou 1997; Gobet-Menozzi 2010), with sigma_n the standard
deviation of the process in the direction normal to the boundary, restores a
bias of order O(dt).  It is enabled by default (``boundary_correction=True``).
"""

import numpy as np

# -zeta(1/2)/sqrt(2): the "continuity correction" constant of discrete barrier
# monitoring (Broadie, Glasserman, Kou 1997).
BETA = 0.5825971579390106


# --------------------------------------------------------------------------
# Domains
# --------------------------------------------------------------------------
class Box:
    """Axis-aligned box  prod_i (lo_i, hi_i)."""

    def __init__(self, lo, hi):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(hi, dtype=float)

    def inside(self, X, A, dt, beta=0.0):
        """True for the paths that are still in the (shrunk) box.

        A is the matrix a(x): shape (n, n) if constant, (M, n, n) otherwise.
        """
        d = np.diagonal(A, axis1=-2, axis2=-1)          # (n,) or (M, n)
        shift = beta * np.sqrt(2.0 * d * dt)            # sigma_n = sqrt(2 a_ii)
        return np.all((X > self.lo + shift) & (X < self.hi - shift), axis=1)

    def project(self, X):
        """Closest point of the closed box (where phi is evaluated at exit)."""
        return np.clip(X, self.lo, self.hi)


class Ball:
    """Open ball of centre ``center`` and radius ``radius``."""

    def __init__(self, center, radius):
        self.center = np.asarray(center, dtype=float)
        self.radius = float(radius)

    def inside(self, X, A, dt, beta=0.0):
        D = X - self.center
        r = np.linalg.norm(D, axis=1)
        nrm = D / np.maximum(r, 1e-300)[:, None]
        if A.ndim == 2:
            q = np.einsum("mi,ij,mj->m", nrm, A, nrm)
        else:
            q = np.einsum("mi,mij,mj->m", nrm, A, nrm)
        return r < self.radius - beta * np.sqrt(2.0 * q * dt)

    def project(self, X):
        D = X - self.center
        r = np.linalg.norm(D, axis=1)
        return self.center + self.radius * D / np.maximum(r, 1e-300)[:, None]


class _CallableDomain:
    """Wraps a plain ``inside(X) -> bool array`` function (no boundary correction)."""

    def __init__(self, fn):
        self.fn = fn

    def inside(self, X, A, dt, beta=0.0):
        return self.fn(X)

    def project(self, X):
        return X


# --------------------------------------------------------------------------
# Solver
# --------------------------------------------------------------------------
def feynman_kac_solver(a, b, g, phi, x_eval, domain, N=1000, dt=1e-3, T=100.0,
                       boundary_correction=True, seed=None, verbose=True):
    """Estimate u(x_eval) for  L u = -g in D,  u = phi on dD.

    Parameters
    ----------
    a : (n, n) array, or callable X(M, n) -> (M, n, n)
        Diffusion matrix, positive definite (a = sigma sigma^T / 2).
    b : (n,) array, or callable X(M, n) -> (M, n)
        Drift.
    g : callable X(M, n) -> (M,), or None (no source term).
    phi : callable X(M, n) -> (M,)
        Dirichlet data.  Evaluated on the boundary point closest to the exit.
    x_eval : (n,) or (P, n) array
        Point(s) where the solution is estimated.  All points are simulated
        together (N independent paths each), so P points cost P*N paths.
    domain : ``Box``, ``Ball`` or callable X(M, n) -> bool (True = inside).
    N : paths per evaluation point.
    dt : Euler-Maruyama time step.
    T : maximum simulated time (paths that have not exited by T are truncated).
    boundary_correction : shift the boundary inward by BETA*sigma_n*sqrt(dt).
        Only available for ``Box`` and ``Ball``.
    seed : seed of the random number generator.

    Returns
    -------
    u_hat, std_error : floats if x_eval is a single point, arrays (P,) otherwise.
        ``std_error`` is the Monte Carlo standard error only; it does not
        include the (deterministic) discretisation bias.
    """
    x_eval = np.asarray(x_eval, dtype=float)
    single = x_eval.ndim == 1
    P = np.atleast_2d(x_eval)
    n_pts, n = P.shape

    if not hasattr(domain, "inside"):
        domain = _CallableDomain(domain)
    beta = BETA if boundary_correction else 0.0

    const_a = not callable(a)
    if const_a:
        A_const = np.asarray(a, dtype=float)
        sigma_const = np.linalg.cholesky(2.0 * A_const)

    if callable(b):
        b_fun = b
    else:
        b_vec = np.asarray(b, dtype=float)
        b_fun = lambda X: np.broadcast_to(b_vec, X.shape)

    rng = np.random.default_rng(seed)
    sqdt = np.sqrt(dt)
    n_steps = int(np.ceil(T / dt))

    M = n_pts * N
    Xa = np.repeat(P, N, axis=0)            # positions of the paths still alive
    acc = np.zeros(M)                       # int_0^t g(X_s) ds of the alive paths
    alive = np.arange(M)                    # their global indices
    value = np.zeros(M)                     # phi(X_tau) + int_0^tau g, per path

    if verbose:
        print(f"Monte Carlo: {n_pts} point(s) x N={N} paths, dt={dt}, "
              f"boundary_correction={boundary_correction}")

    for _ in range(n_steps):
        if alive.size == 0:
            break
        m = Xa.shape[0]
        A = A_const if const_a else np.asarray(a(Xa), dtype=float)

        if g is not None:
            acc += g(Xa) * dt               # left-point rule for the Poisson term

        Z = rng.standard_normal((m, n))
        if const_a:
            noise = Z @ sigma_const.T
        else:
            sigma = np.linalg.cholesky(2.0 * A)
            noise = np.einsum("mij,mj->mi", sigma, Z)
        Xa = Xa + b_fun(Xa) * dt + sqdt * noise

        ins = domain.inside(Xa, A, dt, beta)
        if not ins.all():
            out = ~ins
            value[alive[out]] = phi(domain.project(Xa[out])) + acc[out]
            Xa, acc, alive = Xa[ins], acc[ins], alive[ins]

    if alive.size:                          # not exited by T: truncate
        value[alive] = phi(Xa) + acc
        if verbose:
            print(f"Warning: {alive.size} paths ({100 * alive.size / M:.2f}%) "
                  f"had not exited by T={T}; tau was truncated to T.")

    value = value.reshape(n_pts, N)
    u_hat = value.mean(axis=1)
    std_error = value.std(axis=1, ddof=1) / np.sqrt(N)
    if single:
        return float(u_hat[0]), float(std_error[0])
    return u_hat, std_error
