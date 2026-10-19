"""Reusable BEM core for the R0 revision (bempp-cl, numba backend).

Conventions (bempp-cl): Green's function G = exp(+ik r)/(4 pi r)  (time factor e^{-i w t}).
Momentum equation:  dp/dn = i*w*rho0*v_n   (n = the direction in which v_n is measured).

Functions
---------
closed_neumann_bm(grid, vn, k)          Burton-Miller exterior Neumann solve -> (P, p_coeffs)
open_screen_power(grid, vn, k)          two-sided pressure-jump screen       -> (P, mu_coeffs)
far_field_power_dl(grid, mu, k)         far-field power of a double-layer density mu (screen)
far_field_power_closed(grid, p, q, k)   far-field power of a closed-surface (p, dp/dn) pair
sinc_power_mesh(grid, vn, k)            mesh-based sinc (Rayleigh-type) model power
All vn arguments are P1 *vertex* values (length grid.number_of_vertices, indexed by grid vertex).
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import bempp_cl.api as bem

RHO0 = 1.21
C0 = 343.0


# ----------------------------------------------------------------------------- helpers
def _vertex_dof_map(space, grid):
    """dof index of every grid vertex (-1 if the vertex has no dof)."""
    l2g = np.asarray(space.local2global)          # (n_el, 3) global dof per local vertex
    mult = np.asarray(space.local_multipliers)
    vmap = -np.ones(grid.number_of_vertices, dtype=int)
    el = grid.elements                            # (3, n_el)
    for loc in range(3):
        good = mult[:, loc] != 0
        vmap[el[loc, good]] = l2g[good, loc]
    return vmap


def full_p1(grid):
    """P1 continuous space with ALL vertex dofs and the vertex->dof permutation."""
    space = bem.function_space(grid, "P", 1, include_boundary_dofs=True)
    vmap = _vertex_dof_map(space, grid)
    assert (vmap >= 0).all() and space.global_dof_count == grid.number_of_vertices
    return space, vmap


def boundary_vertices(grid):
    """Boolean mask of vertices lying on the free (open) boundary of the mesh."""
    el = grid.elements.T
    edges = np.sort(np.vstack([el[:, [0, 1]], el[:, [1, 2]], el[:, [2, 0]]]), axis=1)
    uniq, cnt = np.unique(edges, axis=0, return_counts=True)
    bnd = np.zeros(grid.number_of_vertices, dtype=bool)
    bnd[uniq[cnt == 1].ravel()] = True
    return bnd


def mass_matrix_vertex(grid):
    """P1 mass matrix in VERTEX ordering (sparse, real, symmetric)."""
    V = grid.vertices.T
    E = grid.elements.T
    a = np.linalg.norm(np.cross(V[E[:, 1]] - V[E[:, 0]], V[E[:, 2]] - V[E[:, 0]]), axis=1) / 2
    loc = (np.ones((3, 3)) + np.eye(3)) / 12.0
    I = np.repeat(E, 3, axis=1).ravel()
    J = np.tile(E, (1, 3)).ravel()
    vals = (a[:, None, None] * loc[None]).ravel()
    return sp.csr_matrix((vals, (I, J)), shape=(len(V), len(V)))


def _op_vertex_order(op, vmap):
    """dense weak form of a P1-P1 operator, reordered from dof order to vertex order."""
    A = np.asarray(op.weak_form().A)
    inv = np.empty_like(vmap)
    inv[vmap] = np.arange(len(vmap))              # dof -> vertex
    # A_dof[i,j]; vertex matrix B[v,w] = A_dof[vmap[v], vmap[w]]
    return A[np.ix_(vmap, vmap)]


# ----------------------------------------------------------------------------- closed body
def closed_neumann_bm(grid, vn, k, alpha=None, rho0=RHO0, c0=C0, return_ops=False):
    """Exterior Neumann problem, Burton-Miller (direct), P1 for p and for dp/dn.

    Exterior representation (n outward from body):  p = D[p] - S[q],  q = dp/dn = i w rho0 vn.
      (1/2 - K) p + S q = 0                       (conventional BIE)
      and normal derivative:  q = dD/dn[p] - dS/dn[q],  dD/dn[p] = s_H * H p,
                              dS/dn (exterior) = (K' - 1/2)
      => s_H H p + (1/2 - ... ) handled below; s_H fixed by calibration (see SH).
    Burton-Miller:  [(1/2 - K) + alpha*(-s_H*H)... ]  (assembled explicitly in code).
    """
    omega = k * c0
    space, vmap = full_p1(grid)
    M = mass_matrix_vertex(grid)
    S = _op_vertex_order(bem.operators.boundary.helmholtz.single_layer(space, space, space, k), vmap)
    K = _op_vertex_order(bem.operators.boundary.helmholtz.double_layer(space, space, space, k), vmap)
    Kp = _op_vertex_order(bem.operators.boundary.helmholtz.adjoint_double_layer(space, space, space, k), vmap)
    H = _op_vertex_order(bem.operators.boundary.helmholtz.hypersingular(space, space, space, k), vmap)
    Md = M.toarray()
    q = 1j * omega * rho0 * np.asarray(vn, dtype=complex)
    if alpha is None:
        alpha = 1j / (k + 1.0)                    # standard coupling choice ~ i/k, bounded at small k
    # Conventional:  (1/2 M - K) p = -S q
    # Derivative  :  q = dD/dn p - (K' - 1/2 M) q  with dD/dn p = SH * H p   ->  SH*H p = (1/2 M + K') q  (weak, see below)
    SH = SIGN_H
    A1 = 0.5 * Md - K
    b1 = -S @ q
    A2 = SH * H                                    # SH*H p = (1/2 M + K') q
    b2 = (0.5 * Md + Kp) @ q
    A = A1 + alpha * A2
    b = b1 + alpha * b2
    p = sla.solve(A, b)
    P = 0.5 * np.real(np.conj(np.asarray(vn, dtype=complex)) @ (M @ p))
    if return_ops:
        return P, p, dict(S=S, K=K, Kp=Kp, H=H, M=M)
    return P, p


SIGN_H = -1.0  # CALIBRATED (exact point-source test): dD/dn[p] = -H p for bempp-cl hypersingular H


def _sphere_dirs(nth, nph):
    x, w = np.polynomial.legendre.leggauss(nth)
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    ct = x[:, None] * np.ones((1, nph))
    st = np.sqrt(1 - ct**2)
    dirs = np.stack([st * np.cos(ph)[None], st * np.sin(ph)[None], ct], axis=-1).reshape(-1, 3)
    W = (w[:, None] * np.ones((1, nph)) * (2 * np.pi / nph)).ravel()
    return dirs, W


def far_field_power_dl(grid, mu, k, nth=None, nph=None, rho0=RHO0, c0=C0):
    """Far-field radiated power of p = D[mu] (mu given as P1 vertex values)."""
    space, vmap = full_p1(grid)
    L = max(np.ptp(grid.vertices, axis=1))
    n = int(max(48, 4 * k * L + 24))
    dirs, W = _sphere_dirs(nth or n, nph or 2 * n)
    dof = np.zeros(len(vmap), dtype=complex)
    dof[vmap] = mu
    pinf = bem.operators.far_field.helmholtz.double_layer(space, dirs.T, k).evaluate(
        bem.GridFunction(space, coefficients=dof))
    pinf = np.asarray(pinf).ravel()
    return float(np.sum(W * np.abs(pinf) ** 2) / (2 * rho0 * c0))


def far_field_power_closed(grid, p, q, k, nth=None, nph=None, rho0=RHO0, c0=C0):
    """Far-field power for the exterior field p_field = D[p] - S[q] (p, q P1 vertex values)."""
    space, vmap = full_p1(grid)
    L = max(np.ptp(grid.vertices, axis=1))
    n = int(max(48, 4 * k * L + 24))
    dirs, W = _sphere_dirs(nth or n, nph or 2 * n)
    dp = np.zeros(len(vmap), dtype=complex); dp[vmap] = p
    dq = np.zeros(len(vmap), dtype=complex); dq[vmap] = q
    ffd = bem.operators.far_field.helmholtz.double_layer(space, dirs.T, k)
    ffs = bem.operators.far_field.helmholtz.single_layer(space, dirs.T, k)
    pinf = (np.asarray(ffd.evaluate(bem.GridFunction(space, coefficients=dp))).ravel()
            - np.asarray(ffs.evaluate(bem.GridFunction(space, coefficients=dq))).ravel())
    return float(np.sum(W * np.abs(pinf) ** 2) / (2 * rho0 * c0))


# ----------------------------------------------------------------------------- open screen
def open_screen_power(grid, vn, k, rho0=RHO0, c0=C0, return_ops=False):
    """Two-sided open screen, same v_n on both faces; unknown = pressure jump mu = p(+n) - p(-n).

    p = D[mu] (double-layer potential of the jump); Neumann condition on the screen
        dp/dn = i w rho0 vn   <=>   SIGN_H * H mu = i w rho0 vn   (weak form, tested with P1).
    mu has dofs ONLY at interior vertices (boundary-vertex mu = 0: edge condition).
    Returns P = 0.5*Re(int mu conj(vn) dS)*SIGN_P  (SIGN_P fixed by calibration so P>0), mu (vertex values).
    """
    omega = k * c0
    space, vmap = full_p1(grid)
    M = mass_matrix_vertex(grid)
    H = _op_vertex_order(bem.operators.boundary.helmholtz.hypersingular(space, space, space, k), vmap)
    bnd = boundary_vertices(grid)
    ii = np.where(~bnd)[0]
    vnc = np.asarray(vn, dtype=complex)
    rhs = 1j * omega * rho0 * (M @ vnc)
    mu = np.zeros(len(vnc), dtype=complex)
    mu[ii] = sla.solve(SIGN_SCREEN * H[np.ix_(ii, ii)], rhs[ii])
    P = SIGN_P * 0.5 * np.real(np.conj(vnc) @ (M @ mu))
    if return_ops:
        return P, mu, dict(H=H, M=M, ii=ii)
    return P, mu


SIGN_SCREEN = -1.0   # p = D[mu]  =>  dp/dn = dD/dn[mu] = -H mu  (SIGN_H calibrated above)  =>  -H mu = i w rho0 vn
SIGN_P = 1.0         # power into the fluid = 0.5 Re int (p+ - p-) conj(vn) dS = 0.5 Re int mu conj(vn) dS


# ----------------------------------------------------------------------------- sinc model
_TRI3 = (np.array([[1 / 6, 1 / 6, 2 / 3], [1 / 6, 2 / 3, 1 / 6], [2 / 3, 1 / 6, 1 / 6]]), np.ones(3) / 3)


def surface_samples(grid, vn, rule=_TRI3):
    """Quadrature points x_q, weights w_q (areas included) and v at the points."""
    V = grid.vertices.T
    E = grid.elements.T
    a = np.linalg.norm(np.cross(V[E[:, 1]] - V[E[:, 0]], V[E[:, 2]] - V[E[:, 0]]), axis=1) / 2
    bary, wq = rule
    X = np.einsum("qb,ebd->eqd", bary, V[E])              # (n_el, nq, 3)
    v = np.einsum("qb,eb->eq", bary, np.asarray(vn)[E])
    W = a[:, None] * wq[None]
    return X.reshape(-1, 3), W.ravel(), v.ravel()


def sinc_power_mesh(grid, vn, k, rho0=RHO0, c0=C0, chunk=2000):
    """P_sinc = alpha * sum_ij w_i w_j v_i conj(v_j) sinc(k|x_i - x_j|),  alpha = rho0 w^2/(4 pi c0)."""
    omega = k * c0
    alpha = rho0 * omega**2 / (4 * np.pi * c0)
    X, W, v = surface_samples(grid, np.asarray(vn, dtype=complex))
    wv = W * v
    tot = 0.0 + 0.0j
    for s in range(0, len(X), chunk):
        d = np.linalg.norm(X[s:s + chunk, None, :] - X[None, :, :], axis=2)
        tot += np.sum(wv[s:s + chunk, None] * np.sinc(k * d / np.pi) * np.conj(wv)[None, :])
    return float(alpha * np.real(tot))


# ----------------------------------------------------------------------------- meshes
def disc_grid(a, nr):
    """Structured-ring triangulation of a circular disc (radius a), normal +z, nr rings."""
    pts = [(0.0, 0.0)]
    for i in range(1, nr + 1):
        m = 6 * i
        r = a * i / nr
        th = 2 * np.pi * np.arange(m) / m
        pts += list(zip(r * np.cos(th), r * np.sin(th)))
    pts = np.array(pts)
    from scipy.spatial import Delaunay
    tri = Delaunay(pts).simplices
    # drop slivers on the convex boundary (none expected for ring layout, but keep area exact)
    P = pts[tri]
    cr = (P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 1, 1] - P[:, 0, 1]) * (P[:, 2, 0] - P[:, 0, 0])
    tri = np.where((cr > 0)[:, None], tri, tri[:, [0, 2, 1]])   # counter-clockwise => normal +z
    area = np.abs(cr).sum() / 2
    verts = np.column_stack([pts, np.zeros(len(pts))]).T
    g = bem.Grid(np.ascontiguousarray(verts), np.ascontiguousarray(tri.T.astype(np.uint32)))
    return g
