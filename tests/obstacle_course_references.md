# Classical shell obstacle course — definitions, reference values, citations

Geometry + BC + load + target definitions: `tests/obstacle_course_geom.m`
(`problems = obstacle_course_geom()` returns all three as NURBS mid-surfaces, rdim = 3).

## The benchmark set

The "shell obstacle course" is the three-problem set popularised by

- **T. Belytschko, B. L. Wong, H. Stolarski**, *Assumed strain stabilization
  procedure for the 9-node Lagrange shell element*, **Int. J. Numer. Methods
  Engrg. 28 (1989) 385–414** — defines the obstacle course and its target
  values. (An earlier statement of the same set is Belytschko, Stolarski,
  Liu, Carpenter & Ong, CMAME 51 (1985) 221–258.)

and analysed in depth by

- **D. Chapelle, K.-J. Bathe**, *Fundamental considerations for the finite
  element analysis of shell structures*, Computers & Structures 66(1) (1998)
  19–36; and the book **The Finite Element Analysis of Shells – Fundamentals**,
  Springer 2003 (2nd ed. 2011).

The three problems and their load-bearing categories (Chapelle & Bathe):

| # | Problem | Category |
|---|---------|----------|
| 1 | Pinched cylinder with end diaphragms | bending-dominated / mixed (inextensional bending + complex membrane state) |
| 2 | Pinched hemisphere with an 18° hole | inextensional bending (tests rigid-body-motion capability) |
| 3 | Scordelis–Lo roof | membrane-dominated |

## Reference values (confirmed from primary sources)

### 1. Pinched cylinder with end diaphragms
- Geometry: closed circular cylinder, **R = 300, L = 600, t = 3**.
- Material: **E = 3.0e6, ν = 0.3**.
- BCs: rigid diaphragms on both end circles (in-plane/transverse displacements
  zero, axial free). Canonically two radially **opposing unit point loads
  P = 1** at mid-length (z = L/2). Eighth model by symmetry.
- Target: radial displacement under the load point.
- **Reference value: u = 1.8248e-5.**
  - Normalisation: exact `w = 164.24 · P/(E·t)` = 0.18248e-4.
  - Sources: Belytschko/Wong/Stolarski 1989 (p. 404); Bouclier, Elguedj &
    Combescure, *Efficient isogeometric NURBS-based solid-shell elements*,
    CMAME **267 (2013) 86–110** (quoted at 1.8248e-5); Kratos wiki (citing
    Bouclier et al.); Kiendl et al. bending-strip paper, CMAME 199 (2010)
    2403–2416; Hosseini et al. arXiv:1803.02721 (ref. 1.8248e-5);
    Dvorkin & Bathe 1984 / Flügge series (Abaqus benchmark) give 0.1825e-4.

### 2. Pinched hemisphere with an 18° circular hole
- Geometry: hemisphere **R = 10, t = 0.04**, pole cut to a polar angle of 18°
  (hole-edge radius R·sin18° = 3.0902). Equator otherwise free.
- Material: **E = 6.825e7, ν = 0.3**.
- Load: four equal radial point loads at the four cardinal equator points,
  alternating inward/outward. Canonical magnitude |P| = 1 (some authors use 2).
  Quarter model by symmetry.
- Target: radial displacement under a load point.
- **Reference value: |u| = 0.0938** (MacNeal & Harder 1985).
  - **Range across sources / normalisations:** 0.0938 (MacNeal & Harder 1985);
    **0.093** — Simo, Fox & Rifai, CMAME 73 (1989) 53–92 analytic value, used
    e.g. by the HAL prismatic solid-shell paper; **0.094** — used by the
    ORBiLu IGA Reissner–Mindlin paper (u_ref = 0.0940 m); 0.0924 — used by the
    Element-Free-Galerkin paper (Krysl & Belytschko) to normalise.
    The number depends on which normalisation/analytic value the author adopts;
    the obstacle-course value is **0.0938**.

### 3. Scordelis–Lo roof
- Geometry: cylindrical section, **R = 25, L = 50, half-angle 40° (total
  opening 80°), t = 0.25**.
- Material: **E = 4.32e8, ν = 0.0**.
- BCs: rigid diaphragms on the two end arcs (in-plane displacements zero);
  the two straight longitudinal (free) edges are free. Quarter model by symmetry.
- Load: self-weight, uniform vertical load **g = 90 per unit area** (downward).
- Target: vertical displacement u_z at the **midpoint of a free edge**
  (physical point (R·cos50°, L/2, R·sin50°) = (16.0697, 25, 19.1511) in the
  GeoPDEs orientation).
- **Reference value: u_z = 0.3024** (downward).
  - Sources: Belytschko/Wong/Stolarski 1989; **Lee & Bathe**, *On the
    asymptotic behavior of shell structures…*, Computers & Structures 80
    (2002) 1441–1455 — "0.3024 when t = 0.25"; **Ko, Lee & Bathe**, Computers
    & Structures 193 (2017) 187–206 — reference −0.3024; RWTH preprint
    (uref = 0.3024, refs 29,30); arXiv:1812.05596 (0.3024).
  - Provenance caveat: **Sanz et al., *Benchmarking Computational Shell
    Models*, Arch. Comput. Methods Eng. (2022), doi:10.1007/s11831-022-09798-5**
    traces 0.3024 to Belytschko et al. (1985) (superseding Carr's 0.3080/0.3086)
    and notes it is often mis-attributed. Note also that the GeoPDEs KL-shell
    example `ex_KL_shell_scordelis_lo_roof.m` quotes its own converged NURBS
    value −0.3005925, close to but not identical with 0.3024.

## IGA treatments of the same problems

- **T. J. R. Hughes, J. A. Cottrell, Y. Bazilevs**, *Isogeometric analysis:
  CAD, finite elements, NURBS, exact geometry and mesh refinement*, CMAME
  **194 (2005) 4135–4195**, §4.5 "Shell obstacle course" — all three problems,
  convergence compared to the Belytschko et al. benchmarks (NURBS elements
  converge to a "very slightly softer" solution for the pinched cylinder).
- **J. Kiendl, K.-U. Bletzinger, J. Linhard, R. Wüchner**, *Isogeometric shell
  analysis with Kirchhoff–Love elements*, CMAME **198 (2009) 3902–3914** —
  Scordelis–Lo roof, pinched cylinder, pinched hemisphere.
- **J. Kiendl, Y. Bazilevs, M.-C. Hsu, R. Wüchner, K.-U. Bletzinger**, *The
  bending strip method for IGA of Kirchhoff–Love shells…*, CMAME **199 (2010)
  2403–2416** — pinched cylinder (F = 1, ends on rigid diaphragms), hemispherical shell.
- **D. J. Benson, Y. Bazilevs, M.-C. Hsu, T. J. R. Hughes**, *Isogeometric
  shell analysis: the Reissner–Mindlin shell*, CMAME **199 (2010) 276–289** —
  obstacle course with a Reissner–Mindlin IGA shell.
- **R. Bouclier, T. Elguedj, A. Combescure**, CMAME **267 (2013) 86–110** —
  NURBS solid-shell; quotes cylinder 1.8248e-5 and roof 0.3024.

Reported reference solutions used by these IGA papers coincide with the values
above: cylinder **1.8248e-5**, hemisphere **0.0938** (or 0.094), roof **0.3024**.

## Interface the rdim=3 analysis operator needs from this geometry

The current pipeline (`harmonic_kd_setup_curved` + `quadrule_bend_curved`) is
**2D only** (`msh.ndim == 2`, and `quadrule_bend_curved` reshapes the map
Jacobian as `2×2` and the Hessian as `2×2×2`). The geometry supplied here is a
**surface in R³** (`rdim = 3`, `ndim = 2`). The rdim=3 operator must consume:

- `geometry.map`      : `[3 × npts]` physical mid-surface point x(ξ,η);
- `geometry.map_der`  : `[3 × 2 × npts]` tangent vectors a₁ = ∂x/∂ξ, a₂ = ∂x/∂η;
- `geometry.map_der2` : `[3 × 2 × 2 × npts]` second derivatives ∂²x/∂ξᵅ∂ξᵝ
  (needed for the second fundamental form);
- `geometry.rdim = 3`, `msh.ndim = 2`, `geometry.boundary` (from `geo_load`);

and must form, per quadrature/WQ point, the shell kinematics that replace the
flat 2×2 Jacobian/Hessian of the 2D path:

- surface metric `a_{αβ} = a_α · a_β`, area element `sqrt(det(a))` = |a₁ × a₂|;
- unit normal `a₃ = (a₁ × a₂)/|a₁ × a₂|`;
- second fundamental form `b_{αβ} = a₃ · ∂²x/∂ξᵅ∂ξᵝ` (Christoffel symbols for
  covariant surface derivatives);
- the Kirchhoff–Love (or T1) bending/membrane coefficient page expressed in
  these surface quantities, contracted into the same WQ rule slots (`(k1,k2)`
  in {0,1,2}) that `quadrule_bend_fast` already produces.

The load/BC assembly must likewise be 3-component: point forces/body forces as
`[3 × N]`, and Dirichlet sides constrain the listed global components (1=x, 2=y,
3=z), exactly as specified in `obstacle_course_geom.m`.

## Which problems the current pipeline can express

**None of the three, in the strict sense.** All three mid-surfaces are genuine
3D shells (`rdim = 3`, `ndim = 2`): the 2D operator assumes a 2D physical chart
(R² → R²). Verified failure: on the Scordelis–Lo NURBS surface,
`harmonic_kd_setup_curved` aborts inside `quadrule_bend_curved` at line 79
(`reshape(Jall,2,2,nq1,nq2)` applied to a `3×2` Jacobian). **All three are
blocked on the rdim=3 surface operator.** The geometries, BCs, loads and target
quantities are built, validated (every target point lies on its surface to
machine precision) and ready to be analysed as soon as that operator lands.
