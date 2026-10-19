function [X,Y,Z,Nx,Ny,Nz] = nurbs_surface_grid(nrb, U, V)
% NURBS_SURFACE_GRID  Vectorised evaluation of a NURBS surface and its unit
% normal on the tensor grid U x V (ndgrid order: X(i,j) = x(U(j),V(i))).
%
% WHY.  geopdes' geometry.map / geometry.map_der call nrbeval once PER POINT, so
% on a fine display grid (say 200x200 = 40k points) they take tens of seconds.
% Evaluating the rational surface directly through its tensor-product basis is
% milliseconds, and the unit normal is obtained from central differences in
% parameter space (the surface is smooth and the grid is uniform).
%
% nrb : NURBS toolbox surface struct (.coefs homogeneous (dim+1) x (n1*n2),
%       .knots 1x2 cell, .order, .number [n1 n2], .dim).

sz = size(nrb.coefs); nch = sz(1);           % nch = dim+1 (homogeneous coords)
if ndims(nrb.coefs) >= 3
    C = nrb.coefs; n1 = sz(2); n2 = sz(3);
else
    n1 = nrb.number(1); n2 = nrb.number(2);
    C = reshape(nrb.coefs, nch, n1, n2);
end
dim = nch - 1;
sp1 = struct('knots', nrb.knots{1}, 'degree', nrb.order(1)-1, 'ndof', n1);
sp2 = struct('knots', nrb.knots{2}, 'degree', nrb.order(2)-1, 'ndof', n2);
Nu = basis_row(sp1, U);            % nU x n1
Nv = basis_row(sp2, V);            % nV x n2
W  = reshape(C(nch,:,:), n1, n2);            % weights
den = Nv * W.' * Nu.';                       % nV x nU
X = (Nv * reshape(C(1,:,:),n1,n2).' * Nu.') ./ den;
Y = (Nv * reshape(C(2,:,:),n1,n2).' * Nu.') ./ den;
if dim >= 3
    Z = (Nv * reshape(C(3,:,:),n1,n2).' * Nu.') ./ den;
else
    Z = zeros(size(X));
end
% unit normal from central differences in parameter space
hx = U(2)-U(1); hy = V(2)-V(1);
dXx = dpar(X,hx,1); dYx = dpar(Y,hx,1); dZx = dpar(Z,hx,1);
dXy = dpar(X,hy,2); dYy = dpar(Y,hy,2); dZy = dpar(Z,hy,2);
c1 = [dXx(:).'; dYx(:).'; dZx(:).'];
c2 = [dXy(:).'; dYy(:).'; dZy(:).'];
nn = cross(c1, c2, 1);
nn = nn ./ max(vecnorm(nn,2,1), eps);
nrow = numel(V); ncol = numel(U);
Nx = reshape(nn(1,:), nrow, ncol);
Ny = reshape(nn(2,:), nrow, ncol);
Nz = reshape(nn(3,:), nrow, ncol);
end

% ---------------------------------------------------------------------------
function D = dpar(A, h, dir)
% First derivative of A along columns (dir=1) or rows (dir=2), uniform step h.
D = zeros(size(A));
n = size(A, dir);
if n < 2, return; end
if dir == 1
    D(:,2:end-1) = (A(:,3:end) - A(:,1:end-2)) / (2*h);
    D(:,1)   = (A(:,2) - A(:,1)) / h;
    D(:,end) = (A(:,end) - A(:,end-1)) / h;
else
    D(2:end-1,:) = (A(3:end,:) - A(1:end-2,:)) / (2*h);
    D(1,:)   = (A(2,:) - A(1,:)) / h;
    D(end,:) = (A(end,:) - A(end-1,:)) / h;
end
end

% ---------------------------------------------------------------------------
function N = basis_row(spu, q)
kn = spu.knots(:).'; p = spu.degree; nbf = spu.ndof;
q = q(:); nq = numel(q); N = zeros(nq, nbf);
for i = 1:nq
    spn = findspan(nbf-1, p, q(i), kn);
    idx = numbasisfun(spn, q(i), p, kn) + 1;
    N(i, idx) = basisfun(spn, q(i), p, kn);
end
end
