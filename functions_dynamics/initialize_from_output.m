function [C1,C2,C,dsC1,dsC2,dsC,Psi,X,Z,X0,xintegral,svec1, ...
          U,dsU,Q,dsQ,dsw0,dswL,fs,s0,Qgrid] = initialize_from_output(outdir,npoints)
% Initialize state from the LAST timestep saved in 'outdir'.
% Reads each .dat file individually; uses only the last non-empty line.

req = @(f) fullfile(outdir,f);
has = @(f) exist(req(f),'file')>0;

% --- 1) last-row meridian (radius/axis) -------------------------
xr = read_last_row(req('x.dat'));   % 1×Nx (ragged is fine)
zr = read_last_row(req('z.dat'));   % 1×Nz

% pairwise finite trimming
n  = min(numel(xr), numel(zr));
xr = xr(1:n); zr = zr(1:n);
m  = isfinite(xr) & isfinite(zr);
xr = xr(m);    zr = zr(m);

if numel(xr) < 2
    error('Last row in x.dat/z.dat has < 2 valid samples; cannot interpolate.');
end

% --- 2) matching s-grid for that last row (optional) ------------
sr = [];
if has('svecnew.dat')
    sr = read_last_row(req('svecnew.dat'));
elseif has('s0.dat')
    sr = read_last_row(req('s0.dat'));
end

% align lengths if sr exists
if ~isempty(sr)
    k  = min([numel(sr), numel(xr), numel(zr)]);
    sr = sr(1:k); xr = xr(1:k); zr = zr(1:k);
end

% if sr missing/invalid, build from arclength (coarse)
if isempty(sr) || numel(sr)<2 || any(~isfinite(sr)) || any(diff(sr)<=0)
    dsarc = hypot(diff(xr), diff(zr));
    sr = [0, cumsum(dsarc)];
end

% --- 3) optional payloads from the same last timestep -----------
psi_row = [];
if has('psi.dat')
    psi_row = read_last_row(req('psi.dat'));   % may be ragged length
end

% --- 4) sanitize: enforce monotone-unique sr and apply to payloads
if ~isempty(psi_row) && numel(psi_row) >= 2
    [sr, xr, zr, psi_row] = sanitize_samples(sr, xr, zr, psi_row);
else
    [sr, xr, zr] = sanitize_samples(sr, xr, zr);
    psi_row = [];  % force recompute from geometry below
end

% --- 5) resample to solver grid ---------------------------------
s   = linspace(0, sr(end), npoints).';
pce = @(y) interp1(sr, y, s, 'pchip','extrap');

X = pce(xr);
Z = pce(zr);

% Psi (tangent angle): from file if available, else from geometry
if ~isempty(psi_row)
    Psi = unwrap(pce(psi_row));             % geometric angle, unwrap after interp
else
    dX = gradient(X, s);
    dZ = gradient(Z, s);
    tnorm = hypot(dX, dZ);  tnorm(tnorm==0)=1;
    tx = dX ./ tnorm;  tz = dZ ./ tnorm;
    Psi = unwrap(atan2(tz, tx));
end

% enforce endpoint convention
Psi(1)   = 0;
Psi(end) = pi;

% --- 6) curvatures from geometry --------------------------------
dX = gradient(X, s);   dZ = gradient(Z, s);
tn = hypot(dX,dZ);     tn(tn==0)=1;
tx = dX./tn;           tz = dZ./tn;

psi_geom  = atan2(tz, tx);
kappa_s   = gradient(psi_geom, s);           % meridional curvature
r         = max(1e-9, abs(X));
kappa_phi = sin(psi_geom)./r;                % azimuthal curvature (axisymmetric)

C1 = kappa_s;
C2 = kappa_phi;
C  = C1 + C2;                                % match internal convention

% derivatives & auxiliaries used later
dsC1 = gradient(C1, s);
dsC2 = gradient(C2, s);
dsC  = gradient(C , s);

X0        = X(1);
xintegral = trapz(s, X);

% --- 7) function-like fields expected elsewhere -----------------
svec1 = s; 
Qgrid = s;

% identity material map and zero profiles as griddedInterpolants
s0  = griddedInterpolant(svec1, svec1, 'linear','nearest');   % identity: s0(s) = s
U   = griddedInterpolant(svec1, zeros(size(svec1)), 'linear','nearest');
dsU = griddedInterpolant(svec1, zeros(size(svec1)), 'linear','nearest');
Q   = griddedInterpolant(svec1, zeros(size(svec1)), 'linear','nearest');   % Q=0 (per your choice)
dsQ = griddedInterpolant(svec1, zeros(size(svec1)), 'linear','nearest');
fs  = griddedInterpolant(svec1, zeros(size(svec1)), 'linear','nearest');

dsw0 = 0;
dswL = 0;

% --- 8) ensure columns for numeric arrays only ------------------
X   = X(:);  Z = Z(:);  Psi = Psi(:);
C1  = C1(:); C2 = C2(:); C = C(:);
dsC1= dsC1(:); dsC2 = dsC2(:); dsC = dsC(:);
% (do NOT reshape s0, U, dsU, Q, dsQ, fs — they are griddedInterpolant)
end


