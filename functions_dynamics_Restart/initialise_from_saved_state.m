function [C1, C2, C, dsC1, dsC2, dsC, Psi, X, Z, X0, xintegral, svec1, ...
          U, dsU, Q, dsQ, dsw0, dswL, fs, s0, Qgrid, L, L0, V0] = ...
          initialise_from_saved_state(outdir, npoints)
% Load the LAST saved full state from a previous output directory.
%
% Expected files in outdir:
%   x.dat, z.dat, psi.dat, u.dat, q.dat, s0.dat, svecnew.dat
%
% Assumption:
%   each save appends one row per timepoint to those files,
%   and the last non-empty numeric row is the latest state.

    req = @(f) fullfile(outdir, f);
    must_exist(req('x.dat'));
    must_exist(req('z.dat'));
    must_exist(req('psi.dat'));
    must_exist(req('u.dat'));
    must_exist(req('q.dat'));
    must_exist(req('s0.dat'));
    must_exist(req('svecnew.dat'));

    % --- read last saved rows ---
    x_row   = read_last_numeric_row(req('x.dat'));
    z_row   = read_last_numeric_row(req('z.dat'));
    psi_row = read_last_numeric_row(req('psi.dat'));
    u_row   = read_last_numeric_row(req('u.dat'));
    q_row   = read_last_numeric_row(req('q.dat'));
    s0_row  = read_last_numeric_row(req('s0.dat'));
    s_row   = read_last_numeric_row(req('svecnew.dat'));

    % --- trim all arrays to common valid length ---
    n = min([numel(x_row), numel(z_row), numel(psi_row), numel(u_row), ...
             numel(q_row), numel(s0_row), numel(s_row)]);

    x_row   = x_row(1:n);
    z_row   = z_row(1:n);
    psi_row = psi_row(1:n);
    u_row   = u_row(1:n);
    q_row   = q_row(1:n);
    s0_row  = s0_row(1:n);
    s_row   = s_row(1:n);

    mask = isfinite(x_row) & isfinite(z_row) & isfinite(psi_row) & ...
           isfinite(u_row) & isfinite(q_row) & isfinite(s0_row) & isfinite(s_row);

    x_row   = x_row(mask);
    z_row   = z_row(mask);
    psi_row = psi_row(mask);
    u_row   = u_row(mask);
    q_row   = q_row(mask);
    s0_row  = s0_row(mask);
    s_row   = s_row(mask);

    if numel(s_row) < 5
        error('Restart files contain too few valid points.');
    end

    % --- ensure strictly increasing current arc-length grid ---
    [s_row, x_row, z_row, psi_row, u_row, q_row, s0_row] = ...
        sanitize_samples_full(s_row, x_row, z_row, psi_row, u_row, q_row, s0_row);

    % --- restart lengths ---
    L  = s_row(end);
    L0 = s0_row(end);   % assumes s0 stores reference arclength coordinate

    % --- resample to solver grid ---
    svec1 = linspace(0, L, npoints).';

    x_vec   = interp1(s_row, x_row,   svec1, 'pchip', 'extrap');
    z_vec   = interp1(s_row, z_row,   svec1, 'pchip', 'extrap');
    psi_vec = unwrap(interp1(s_row, psi_row, svec1, 'pchip', 'extrap'));
    u_vec   = interp1(s_row, u_row,   svec1, 'pchip', 'extrap');
    q_vec   = interp1(s_row, q_row,   svec1, 'pchip', 'extrap');
    s0_vec  = interp1(s_row, s0_row,  svec1, 'pchip', 'extrap');

    % --- enforce monotone s0 if needed ---
    s0_vec = enforce_monotone(s0_vec);

    % --- geometry as interpolants ---
    X   = griddedInterpolant(svec1, x_vec,   'spline');
    Z   = griddedInterpolant(svec1, z_vec,   'spline');
    Psi = griddedInterpolant(svec1, psi_vec, 'spline');

    % --- mechanics/state as interpolants ---
    U   = griddedInterpolant(svec1, u_vec,  'spline');
    Q   = griddedInterpolant(svec1, q_vec,  'spline');
    s0  = griddedInterpolant(svec1, s0_vec, 'spline');

    ds = L / (npoints - 1);

    du_vec = gradient(u_vec, ds);
    dq_vec = gradient(q_vec, ds);

    dsU = griddedInterpolant(svec1, du_vec, 'spline');
    dsQ = griddedInterpolant(svec1, dq_vec, 'spline');

    % conservative placeholder
    fs = griddedInterpolant(svec1, zeros(size(svec1)), 'linear', 'nearest');

    dsw0 = 0;
    dswL = 0;

    % --- reconstruct curvatures from geometry ---
    dX = gradient(x_vec, ds);
    dZ = gradient(z_vec, ds);
    tnorm = hypot(dX, dZ);
    tnorm(tnorm == 0) = 1;

    tx = dX ./ tnorm;
    tz = dZ ./ tnorm;
    psi_geom = unwrap(atan2(tz, tx));

    if max(abs(psi_geom - psi_vec)) < 0.2
        psi_use = psi_vec;
    else
        psi_use = psi_geom;
        Psi = griddedInterpolant(svec1, psi_use, 'spline');
    end

    C1_vec = gradient(psi_use, ds);
    r = max(1e-9, abs(x_vec));
    C2_vec = sin(psi_use) ./ r;
    C_vec  = C1_vec + C2_vec;

    dsC1_vec = gradient(C1_vec, ds);
    dsC2_vec = gradient(C2_vec, ds);
    dsC_vec  = gradient(C_vec, ds);

    C1   = griddedInterpolant(svec1, C1_vec,   'spline');
    C2   = griddedInterpolant(svec1, C2_vec,   'spline');
    C    = griddedInterpolant(svec1, C_vec,    'spline');
    dsC1 = griddedInterpolant(svec1, dsC1_vec, 'spline');
    dsC2 = griddedInterpolant(svec1, dsC2_vec, 'spline');
    dsC  = griddedInterpolant(svec1, dsC_vec,  'spline');

    % --- integrals / observables ---
    xintegral = trapz(svec1, x_vec);
    X0 = trapz(svec1, x_vec .* z_vec) / max(xintegral, 1e-12);

    dzds = gradient(z_vec, ds);
    V0 = abs(pi * trapz(svec1, x_vec.^2 .* dzds));

    Qgrid = svec1;
end

function must_exist(fname)
    if exist(fname, 'file') ~= 2
        error('Missing restart file: %s', fname);
    end
end

function row = read_last_numeric_row(fname)
    fid = fopen(fname, 'r');
    if fid == -1
        error('Could not open file: %s', fname);
    end

    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    last_good = [];

    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end

        if isempty(strtrim(tline))
            continue;
        end

        vals = sscanf(tline, '%f').';
        if ~isempty(vals)
            last_good = vals;
        end
    end

    if isempty(last_good)
        error('No numeric row found in file: %s', fname);
    end

    row = last_good;
end

function [s, x, z, psi, u, q, s0] = sanitize_samples_full(s, x, z, psi, u, q, s0)

    % force everything to columns first
    s   = s(:);
    x   = x(:);
    z   = z(:);
    psi = psi(:);
    u   = u(:);
    q   = q(:);
    s0  = s0(:);

    % sort by s
    [s, idx] = sort(s);
    x   = x(idx);
    z   = z(idx);
    psi = psi(idx);
    u   = u(idx);
    q   = q(idx);
    s0  = s0(idx);

    % unique by s
    [s, ia] = unique(s, 'stable');
    x   = x(ia);
    z   = z(ia);
    psi = psi(ia);
    u   = u(ia);
    q   = q(ia);
    s0  = s0(ia);

    % keep only finite entries
    good = isfinite(s) & isfinite(x) & isfinite(z) & isfinite(psi) & ...
           isfinite(u) & isfinite(q) & isfinite(s0);

    s   = s(good);
    x   = x(good);
    z   = z(good);
    psi = psi(good);
    u   = u(good);
    q   = q(good);
    s0  = s0(good);

    if numel(s) < 5
        error('Too few valid restart points after sanitization.');
    end

    if any(diff(s) <= 0)
        error('Restart s-grid is not strictly increasing.');
    end
end

function y = enforce_monotone(y)
    y = y(:).';
    for k = 2:numel(y)
        if y(k) <= y(k-1)
            y(k) = y(k-1) + 1e-9;
        end
    end
    y = y(:);
end