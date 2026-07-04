function [sr,xr,zr,varargout] = sanitize_samples(sr,xr,zr,varargin)
% Sanitize a ragged last-row meridian and any number of co-sampled payloads.
% - Ensures sr is strictly increasing & unique (builds from arclength if needed)
% - Drops NaNs consistently across all vectors
% - Removes zero-length segments (duplicate (x,z) points)
% - Sorts by sr and deduplicates within a tolerance
%
% Inputs:  sr ([], or 1×N), xr (1×N), zr (1×N), payload1 (1×N), payload2 (1×N), ...
% Outputs: sr, xr, zr, payload1, payload2, ...

    % --- shape as row vectors
    xr = xr(:).'; zr = zr(:).';
    if ~isempty(sr), sr = sr(:).'; end
    P  = numel(varargin);
    pay = cell(1,P);
    for k=1:P, pay{k} = varargin{k}(:).'; end

    % --- align lengths to the shortest among inputs provided
    lens = [numel(xr), numel(zr)];
    if ~isempty(sr), lens(end+1) = numel(sr); end
    for k=1:P, lens(end+1) = numel(pay{k}); end
    L = min(lens);
    xr = xr(1:L); zr = zr(1:L);
    if ~isempty(sr), sr = sr(1:L); end
    for k=1:P, pay{k} = pay{k}(1:L); end

    % --- finite mask across all vectors
    m = isfinite(xr) & isfinite(zr);
    if ~isempty(sr), m = m & isfinite(sr); end
    for k=1:P, m = m & isfinite(pay{k}); end
    xr = xr(m); zr = zr(m);
    if ~isempty(sr), sr = sr(m); end
    for k=1:P, pay{k} = pay{k}(m); end

    % --- (re)build sr if missing or non-monotone
    if isempty(sr) || numel(sr)<2 || any(diff(sr)<=0)
        dx = diff(xr); dz = diff(zr);
        ds = hypot(dx,dz);
        keep = [true, ds>0];      % remove zero-length segments
        xr = xr(keep); zr = zr(keep);
        for k=1:P, pay{k} = pay{k}(keep); end
        ds = ds(ds>0);
        sr = [0, cumsum(ds)];
    end

    % --- sort by sr and deduplicate within tolerance
    [sr, idx] = sort(sr);
    xr = xr(idx); zr = zr(idx);
    for k=1:P, pay{k} = pay{k}(idx); end

    tol  = max(1e-12, 1e-9*max(sr(end),1));   % length-scale tolerant
    keep = [true, diff(sr) > tol];
    sr = sr(keep); xr = xr(keep); zr = zr(keep);
    for k=1:P, pay{k} = pay{k}(keep); end

    % --- ensure we still have at least 2 points
    if numel(sr) < 2
        error('sanitize_samples:NotEnoughPoints', ...
              'Not enough distinct samples after cleaning.');
    end

    varargout = pay;
end
