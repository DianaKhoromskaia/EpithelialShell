function v = read_last_row(fname)
% Return the numeric vector on the LAST non-empty line of a whitespace-delimited text file.
% No padding/rectangular assumptions; supports ragged rows.

fid = fopen(fname,'r');
if fid < 0
    error('Cannot open file: %s', fname);
end
cleanup = onCleanup(@() fclose(fid));

v = [];
last = '';
while true
    t = fgetl(fid);
    if ~ischar(t), break; end
    if ~isempty(strtrim(t))
        last = t;  % keep only the last non-empty line
    end
end

if isempty(last)
    error('File %s has no non-empty lines.', fname);
end

v = sscanf(last, '%f').';
% strip NaNs/Infs if any slipped in
v = v(isfinite(v));
end
