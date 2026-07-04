function [s_last, y_last] = load_last_bvp_block(fname)
% Load the last saved BVP solution block from bvpsolution.dat
% Each block is assumed to contain:
%   1 row: s grid
%  10 rows: sol.y

    fid = fopen(fname, 'r');
    if fid == -1
        error('Could not open file: %s', fname);
    end
    cleaner = onCleanup(@() fclose(fid));

    rows = {};
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
            rows{end+1} = vals; %#ok<AGROW>
        end
    end

    if isempty(rows)
        error('File is empty: %s', fname);
    end

    nrows = numel(rows);
    blocksize = 11;

    if mod(nrows, blocksize) ~= 0
        error('Unexpected number of numeric rows in %s: got %d, expected multiple of %d.', ...
              fname, nrows, blocksize);
    end

    block = rows(end-blocksize+1:end);

    s_last = block{1};
    y_last = zeros(10, numel(s_last));

    for k = 1:10
        rowk = block{k+1};
        nk = min(numel(rowk), numel(s_last));
        y_last(k,1:nk) = rowk(1:nk);
        if nk < numel(s_last)
            y_last(k,nk+1:end) = NaN;
        end
    end

    good = isfinite(s_last);
    for k = 1:10
        good = good & isfinite(y_last(k,:));
    end

    s_last = s_last(good);
    y_last = y_last(:,good);
    
    [s_last, ia] = unique(s_last, 'stable');
    y_last = y_last(:, ia);
    
    s_last = make_strictly_increasing(s_last);
    s_last = s_last(:).';
    
    if numel(s_last) < 5
        error('Last BVP block in %s has too few valid points.', fname);
    end
end