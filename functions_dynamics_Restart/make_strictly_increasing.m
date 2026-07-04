function x = make_strictly_increasing(x)
    x = x(:).';
    for k = 2:numel(x)
        if x(k) <= x(k-1)
            x(k) = x(k-1) + 1e-12;
        end
    end
end