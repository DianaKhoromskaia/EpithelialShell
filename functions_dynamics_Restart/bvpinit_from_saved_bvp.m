function solinit = bvpinit_from_saved_bvp(bvpfile, parguess)
    [s_last, y_last] = load_last_bvp_block(bvpfile);

    guessfun = @(xq, varargin) eval_saved_bvp_guess(xq, s_last, y_last);

    solinit = bvpinit(s_last, guessfun, parguess);
end

function yq = eval_saved_bvp_guess(xq, s_last, y_last)
    yq = zeros(size(y_last,1), numel(xq));
    for i = 1:size(y_last,1)
        yq(i,:) = interp1(s_last, y_last(i,:), xq, 'pchip', 'extrap');
    end
end