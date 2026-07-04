function y = yguessfun(s, sol, sfun)

    sq = sfun(s);

    % clamp to valid domain of sol
    sq = max(min(sq, sol.x(end)), sol.x(1));

    y = deval(sol, sq);

end
