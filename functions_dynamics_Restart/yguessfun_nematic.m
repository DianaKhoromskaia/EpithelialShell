function y = yguessfun_nematic(s, sfun, solnem)
    sq = sfun(s);
    
    % clamp to valid interval
    sq = max(min(sq, solnem.x(end)), solnem.x(1));
    
    y = deval(solnem, sq);
     
    %y(1) = solnem.y(sfun(s));
    %y(2) = dsq(sfun(s))/zetanemconst;
 
end

