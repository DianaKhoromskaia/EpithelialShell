function y = yguessfun_nematic(s, sfun, solnem)
 y = deval(solnem, sfun(s));
 
%y(1) = solnem.y(sfun(s));
%y(2) = dsq(sfun(s))/zetanemconst;
 
end

