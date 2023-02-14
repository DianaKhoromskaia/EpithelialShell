function y = yguessfun_nematic_t0(s, sfun, q, dsq)
 %y = deval(solnem, sfun(s));
 
y(1) = q(sfun(s));
y(2) = dsq(sfun(s));
 
end

