function y = yguessfun_initial(s, sol, factor)
 y = deval(sol, s);
 
 y(1:8) = factor*y(1:8);
 
end

