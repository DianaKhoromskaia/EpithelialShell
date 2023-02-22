function y = yguessfun_nematic_sphere(s,L0)
 %y = deval(solnem, sfun(s));
 
y(1) = 0.5.*(1-cos(pi*s/L0));
y(2) = 0.5*pi*sin(pi*s/L0)/L0;
 
end