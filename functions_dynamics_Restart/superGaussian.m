function y = superGaussian(x, x0, A, Aconst, s, P)
%super-Gaussian function centred at x0, width s, power P, offset Aconst

y = A*exp(-(((x-x0*ones(1,length(x))).*(x-x0*ones(1,length(x))))/(2*s^2)).^P)+Aconst;

end

