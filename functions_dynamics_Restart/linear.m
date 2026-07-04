function y = linear(x, x0, A, Aconst, s, sigma)
%function linear 

    y=Aconst+A*(x-x0+0.5*sign(A)*s).*rect(x,x0,1,0,s,sigma)/s;

end