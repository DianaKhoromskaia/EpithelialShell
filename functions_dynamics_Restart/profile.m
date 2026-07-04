function y = sigmoidal(x,x0,sigma)
%sigmoidal with half-maximum at x0, width sigma
y = 1-(1+exp(-(x-x0)./sigma)).^(-1);
end

