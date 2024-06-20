function y = exponential(x,x0,sigma)
%sigmoidal with half-maximum at x0, width sigma
y = exp(-(x-x0)./sigma);
end

