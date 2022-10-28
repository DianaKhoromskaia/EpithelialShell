function x = intensity(t,th,w)
%sigmoidal with half-maximum at th, width 1/w
x = 1-(1+exp(-w.*(t-th))).^(-1);
end

