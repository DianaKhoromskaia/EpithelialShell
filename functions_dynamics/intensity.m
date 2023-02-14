function x = intensity(t,th,w)
x = (1+exp(-w.*(t-th))).^(-1);
end

