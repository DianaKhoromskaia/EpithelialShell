function x = dxprofile(t,th,w)
%derivative of sigmoidal with half-maximum at th, width 1/w
x = -w*exp(-w.*(t-th))./(1+exp(-w.*(t-th))).^2;
end