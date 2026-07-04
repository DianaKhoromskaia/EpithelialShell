function y = rect(x, x0, A, Aconst, s, sigma)
%function rectangle at x0, width s, height A, offset Aconst approximates thanks to
%sigmoidals of width sigma

    y=A*(sigmoidal(x,x0+0.5*s,sigma)-sigmoidal(x,x0-0.5*s,sigma))+Aconst;

end