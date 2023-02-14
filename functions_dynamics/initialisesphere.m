function [C1, C2, C, dsC1, dsC, Psi, X, Z, X0, xintegral, svec1, U, dsU, Q, dsQ, dsw0, dswL, s0] = initialisesphere(L0, R0, z0, npoints)
svec1 = linspace(0, L0, npoints);
sone = ones(1, npoints);
c1vec = sone/R0;   %C_{\phi}^{\phi} = sin(psi)/x
c2vec = c1vec;          %C_s^s = \partial_s psi
psivec = svec1/R0;
xvec = R0*sin(svec1/R0);
zvec = z0*sone+R0*(sone-cos(svec1/R0));

Qinit = dlmread('Qinit.dat','\t');
dsw0 = Qinit(1,1);
dswL = dsw0;
Q = griddedInterpolant(L0*Qinit(2,:),Qinit(3,:),'spline' );
dsQ = griddedInterpolant(L0*Qinit(2,:),Qinit(4,:), 'spline');

C1 = griddedInterpolant(svec1, c1vec, 'linear');
C2 = griddedInterpolant(svec1, c2vec, 'linear');
C =  griddedInterpolant(svec1, 2*c1vec, 'linear');
dsC1 = griddedInterpolant(svec1, zeros(npoints,1),'linear');
dsC = griddedInterpolant(svec1, zeros(npoints,1),'linear');
Psi = griddedInterpolant(svec1, psivec, 'spline');
X = griddedInterpolant(svec1, xvec, 'spline');
Z = griddedInterpolant(svec1, zvec, 'spline');

U = griddedInterpolant(svec1, 0*svec1, 'linear');
dsU = griddedInterpolant(svec1, 0*svec1, 'linear');
s0 = griddedInterpolant(svec1, svec1, 'linear');

%% calculate centre of shape:

xintegral = 2*R0^2;  
X0 = z0+R0; 

end

