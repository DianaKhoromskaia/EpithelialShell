function [C1, C2, C, dsC1, dsC, Psi, X, Z, X0, xintegral, svec1, U, dsU, Q, dsQ, dsw0, s0] = initialisesphere(L0, R0, z0, npoints, lc, optode_nem)
svec1 = linspace(0, L0, npoints);
sone = ones(1, npoints);
c1vec = sone/R0;   %C_{\phi}^{\phi} = sin(psi)/x
c2vec = c1vec;          %C_s^s = \partial_s psi
psivec = svec1/R0;
xvec = R0*sin(svec1/R0);
zvec = z0*sone+R0*(sone-cos(svec1/R0));

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

% Solve Franck-Landau-de Gennes equation for nematic on initial shape

sgrid_nem=linspace(0, L0, 600);

parguess = 0.5*(pi/L0).^2;
solinit_nem = bvpinit(sgrid_nem, @yguessfun_nematic_sphere, parguess, L0);
   
solnem = bvp4c( @ode_nematic, ...
        @bc_nematic, ...
        solinit_nem, optode_nem, Psi, X, L0, lc);

sgridnem=solnem.x;
Vint=solnem.y;
Vpint=solnem.yp;

dsw0 = Vpint(2,1);
%dswL = dsw0;
Q = griddedInterpolant(sgridnem,Vint(1,:),'spline' );
dsQ = griddedInterpolant(sgridnem,Vint(2,:), 'spline');

% figure()
% plot(sgrid_nem,Q(sgrid_nem))
% 
% figure()
% plot(sgrid_nem,dsQ(sgrid_nem))


%% calculate centre of shape:

xintegral = R0^2*(1-cos(L0/R0));  
X0 = z0+R0^3*(0.75-cos(L0/R0)+0.25*cos(2*L0/R0))/xintegral;

end

