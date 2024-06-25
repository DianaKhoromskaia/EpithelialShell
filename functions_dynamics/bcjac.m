function [dbcdya, dbcdyb, dbcdp] = bcjac(ya, yb, varpar, U, dsU, C1, C2, C, C0, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, eta, etab,etacb, etap, xintegral, fext, kappa, K, xi, kb, x0b, kp, psi0b, FixedPar, t, P0, thalf_P, tsigma)
% Jacobian of the residual function for bvp4c
% size of matrix = (nBC, nreg x nequations)
%xi=1;

dbcdya = zeros(12,10);
dbcdyb = zeros(12,10);
dbcdp = zeros(12,2);

%partial derivatives wrt ya
dbcdya(1,3) = 1;
dbcdya(2,5) = 1;
dbcdya(3,6) = 1;
dbcdya(4,7) = 1;
dbcdya(5,8) = 1;
dbcdya(6,9) = 1;
dbcdya(11,10) = 1;

%partial derivatives wrt yb
dbcdyb(7,4) = 1;
dbcdyb(8,6) = 1;
dbcdyb(9,7) = 1;

dbcdyb(10,1) = eta+etab;
dbcdyb(10,8) = (etab-eta)*cos(Psi(L))./X(L);
dbcdyb(10,2) = etab*C(L)+eta*(C2(L)-C1(L));

%dbcdyb(9,10) = 1/xi;
%dbcdyb(12,10)= 1;
dbcdyb(12,5) = 1;
%dbcdyb(13,10)= 1;
    
%partial derivatives wrt parameters
%dbcdp(9,2) = -xintegral/xi;

if strcmp(FixedPar,'P')
   dbcdp(8,1) = -1;
end

end

