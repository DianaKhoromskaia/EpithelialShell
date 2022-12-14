function [dbcdya, dbcdyb, dbcdp] = bcjac(ya, yb, varpar, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab,etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt)
% Jacobian of the residual function specified in bc.m
% size of matrix = (nBC, nregions x nequations)

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
dbcdyb(7,3) = 1;
dbcdyb(8,6) = 1;
dbcdyb(9,7) = 1;
dbcdyb(10,8) = 1;
dbcdyb(12,5) = 1;
  
%partial derivatives wrt parameters

if strcmp(FixedPar,'P')
   dbcdp(8,1) = -1;
end

end

