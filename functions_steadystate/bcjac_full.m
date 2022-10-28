function [dbcdya, dbcdyb, dbcdp] = bcjac_full(ya, yb, varpar, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar)
% Jacobian of the residual function for solution with homogeneous profile (la=L0), in non-parametric mode
% size of matrix = (nBC, nreg x nequations)

neq = 12;
nreg = 1;
nparam = 5;
nbc = nreg*neq+nparam;

dbcdya = zeros(nbc,nreg*neq);
dbcdyb = zeros(nbc,nreg*neq);
dbcdp = zeros(nbc,nparam);

%% partial derivatives wrt left boundary values ya = (ya_firstinterval, ya_secondinterval, ...)

% south pole:
dbcdya(1,2) = 1;
dbcdya(2,4) = 1;
dbcdya(3,5) = 1;
dbcdya(4,6) = 1;
dbcdya(5,7) = 1;
dbcdya(6,8) = 1;
dbcdya(7,9) = 1;
dbcdya(8,10) = 1;
dbcdya(9,12) = 1;

% dsq = 0 on SP
dbcdya(16,11) = 1;

%% partial derivatives wrt right boundary values yb = (yb_firstinterval, yb_secondinterval, ...)

% north pole:
dbcdyb(10,4) = 1;
dbcdyb(11,5) = 1;
dbcdyb(12,7) = 1;
dbcdyb(13,9) = 1;
dbcdyb(14,10) = 1;
dbcdyb(15,12) = 1;

% dsq = 0 on NP:
dbcdyb(17,11) = 1;

%% partial derivatives wrt parameters
switch FixedPar
    case 'V'

    case 'P'
        dbcdp(12,1) = -1;
end
end