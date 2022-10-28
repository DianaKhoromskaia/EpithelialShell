function [dbcdya, dbcdyb, dbcdp] = bcjac(ya, yb, varpar, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar)
% Jacobian of the residual function for step profile, non-parametric mode
% size of matrix = (nBC, nreg x nequations)

neq = 12;
nreg = 2;
nparam = 6;
nbc = nreg*neq + nparam;

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

% matching boundary:
dbcdya(10,neq+1) = 1;
dbcdya(11,neq+2) = 1;
dbcdya(12,neq+3) = 1;
dbcdya(12,neq+10) = 0.25*zetacnemred;
dbcdya(13,neq+4) = 1;
dbcdya(14,neq+5) = 1;
dbcdya(15,neq+6) = 1;
dbcdya(16,neq+7) = 1;
dbcdya(17,neq+8) = 1;
dbcdya(18,neq+9) = 1;
dbcdya(19,neq+10) = 1;
dbcdya(20,neq+11) = 1;
dbcdya(21,neq+12) = 1;

% dsq = 0 on SP
dbcdya(29,11) = 1;

%% partial derivatives wrt right boundary values yb = (yb_firstinterval, yb_secondinterval, ...)

% matching boundary:
dbcdyb(10,1) = -1;
dbcdyb(11,2) = -1;
dbcdyb(12,3) = -1;
dbcdyb(12,10) = 0.25*zetacnemred;
dbcdyb(13,4) = -1;
dbcdyb(14,5) = -1;
dbcdyb(15,6) = -1;
dbcdyb(16,7) = -1;
dbcdyb(17,8) = -1;
dbcdyb(18,9) = -1;
dbcdyb(19,10) = -1;
dbcdyb(20,11) = -1;
dbcdyb(21,12) = -1;

% north pole:
dbcdyb(22,neq+4) = 1;
dbcdyb(23,neq+5) = 1;
dbcdyb(24,neq+7) = 1;
dbcdyb(25,neq+9) = 1;
dbcdyb(26,neq+10) = 1;
dbcdyb(27,neq+12) = 1;

%boundary at so=la
dbcdyb(28,9) = 1;

% dsq = 0 on NP:
dbcdyb(30,neq+11) = 1;

%% partial derivatives wrt parameters
switch FixedPar
    case 'V'

    case 'P'
        dbcdp(24,1) = -1; %wrt volume
end
end