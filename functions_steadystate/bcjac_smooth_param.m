function [dbcdya, dbcdyb, dbcdp] = bcjac_smooth_param(ya, yb, varpar, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar, Profile, width)
% Jacobian of the residual function for solution with smooth profile, in parametric mode
% size of matrix = (nBC, nreg x nequations)

switch FixedPar
    case 'V'
        P = varpar(1);
        V = fixedpar;
    case 'P'
        V = varpar(1);
        P = fixedpar;
end

neq = 12;
nreg = 1;
nparam = 5;
nbc = nreg*neq+nparam;

dbcdya = zeros(nbc,nreg*neq);
dbcdyb = zeros(nbc,nreg*neq);
dbcdp = zeros(nbc,nparam);

%% partial derivatives wrt ya = (ya_firstinterval, ya_secondinterval, ...)

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

dbcdya(16,11) = 1; 

%% partial derivatives wrt yb = (yb_firstinterval, yb_secondinterval, ...)

% north pole:
dbcdyb(10,4) = 1;
dbcdyb(11,5) = 1;
dbcdyb(12,7) = 1;
dbcdyb(13,9) = 1;
dbcdyb(14,10) = 1;
dbcdyb(15,12) = 1;

dbcdyb(17,11) = 1;

%% partial derivatives wrt parameters
switch FixedPar
    case 'V'

    case 'P'
        dbcdp(12,1) = -1; %wrt volume
end

end