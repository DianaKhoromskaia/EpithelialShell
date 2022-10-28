function [dbcdya, dbcdyb, dbcdp] = bcjac_nematic(ya, yb, varpar, Psi, X, L, lc)
% Jacobian of the residual function specified in bc_nematic.m
% size of matrix = (nBC, nregions x nequations)

dbcdya = zeros(4,2);
dbcdyb = zeros(4,2);
dbcdp = zeros(4,2);

%partial derivatives wrt ya
dbcdya(1,1) = 1;
dbcdya(2,2) = 1;

%partial derivatives wrt yb
dbcdyb(3,1) = 1;
dbcdyb(4,2) = 1;

%partial derivatives wrt parameters

end

