function dvds = ode_nematic(s, v, varpar, Psi, X, L, lc)
% RHS of the system of ODE's representing the Euler-Lagrange equations on
% the displaced surface at time t=t+dt,
% written for the vector v = (q, \partial_s q)

dsw0 = varpar(1);
dswL = varpar(2);
onevec = ones(size(s));

q = v(1,:);
dsq = v(2,:);

% force balance within the interval:
dvds = [    dsq;
    (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(Psi(s))./X(s)).*(4*(cos(Psi(s))./X(s)).*q - dsq)];

% force balance at SP:
indices = find(~s);
zervec = zeros(size(s(indices)));
onevec = ones(size(s(indices)));

dvds(:,indices) =   [   zervec;
    dsw0*onevec];

% force balance at NP:
indices = find(~(s-L));
zervec = zeros(size(s(indices)));
onevec = ones(size(s(indices)));

dvds(:,indices) =   [   zervec;
    dswL*onevec];

end

