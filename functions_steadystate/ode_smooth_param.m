function dvds = ode_smooth_param(s, v, varpar, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar, Profile, width)
% RHS of ode for smooth profile, parametric mode

%v = (   tss,  tns,    C,  psi,    x,    z, v(s), a(s),   s0,      q,   dsq,    I1,    L,  D)
%  = (  v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9),  v(10), v(11), v(12), v(13), v(14))

switch FixedPar
    case 'V'
        P = varpar(1);
        V = fixedpar;
    case 'P'
        V = varpar(1);
        P = fixedpar;
end
fc = varpar(2);
dsw_0 = varpar(3);
dsw_L = varpar(4);
L = varpar(5);

switch ControlPar
    case 'Nematic'
        zetanem = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'Bending'
        C0 = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'Force'
        fconst = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'Tension'
        zeta = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'BendingNematic'
        zetacnem = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
end

tss = v(1,:);
tns = v(2,:);
mss = v(3,:);
psi = v(4,:);
x = v(5,:);
z = v(6,:);
vol = v(7,:);
area = v(8,:);
s0 = v(9,:);
q = v(10,:);
dsq = v(11,:);
I1 = v(12,:);

onevec = ones(size(s));

% evaluate smooth profiles, on S0:
switch Profile
    case 'Sigmoidal'
        profile = intensity(s0, la, width);
    case 'Gaussian'
        profile = superGaussian(s, 0, -1, 0, la, 1) + superGaussian(s, 1, 1, 0, la, 1);
end

C = 0.5*(mss - C0*profile + zetacnem*profile.*q);
c1 = sin(psi)./x;
c2 = C-c1;

T = zetanem*q.*profile;
t = tss+T;
u = (t - zeta*profile)/(2*K);

fextz = fc*onevec + fconst*profile;
fexts = sin(psi).*fextz;
fextn = -cos(psi).*fextz;

%% RHS in main interval:
dvds = pi*L.*[  2*cos(psi).*T./x - c2.*tns - fexts; ...
                2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                (tns + 2*zetacnem*profile.*cos(psi).*q./x);
                c2;...
                cos(psi);...
                sin(psi);...
                (3/4)*x.*x.*sin(psi);...
                (1/2)*x;...
                (1/pi)*(x./x0(s0)).*(onevec./(onevec+u));...
                dsq;...
                (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                x.*fextz];

%% specify force balance at SP:
indices = find(~s);
zervec = zeros(size(s(indices)));
onevec = ones(size(s(indices)));
u0 = (tss(indices) - zeta*profile(indices))/(2*K);
C = 0.5*(mss - C0*profile);

dvds(:,indices) =   pi*L*[  zervec;
                            0.5*(C(indices).*tss(indices)-P+fextz(indices));
                            zervec;
                            0.5*C(indices);
                            onevec;
                            zervec;
                            zervec;
                            zervec;
                            (1/pi)*sqrt(onevec./(onevec+u0));
                            zervec;
                            dsw_0*onevec;
                            zervec];

%% specify force balance at NP:
indices = find(~(s-1));
zervec = zeros(size(s(indices)));
onevec = ones(size(s(indices)));
uL = (tss(indices) - zeta*profile(indices))/(2*K);
C = 0.5*(mss - C0*profile);

dvds(:,indices) =   pi*L.*[-0.5*C(indices).*tns(indices);
                            0.5*(C(indices).*tss(indices)-P-fextz(indices));
                            tns(indices);
                            0.5*C(indices);
                           -onevec;
                            zervec;
                            zervec;
                            zervec;
                            (1/pi)*sqrt(onevec./(onevec+uL));
                            zervec;
                            dsw_L*onevec;
                            zervec];
end