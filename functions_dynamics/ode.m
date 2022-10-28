function dvds = ode(s, v, varpar, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt)
% RHS of the system of ODE's representing the force balance on the surface at time t together with constraints, 
% written for the vector 
% v = (\partial_s v_s, v_n, \partial_s v_n, \bar{m}_s^s, t^n_s, r_v,
% r_c/xintegral, v_s, \partial_s s_{new}, I)
% and solved on the interval [0,L]


if strcmp(FixedPar,'V')
    P = varpar(1);
elseif strcmp(FixedPar,'P')
    P = -etap*varpar(1); % varpar(1)=dV/dt
end

fc=varpar(2);

dsvs = v(1,:);
vn = v(2,:);
dsvn = v(3,:);
mss = v(4,:);
tns = v(5,:);
%v0 = v(6,:);
%x0 = v(7,:);
vs = v(8,:);
dsnew = v(9,:);
I1 = v(10,:);

T = zetanem(s);
V = eta*(cos(Psi(s)).*vs./X(s)-dsvs+(C1(s)-C2(s)).*vn);

fextn = -xi*vn;

onevec = ones(size(s));
zervec = zeros(size(s));

ds2vn = -cos(Psi(s)).*dsvn./X(s) - vn.*(C1(s).^2 + C2(s).^2) + vs.*dsC(s) - (mss - 2*kappa*C(s) - zetac(s) + zetacnem(s))/etacb;

tss = 2*K*U(s) + zeta(s) - zetanem(s) + (eta+etab)*dsvs + (etab-eta)*cos(Psi(s)).*vs./X(s) + (etab*C(s)+eta*(C2(s)-C1(s))).*vn;
dstns = 2*C1(s).*(T+V) + (C1(s)+C2(s)).*tss - cos(Psi(s)).*tns./X(s) - P + cos(Psi(s)).*fc  - fextn;
dsmss = tns + 2*cos(Psi(s)).*zetacnem(s)./X(s);

ds2vs = -cos(Psi(s)).*(dsvs-cos(Psi(s)).*vs./X(s))./X(s) - (eta-etab)*C1(s).*C2(s).*vs/(eta+etab) - dsC(s).*vn - (etab*C(s)+eta*(C2(s)-C1(s))).*dsvn/(eta+etab) - (2*K*dsU(s)+dszeta(s))/(eta+etab)  + (dszetanem(s) + 2*cos(Psi(s)).*zetanem(s)./X(s))/(eta+etab) - C2(s).*tns./(eta+etab) - sin(Psi(s)).*fc/(eta+etab);

% force balance within the interval:
dvds = [ds2vs;
    dsvn;
    ds2vn;
    dsmss;
    dstns;
    2*pi*X(s).*vn;
    X(s).*vn.*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;
    dsvs;
    dsvs + C2(s).*vn;
    X(s).*(fc.*onevec - cos(Psi(s)).*fextn)];

% force balance at SP:
        indices = find(~s);
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        
        dvds(:,indices) =   [zervec;
            zervec;
            0.5*(-vn(indices).*(C1(0).^2 + C2(0).^2) - (mss(indices)-2*kappa*C(0)- zetac(0))/etacb);
            tns(indices);
            -0.5*P*onevec+C2(0).*(2*K*U(0) + zeta(0) - zetanem(0) + (eta+etab)*dsvs(indices) + etab*C(0)*vn(indices))+0.5*xi*vn(indices)+0.5*fc*onevec;
            zervec;
            zervec;
            dsvs(indices);
            dsvs(indices) + C2(0).*vn(indices);
            zervec];

 % force balance at NP:       
        indices = find(~(s-L));
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        
        dvds(:,indices) =   [zervec;
            zervec;
            0.5*(-vn(indices).*(C1(L).^2 + C2(L).^2) - (mss(indices)-2*kappa*C(L)- zetac(L))/etacb);
            tns(indices);
            -0.5*P*onevec+C2(L).*(2*K*U(L) + zeta(L) - zetanem(L) + (eta+etab)*dsvs(indices) + etab*C(L)*vn(indices))+0.5*xi*vn(indices)-0.5*fc.*onevec;
            zervec;
            zervec;
            dsvs(indices);
            dsvs(indices) + C2(L).*vn(indices);
            zervec]; 
        
end

