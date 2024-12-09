function dvds = ode(s, v, varpar, U, dsU, C1, C2, C, C0, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, dskappa, K, xi, FixedPar, t, P0, thalf_P, tsigma, zc0_on)
% rhs for ode for the vector v = (dsvs, vn, dsvn, mss, tns, dV(s), dX(s),
% vs, dsnew(s), I(s)), solved on interval [0,L]
%figure(60)
%plot(s,v)
%legend('dsvs', 'vn', 'dsvn', 'mss', 'tns', 'dV(s)', 'dX(s)','vs', 'dsnew(s)', 'I')

if strcmp(FixedPar,'V')
    P = varpar(1);
elseif strcmp(FixedPar,'P')
    P = P0*(1-sigmoidal(t,thalf_P,tsigma))-etap*varpar(1); %here varpar(1)=dV/dt
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
% dsnew = v(9,:);
% I1 = v(10,:);

T = zetanem(s);
V = eta*(cos(Psi(s)).*vs./X(s)-dsvs+(C1(s)-C2(s)).*vn);

fextn = -xi*vn;

onevec = ones(size(s));

ds2vn = -cos(Psi(s)).*dsvn./X(s) - vn.*(C1(s).^2 + C2(s).^2) + vs.*dsC(s) - (mss - 2*kappa(s).*(C(s)-C0) - zetac(s) + zetacnem(s))/etacb;

tss = 2*K*U(s) + zeta(s) - zetanem(s) -(kappa(s).*(C(s)-C0)+zc0_on*0.5*zetac(s)).*(C2(s)-C1(s)) + (zc0_on*0.5*zetac(s)-C0*kappa(s)).*(C(s)-C0+zc0_on*0.5*zetac(s)./kappa(s)) + (eta+etab)*dsvs + (etab-eta)*cos(Psi(s)).*vs./X(s) + (etab*C(s)+eta*(C2(s)-C1(s))).*vn;
dstns = 2*C1(s).*(T+V+(kappa(s).*(C(s)-C0)+zc0_on*0.5*zetac(s)).*(C2(s)-C1(s))) + C(s).*tss - cos(Psi(s)).*tns./X(s) - P + cos(Psi(s)).*fc  - fextn;
dsmss = tns + 2*cos(Psi(s)).*zetacnem(s)./X(s);

ds2vs = -cos(Psi(s)).*(dsvs-cos(Psi(s)).*vs./X(s))./X(s) - (eta-etab)*C1(s).*C2(s).*vs/(eta+etab) - dsC(s).*vn - (etab*C(s)+eta*(C2(s)-C1(s))).*dsvn/(eta+etab) - (2*K*dsU(s)+dszeta(s))/(eta+etab)  + (dszetanem(s) + 2*cos(Psi(s)).*zetanem(s)./X(s))/(eta+etab) - C2(s).*tns./(eta+etab) - sin(Psi(s)).*fc/(eta+etab) + (dskappa(s).*C(s).*(C2(s)-C1(s)) + 2*C2(s).*kappa(s).*dsC(s)+ zc0_on*C1(s).*dszetac(s) + (zc0_on*0.5*zetac(s)/kappa(s)-C0).*(zc0_on*dszetac(s)+(C0-zc0_on*0.5*zetac(s)/kappa(s)).*dskappa(s)/kappa(s)))./(eta+etab);

% force balance within the interval:
dvds = [ds2vs;
    dsvn;
    ds2vn;
    dsmss;
    dstns;
    2*pi*X(s).*vn;
    X(s).*vn.*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;
    dsvs;
    dsvs + C2(s).*vn;%sqrt(onevec+2*dt*(dsvs + C2(s).*vn)+(dt^2)*((dsvs + C2(s).*vn).^2+(dsvn-C2(s).*vs).^2));%dsvs + C2(s).*vn;
    X(s).*(fc.*onevec - cos(Psi(s)).*fextn)];

% specify force balance at SP:
        indices = find(~s);
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        
        dvds(:,indices) =   [(2/3)*(-dszeta(0) + zc0_on*dszetac(0).*C1(0) + (zc0_on*0.5*zetac(0)./kappa(0)-C0).*(zc0_on*dszetac(0)+(C0-zc0_on*0.5*zetac(0)./kappa(0)).*dskappa(0)./kappa(0)))/(eta+etab);%-C2(0).*tns(indices);%zervec;
            zervec;
            0.5*(-vn(indices).*(C1(0).^2 + C2(0).^2) - (mss(indices) - 2*kappa(0).*(C(0)-C0) - zetac(0))/etacb);
            tns(indices);
            -0.5*P*onevec+C2(indices).*(2*K*U(indices) + zeta(indices) + (zc0_on*0.5*zetac(indices)-kappa(indices)*C0).*(C(indices)-C0+zc0_on*0.5*zetac(indices)./kappa(indices)) + 2*etab*dsvs(indices) + etab*C(indices).*vn(indices))+0.5*xi*vn(indices)+0.5*fc*onevec;%-0.5*P*onevec+C2(0).*tss(indices)+0.5*xi*vn(indices)+0.5*fc*onevec;
            zervec;
            zervec;
            dsvs(indices);%0.5*((tss(indices) - 2*K*U(0) - zeta(0))/etab - C(0).*vn(indices));
            dsvs(indices) + C2(indices).*vn(indices);%sqrt(onevec+2*dt*(dsvs(indices) + C2(0).*vn(indices))+(dt^2)*((dsvs(indices) + C2(0).*vn(indices)).^2+(dsvn(indices)-C2(0).*vs(indices)).^2));%sqrt(onevec+dt*(tss(indices) - 2*K*U(0) - zeta(0))/etab+(dt^2)*(tss(indices) - 2*K*U(0) - zeta(0)).^2/(4*etab^2));%(0.5*((tss(indices) - 2*K*U(0) - zeta(0))/etab - C(0).*vn(indices))) + C2(0)*vn(indices);
            zervec];

 % specify force balance at NP:       
        indices = find(~(s-L));
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        
        dvds(:,indices) =   [(2/3)*(-dszeta(indices) + zc0_on*dszetac(indices).*C1(indices) + (zc0_on*0.5*zetac(indices)./kappa(indices)-C0).*(zc0_on*dszetac(indices)+(C0-zc0_on*0.5*zetac(indices)./kappa(indices)).*dskappa(indices)./kappa(indices)))/(eta+etab);%-C2(0).*tns(indices);%zervec;
            zervec;
            0.5*(-vn(indices).*(C1(L).^2 + C2(L).^2) - (mss(indices)-2*kappa(L).*(C(L)-C0) - zetac(L))/etacb);
            tns(indices);
            -0.5*P*onevec+C2(L).*(2*K*U(L) + zeta(L) + (zc0_on*0.5*zetac(L)-kappa(L)*C0).*(C(L)-C0+zc0_on*0.5*zetac(L)./kappa(L)) + 2*etab*dsvs(indices) + etab*C(L).*vn(indices))+0.5*xi*vn(indices)-0.5*fc*onevec;%-0.5*P*onevec+C2(0).*tss(indices)+0.5*xi*vn(indices)+0.5*fc*onevec;
            zervec;
            zervec;
            dsvs(indices);%0.5*((tss(indices) - 2*K*U(L) - zeta(L))/etab - C(L).*vn(indices));
            dsvs(indices) + C2(L).*vn(indices);%sqrt(onevec+2*dt*(dsvs(indices) + C2(L).*vn(indices))+(dt^2)*((dsvs(indices) + C2(L).*vn(indices)).^2+(dsvn(indices)-C2(L).*vs(indices)).^2));%sqrt(onevec+dt*(tss(indices) - 2*K*U(L) - zeta(L))/etab+(dt^2)*(tss(indices) - 2*K*U(L) - zeta(L)).^2/(4*etab^2));%-onevec+sqrt(onevec+2*dt*((0.5*((tss(indices) - 2*K*U(L) - zeta(L))/etab - C(L).*vn(indices))) + C2(L)*vn(indices))+(dt)^2*((0.5*((tss(indices) - 2*K*U(L) - zeta(L))/etab - C(L).*vn(indices))) + C2(L)*vn(indices)));%(0.5*((tss(indices) - 2*K*U(L) - zeta(L))/etab - C(L).*vn(indices))) + C2(L)*vn(indices);
            zervec]; 
        
end