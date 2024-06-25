function res = bc(yleft,yright, varpar,  U, dsU, C1, C2, C, C0, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, dskappa, K, xi, kb, x0b, kp, psi0b, FixedPar, t, P0, thalf_P, tsigma) 
%v = (dsvs, vn, dsvn, mss, tns, dV(s), dX(s), vs, dsnew(s), I(s))

res = [ yleft(3) - 0.; %dsvn(0)=0
        yleft(5) - 0; %tns(0)=0
        yleft(6) - 0.; %v0(0)=0 
        yleft(7) - 0.; %dX0(0)=0
        yleft(8) - 0.; %vs(0)=0
        yleft(9) - 0,; %dsnew(0) = 0
        yright(4) + kp*(Psi(L)-psi0b); %mss(L) open boundary
        yright(6) - 0; %dV(L)=0
        yright(7) - 0;%(-1/xi)*(varpar(2)*xintegral); %dX0(L)=0
        2*K*U(L) + zeta(L) - zetanem(L) - (2*kappa(L)*(C(L)-C0)+zetac(L)).*(C2(L)-0.5*C(L)) - (kappa(L)*C0-0.5*zetac(L))*C(L) + (eta+etab)*yright(1) + (etab-eta)*cos(Psi(L)).*yright(8)./X(L) + (etab*C(L)+eta*(C2(L)-C1(L))).*yright(2)+2*pi*kb*(1-x0b/X(L))*cos(Psi(L))+0.5*kp*(Psi(L)-psi0b)^2*cos(Psi(L))/X(L); %tss open boundary
        yleft(10); %I1(0)=0
        yright(5)+2*pi*kb*(1-x0b/X(L))*sin(Psi(L))+0.5*kp*(Psi(L)-psi0b)^2*sin(Psi(L))/X(L)]; %tns(L) open boundary


if strcmp(FixedPar,'P')
    res(8) = yright(6)-varpar(1);
end
        
end


