function res = bc(yleft,yright, varpar,  U, dsU, C1, C2, C, C0, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, kb, x0b, FixedPar, t, P0, thalf_P, tsigma) 
%v = (dsvs, vn, dsvn, mss, tns, dV(s), dX(s), vs, dsnew(s), I(s))

res = [ yleft(3) - 0.; %dsvn(0)=0
        yleft(5) - 0; %tns(0)=0
        yleft(6) - 0.; %v0(0)=0 
        yleft(7) - 0.; %dX0(0)=0
        yleft(8) - 0.; %vs(0)=0
        yleft(9) - 0,; %dsnew(0) = 0
        yright(4) - 0; %mss(L)=0
        yright(6) - 0; %dV(L)=0
        yright(7) - 0;%(-1/xi)*(varpar(2)*xintegral); %dX0(L)=0
        2*K*U(L) + zeta(L) - zetanem(L) + (eta+etab)*yright(1) + (etab-eta)*cos(Psi(L)).*yright(8)./X(L) + (etab*C(L)+eta*(C2(L)-C1(L))).*yright(2)+2*pi*kb*(1-x0b/X(L))*cos(Psi(L)); %tss(L)=2pi*k*(1-x0/x) cos(psi)
        yleft(10); %I1(0)=0
        yright(5)+2*pi*kb*(1-x0b/X(L))*sin(Psi(L))]; %tns(L)=2pi*k*(1-x0/x) sin(psi)

if strcmp(FixedPar,'P')
    res(8) = yright(6)-varpar(1);
end
        
end


