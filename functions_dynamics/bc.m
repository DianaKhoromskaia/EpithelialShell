function res = bc(yleft,yright, varpar,  U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, t, P0, thalf_P, tsigma) 
%v = (tss, vn, dsvn, mss, tns, dV(s), dX(s), vs, dsnew(s))
%if mod(i,2)==0
res = [ yleft(3) - 0.; %dsvn(0)=0
        yleft(5) - 0; %tns(0)=0
        yleft(6) - 0.; %v0(0)=0 
        yleft(7) - 0.; %dX0(0)=0
        yleft(8) - 0.; %vs(0)=0
        yleft(9) - 0,; %dsnew(0) = 0
        yright(3) - 0; %dsvn(L)=0
        yright(6) - 0; %dV(L)=0
        yright(7) - 0;%(-1/xi)*(varpar(2)*xintegral); %dX0(L)=0
        yright(8) - 0;%vs(L)=0
        yleft(10); %I1(0)=0
        yright(5)]; %tns(L)=0 instead of I1(L)=0 
        %yright(10)];
% elseif mod(i,2)==1
%     res = [ yleft(3) - 0.; %dsvn(0)=0
%         yleft(5) - 0.; %tns(L)=0
%         yleft(6) - 0.; %v0(0)=0 
%         yleft(7) - 0.; %dX0(0)=0
%         yleft(8) - 0.; %vs(0)=0
%         yleft(9) - 0,; %dsnew(0) = 0
%         yright(3) - 0; %dsvn(L)=0
%         yright(6) - 0.; %dV=0
%         yright(7) - 0.; %dX0=0
%         yright(8) - 0;%vs(0)=0
%         yleft(10); %I1(0)=0
%         yright(10)]; %tns(L)=0 instead of I1(L)=0
% end    
    
if strcmp(FixedPar,'P')
    res(8) = yright(6)-varpar(1);
end
        
end


