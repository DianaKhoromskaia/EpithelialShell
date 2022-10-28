function res = bc(yleft,yright, varpar,  U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt) 
% residual function specifying the boundary conditions for the system of
% ODE's set in ode.m

res = [ yleft(3) - 0.; %dsvn(0)=0
        yleft(5) - 0; %tns(L)=0
        yleft(6) - 0.; %v0(0)=0 
        yleft(7) - 0.; %dX0(0)=0
        yleft(8) - 0.; %vs(0)=0
        yleft(9) - 0,; %dsnew(0) = 0
        yright(3) - 0; %dsvn(L)=0
        yright(6) - 0; %dV=0
        yright(7) - 0; %dX0=0
        yright(8) - 0;%vs(0)=0
        yleft(10); %I1(0)=0
        yright(5)]; %tns(L)=0
      
if strcmp(FixedPar,'P')
    res(8) = yright(6)-varpar(1);
end
        
end


