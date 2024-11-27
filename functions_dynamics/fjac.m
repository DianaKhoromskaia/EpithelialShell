function [dfdv,dfdvar] = fjac(s, v, varpar, U, dsU, C1, C2, C, C0, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, dskappa, K, xi, kb, x0b, kp, psi0b, FixedPar, t, P0, thalf_P, tsigma, zc0_on)
%Jacobian of RHS for the bvp4c solver
%v = (dsvs, vn, dsvn, mss, tns, dV(s), dX(s), vs, dsnew(s), I(s))

dfdv = zeros(10,10);
dfdvar = zeros(10,2);

dsvs = v(1,:);
vn = v(2,:);
dsvn = v(3,:);
vs = v(8,:);

switch s
    case 0.
        %partial derivatives wrt v-components 
        
        dfdv(3,2) = -0.5*(C1(0).^2 + C2(0).^2);
        dfdv(3,4) = -0.5/etacb;
        
        dfdv(4,5) = 1;
        
        dfdv(5,1) = C2(0)*(2*etab);
        dfdv(5,2) = 0.5*xi + C2(0)*(etab*C(0));
        
        dfdv(8,1) = 1;
         
        dfdv(9,1) = 1;
        dfdv(9,2) = C2(0);
        
        %partial derivatives wrt to parameters
        if strcmp(FixedPar,'V')
            dfdvar(5,1) = -0.5;
        elseif strcmp(FixedPar,'P')
            dfdvar(5,1) = -0.5*(-etap);
        end
        dfdvar(5,2) = 0.5;

    otherwise
        
        %partial derivatives wrt v-components
        dfdv(1,1) = -cos(Psi(s))./X(s);
        dfdv(1,2) = -dsC(s);
        dfdv(1,3) = -(etab*C(s) + eta*(C2(s)-C1(s)))/(eta+etab);
        dfdv(1,5) = - C2(s)/(eta+etab);
        dfdv(1,8) = (cos(Psi(s))./X(s)).^2 - ((eta-etab)/(eta+etab))*C1(s).*C2(s);

        dfdv(2,3) = 1;
        
        dfdv(3,2) = -(C1(s).^2 + C2(s).^2);
        dfdv(3,3) = -cos(Psi(s))./X(s);
        dfdv(3,4) = -1/etacb;
        dfdv(3,8) = dsC(s);
        
        dfdv(4,5) = 1;
        
        dfdv(5,1) = (C1(s)+C2(s))*(eta+etab)+2*C1(s)*(-eta);
        dfdv(5,2) = xi+(C1(s)+C2(s))*((etab*C(s)+eta*(C2(s)-C1(s))))+2*C1(s)*eta*((C1(s)-C2(s)));
        dfdv(5,5) = -cos(Psi(s))./X(s);
        dfdv(5,8) = (C1(s)+C2(s))*((etab-eta)*cos(Psi(s))./X(s))+2*C1(s)*eta*cos(Psi(s))./X(s);
  
        dfdv(6,2) = 2*pi*X(s);
        
        dfdv(7,2) = X(s).*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;

        dfdv(8,1) = 1;

        dfdv(9,1) = 1;
        dfdv(9,2) = C2(s);
        
        dfdv(10,2) = xi*X(s).*cos(Psi(s));
        
        %partial derivatives wrt to parameters
        if strcmp(FixedPar,'V')
            dfdvar(5,1) = -1;%-1;%-0.5;
        elseif strcmp(FixedPar,'P')
            dfdvar(5,1) = -1*(-etap);%-1*(-etap);%-0.5*(-etap);
        end

        dfdvar(1,2) = - sin(Psi(s))/(eta+etab);
        dfdvar(5,2) = cos(Psi(s));
        dfdvar(10,2) = X(s);
end


end

