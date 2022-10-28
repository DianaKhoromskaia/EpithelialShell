function [dfdv,dfdvar] = fjac_param(s, v, k, varpar, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar)
% Jacobian of RHS of ode for step profile, in parametric mode

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
L1 = varpar(5);
L2 = varpar(6);

dfdv = zeros(12,12);
dfdvar = zeros(12,6);

tss = v(1,:);
tns = v(2,:);
C = v(3,:);
psi = v(4,:);
x = v(5,:);
z = v(6,:);
vol = v(7,:);
area = v(8,:);
s0 = v(9,:);
q = v(10,:);
dsq = v(11,:);
I1 = v(12,:);

c1 = sin(psi)./x;
c2 = C-c1;

onevec = ones(size(s));
zervec = zeros(size(s));

switch ControlPar
    case 'Nematic'
        zetanemred = guess(1)-dot([varpar(1) varpar(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'Bending'
    case 'Force'
    case 'Tension'
        zetared = guess(1)-dot([varpar(1) varpar(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'BendingNematic'
        zetacnemred = guess(1)-dot([varpar(1) varpar(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
end

switch k
    case 1 %% red region
        
        zeta = zeta + zetared;
        zetanem = zetanem + zetanemred;
        fconst = fconst + fconstred;
        zetacnem = zetacnem + zetacnemred;
        
        T = zetanem*q;
        t = tss+T;
        u = (t - zeta*onevec)/(2*K);
        
        fexts = fc*sin(psi) + fconst*onevec;
        fextn = -cos(psi)*fc;
        fextz = sin(psi)*fconst;
        
        switch s
            case 0.
                %% partial derivatives wrt v-components
                u0 = (tss - zeta*onevec)/(2*K);
                
                dfdv(2,1) = pi*L1.*(0.5*C);
                dfdv(2,3) = pi*L1.*(0.5*tss);
                
                dfdv(4,3) = 0.5*pi*L1;
                
                dfdv(9,1) = L1.*(-onevec./(4*K*(onevec+u0).^(3/2)));
                
                %% partial derivatives wrt parameters
                % wrt L1:
                dfdvar(:,5) =  pi*[-fconst*onevec;
                                    0.5*(C.*tss-P+fc);
                                    zervec;
                                    0.5*C;
                                    onevec;
                                    zervec;
                                    zervec;
                                    zervec;
                                    (1/pi)*sqrt(onevec./abs(onevec+u0));
                                    zervec;
                                    dsw_0;
                                    zervec];
                
                % wrt P:
                if strcmp(FixedPar,'V')
                    dfdvar(2,1) = -pi*0.5*L1;
                end
                
                % wrt fc:
                dfdvar(2,2)  = pi*0.5*L1;
                
                % wrt dsw_0:
                dfdvar(11,3) = pi*L1;
                
            otherwise              
               %% partial derivatives wrt v-components
                dfdv(1,2)   =  pi*L1*(-C+sin(psi)./x);
                dfdv(1,3)   = -pi*L1*tns;
                dfdv(1,4)   =  pi*L1*((tns.*cos(psi)-2*zetanem*sin(psi).*q)./x - fc*cos(psi));
                dfdv(1,5)   = -pi*L1*((tns.*sin(psi)+2*zetanem*cos(psi).*q)./(x.*x));
                dfdv(1,10)  =  pi*L1*(2*zetanem*cos(psi)./x);                
                
                dfdv(2,1)   =  pi*L1*(C-sin(psi)./x);
                dfdv(2,3)   =  pi*L1*tss;
                dfdv(2,4)   =  pi*L1*(cos(psi).*(-tss+2*zetanem*q)./x - sin(psi)*fc);
                dfdv(2,5)   =  pi*L1*(2*I1+x.*sin(psi).*(tss-2*zetanem*q))./(x.*x.*x);
                dfdv(2,10)  =  pi*L1*(2*zetanem*sin(psi)./x);
                dfdv(2,12)  = -pi*L1/(x.*x);
                
                dfdv(3,2)   =  0.5*pi*L1;
                dfdv(3,4)   = -pi*L1*zetacnem*sin(psi).*q./x;
                dfdv(3,5)   = -pi*L1*zetacnem*cos(psi).*q./(x.*x);
                dfdv(3,10)  =  pi*L1*zetacnem*cos(psi)./x;
                dfdv(3,11)  =  0.5*pi*L1*zetacnem;
                
                dfdv(4,3)   =  pi*L1;
                dfdv(4,4)   = -pi*L1*cos(psi)./x;
                dfdv(4,5)   =  pi*L1*sin(psi)./(x.*x);
                
                dfdv(5,4)   = -pi*L1*sin(psi);
                
                dfdv(6,4)   =  pi*L1*cos(psi);
                
                dfdv(7,4)   =  pi*L1*(3/4)*x.*x.*cos(psi);
                dfdv(7,5)   =  pi*L1*(3/2)*x.*sin(psi);
                
                dfdv(8,5)   =  0.5*pi*L1;
                
                dfdv(9,1)   = -L1*2*K*x./((2*K+tss-zeta+zetanem*q).^2 *x0(s0));
                dfdv(9,5)   =  L1*2*K./((2*K+tss-zeta+zetanem*q)*x0(s0));
                dfdv(9,9)   = -L1*2*K*x.*ds0x0(s0)./((2*K+tss-zeta+zetanem*q)*x0(s0)^2);
                dfdv(9,10)  = -L1*2*K*zetanem*x./((2*K+tss-zeta+zetanem*q)^2 *x0(s0));
                
                dfdv(10,11) =  pi*L1;
                
                dfdv(11,4)  =  pi*L1*sin(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x);
                dfdv(11,5)  =  pi*L1*cos(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x.*x);
                dfdv(11,10) =  pi*L1*((-1+3*q.*q)./(2*lc^2)+4*cos(psi).*cos(psi)./(x.*x));
                dfdv(11,11) = -pi*L1.*cos(psi)./x;
                
                dfdv(12,4)  =  pi*L1*x.*cos(psi)*fconst;
                dfdv(12,5)  =  pi*L1*(fc+fconst*sin(psi));
                
                %% partial derivatives parameters

                % wrt L1:
                dfdvar(:,5) = pi*[  2*cos(psi).*T./x - c2.*tns - fexts; ...
                                    2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                                    0.5*(tns + zetacnem*(dsq+2*cos(psi).*q./x));...
                                    c2;...
                                    cos(psi);...
                                    sin(psi);...
                                    (3/4)*x.*x.*sin(psi);...
                                    (1/2)*x;...
                                    (1/pi)*(x./x0(v(9,:))).*abs(onevec./(onevec+u));...
                                    I1;...
                                    (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                                    x.*(fc*onevec + fextz)];
                
                % wrt to P (if applicable)
                if strcmp(FixedPar,'V')
                    dfdvar(2,1) = -pi*0.5*L1;
                end
                
                % wrt fc:
                dfdvar(1,2)  = -pi*L1*sin(psi);
                dfdvar(2,2)  =  pi*L1*cos(psi);
                dfdvar(12,2) =  pi*L1*x;
               
                %% new contributions to dF/dpar due to parametric version:

                switch ControlPar
                    case 'Nematic'
                        zetanemred_noL1 = guess(1)-dot([varpar(1) 0.]-guess(2:end),tangentvec(2:end))/tangentvec(1);
                        % wrt P or V:
                        dfdvar(1,1) = dfdvar(1,1) + pi*L1.*q.*(-tangentvec(2)/tangentvec(1)).*(2*cos(psi)./x);
                        dfdvar(2,1) = dfdvar(2,1) + pi*L1.*q.*(-tangentvec(2)/tangentvec(1)).*(2*sin(psi)./x);
                        dfdvar(9,1) = dfdvar(9,1) + (x./x0(s0))*(-L1*(q.*(-tangentvec(2)/tangentvec(1)))/(2*K*(1+u)^2));
                        % wrt L1: (this is instead of the above)
                        dfdvar(1,5) = pi*(2*cos(psi).*(zetanemred_noL1*q)./x - c2.*tns - fexts) -4*pi*L1.*q.*(tangentvec(3)/tangentvec(1)).*cos(psi)./x;
                        dfdvar(2,5) = pi*(2*c1.*(zetanemred_noL1*q) + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn) -4*pi*L1.*q.*(tangentvec(3)/tangentvec(1)).*sin(psi)./x;
                        dfdvar(9,5) = (x./x0(s0))*(1/(1+u) - L1*(q.*(-tangentvec(3)/tangentvec(1)))/(2*K*(1+u)^2));
                    case 'Bending'
                    case 'Force'
                    case 'Tension'
                    case 'BendingNematic'
                        zetacnemred_noL1 = guess(1)-dot([varpar(1) 0.]-guess(2:end),tangentvec(2:end))/tangentvec(1);
                        % wrt P or V:
                        dfdvar(3,1) = dfdvar(3,1) + pi*L1*0.5*(dsq+2*cos(psi).*q./x).*(-tangentvec(2)/tangentvec(1));
                        % wrt L1:
                        dfdvar(3,5) = pi*0.5*(tns + zetacnemred_noL1*(dsq+2*cos(psi).*q./x)) + pi*L1*(dsq+2*cos(psi).*q./x).*(-tangentvec(3)/tangentvec(1)); %this is instead of the above
                end

        end
        
    case 2 %% blue region
                
        T = zetanem*q;
        t = tss+T;
        u = (t - zeta*onevec)/(2*K);
        
        fexts = fc*sin(psi) + fconst*onevec;
        fextn = -cos(psi)*fc;
        fextz = sin(psi)*fconst;
        
        switch s
            case 2
                %% partial derivatives wrt v-components
                uL = (tss - zeta*onevec)/(2*K);
                
                dfdv(1,2) = -pi*L2.*0.5*C;
                dfdv(1,3) = -pi*L2.*0.5*tns;
                
                dfdv(2,1) = pi*L2.*0.5*C;
                dfdv(2,3) = pi*L2.*0.5*tss;
                
                dfdv(3,2) = 0.5*pi*L2;
                
                dfdv(4,3) = 0.5*pi*L2;
                
                dfdv(9,1) = L2.*(-onevec./(4*K*(onevec+uL).^(3/2)));
                
                %% partial derivatives wrt parameters
                % wrt L2:
                dfdvar(:,6) =  pi*[- 0.5*C.*tns-fconst*onevec;
                    0.5*(C.*tss-P-fc);
                    0.5*tns;
                    0.5*C;
                    -onevec;
                    zervec;
                    zervec;
                    zervec;
                    (1/pi)*sqrt(onevec./abs(onevec+uL));
                    zervec;
                    dsw_L;
                    zervec];
                
                % wrt P (if applicable):
                if strcmp(FixedPar,'V')
                    dfdvar(2,1) = -pi*0.5*L2;
                end
                
                % wrt fc:
                dfdvar(2,2) = -pi*0.5*L2;
                
                % wrt dsw_L:
                dfdvar(11,4) = pi*L2;
                
            otherwise
                %% partial derivatives wrt v-components
                dfdv(1,2)   =  pi*L2*(-C+sin(psi)./x);
                dfdv(1,3)   = -pi*L2*tns;
                dfdv(1,4)   =  pi*L2*((tns.*cos(psi)-2*zetanem*sin(psi).*q)./x - fc*cos(psi));
                dfdv(1,5)   = -pi*L2*((tns.*sin(psi)+2*zetanem*cos(psi).*q)./(x.*x));
                dfdv(1,10)  =  pi*L2*(2*zetanem*cos(psi)./x);                
                
                dfdv(2,1)   =  pi*L2*(C-sin(psi)./x);
                dfdv(2,3)   =  pi*L2*tss;
                dfdv(2,4)   =  pi*L2*(cos(psi).*(-tss+2*zetanem*q)./x - sin(psi)*fc);
                dfdv(2,5)   =  pi*L2*(2*I1+x.*sin(psi).*(tss-2*zetanem*q))./(x.*x.*x);
                dfdv(2,10)  =  pi*L2*(2*zetanem*sin(psi)./x);
                dfdv(2,12)  = -pi*L2/(x.*x);
                
                dfdv(3,2)   =  0.5*pi*L2;
                dfdv(3,4)   = -pi*L2*zetacnem*sin(psi).*q./x;
                dfdv(3,5)   = -pi*L2*zetacnem*cos(psi).*q./(x.*x);
                dfdv(3,10)  =  pi*L2*zetacnem*cos(psi)./x;
                dfdv(3,11)  =  0.5*pi*L2*zetacnem;
                
                dfdv(4,3)   =  pi*L2;
                dfdv(4,4)   = -pi*L2*cos(psi)./x;
                dfdv(4,5)   =  pi*L2*sin(psi)./(x.*x);
                
                dfdv(5,4)   = -pi*L2*sin(psi);
                
                dfdv(6,4)   =  pi*L2*cos(psi);
                
                dfdv(7,4)   =  pi*L2*(3/4)*x.*x.*cos(psi);
                dfdv(7,5)   =  pi*L2*(3/2)*x.*sin(psi);
                
                dfdv(8,5)   =  0.5*pi*L2;
                
                dfdv(9,1)   = -L2*2*K*x./((2*K+tss-zeta+zetanem*q).^2 *x0(s0));
                dfdv(9,5)   =  L2*2*K./((2*K+tss-zeta+zetanem*q)*x0(s0));
                dfdv(9,9)   = -L2*2*K*x.*ds0x0(s0)./((2*K+tss-zeta+zetanem*q)*x0(s0)^2);
                dfdv(9,10)  = -L2*2*K*zetanem*x./((2*K+tss-zeta+zetanem*q)^2 *x0(s0));
                
                dfdv(10,11) =  pi*L2;
                
                dfdv(11,4)  =  pi*L2*sin(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x);
                dfdv(11,5)  =  pi*L2*cos(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x.*x);
                dfdv(11,10) =  pi*L2*((-1+3*q.*q)./(2*lc^2)+4*cos(psi).*cos(psi)./(x.*x));
                dfdv(11,11) = -pi*L2.*cos(psi)./x;
                
                dfdv(12,4)  =  pi*L2*x.*cos(psi)*fconst;
                dfdv(12,5)  =  pi*L2*(fc+fconst*sin(psi));
                
                %% partial derivatives wrt parameters
                % wrt to L2:                             
                dfdvar(:,6) = pi*[  2*cos(psi).*T./x - c2.*tns - fexts; ...
                                    2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                                    0.5*(tns + zetacnem*(dsq+2*cos(psi).*q./x));...
                                    c2;...
                                    cos(psi);...
                                    sin(psi);...
                                    (3/4)*x.*x.*sin(psi);...
                                    (1/2)*x;...
                                    (1/pi)*(x./x0(v(9,:))).*abs(onevec./(onevec+u));...
                                    I1;...
                                    (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                                    x.*(fc*onevec + fextz)];
                % wrt P (if applicable):
                if strcmp(FixedPar,'V')   
                    dfdvar(2,1) = -pi*0.5*L2;
                end
                
                % wrt to fc:
                dfdvar(1,2)  = -pi*L2*sin(psi);
                dfdvar(2,2)  =  pi*L2*cos(psi);
                dfdvar(12,2) =  pi*L2*x;
        end
        
end
end