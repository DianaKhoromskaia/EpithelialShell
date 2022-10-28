function [dfdv,dfdvar] = fjac_full_param(s, v, varpar, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar)
% Jacobian for RHS of ode for homogeneous profile (la=L0), parametric mode

switch FixedPar
    case 'V'
        P = varpar(1);
        V = fixedpar;
    case 'P'
        P = fixedpar;
        V = varpar(1);
end
fc = varpar(2);
dsw_0 = varpar(3);
dsw_L = varpar(4);
L = varpar(5);

dfdv = zeros(12,12);
dfdvar = zeros(12,5);

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
        zetanem = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'Bending'
    case 'Force'
    case 'Tension'
        zeta = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'BendingNematic'
        zetacnem = guess(1)-dot([varpar(1) L]-guess(2:end),tangentvec(2:end))/tangentvec(1);
end

T = zetanem*q;
t = tss+T;
u = (t - zeta*onevec)/(2*K);

fexts = fc*sin(psi) + fconst*onevec;
fextn = -cos(psi)*fc;
fextz = sin(psi)*fconst;

switch s
    case 0. %% on south pole
        %% partial derivatives wrt v-components
        u0 = (tss - zeta*onevec)/(2*K);
        
        dfdv(2,1) = pi*L.*(0.5*C);
        dfdv(2,3) = pi*L.*(0.5*tss);
        
        dfdv(4,3) = 0.5*pi*L;
        
        dfdv(9,1) = L.*(-onevec./(4*K*(onevec+u0).^(3/2)));
        
        %% partial derivatives wrt parameters
        % wrt L:
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
            dfdvar(2,1) = -pi*0.5*L;
        end
        
        % wrt fc:
        dfdvar(2,2)  = pi*0.5*L;
        
        % wrt dsw_0:
        dfdvar(11,3) = pi*L;
    case 1 %% on north pole
        %% partial derivatives wrt v-components
        uL = (tss - zeta*onevec)/(2*K);
        
        dfdv(1,2) = -pi*L.*0.5*C;
        dfdv(1,3) = -pi*L.*0.5*tns;
        
        dfdv(2,1) = pi*L.*0.5*C;
        dfdv(2,3) = pi*L.*0.5*tss;
        
        dfdv(3,2) = 0.5*pi*L;
        
        dfdv(4,3) = 0.5*pi*L;
        
        dfdv(9,1) = L.*(-onevec./(4*K*(onevec+uL).^(3/2)));
        
        %% partial derivatives wrt parameters
        % wrt L:
        dfdvar(:,5) =  pi*[-0.5*C.*tns-fconst*onevec;
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
            dfdvar(2,1) = -pi*0.5*L;
        end
        
        % wrt fc:
        dfdvar(2,2) = -pi*0.5*L;
        
        % wrt dsw_L:
        dfdvar(11,4) = pi*L;
        
    otherwise
               %% partial derivatives wrt v-components
                dfdv(1,2)   =  pi*L*(-C+sin(psi)./x);
                dfdv(1,3)   = -pi*L*tns;
                dfdv(1,4)   =  pi*L*((tns.*cos(psi)-2*zetanem*sin(psi).*q)./x - fc*cos(psi));
                dfdv(1,5)   = -pi*L*((tns.*sin(psi)+2*zetanem*cos(psi).*q)./(x.*x));
                dfdv(1,10)  =  pi*L*(2*zetanem*cos(psi)./x);                
                
                dfdv(2,1)   =  pi*L*(C-sin(psi)./x);
                dfdv(2,3)   =  pi*L*tss;
                dfdv(2,4)   =  pi*L*(cos(psi).*(-tss+2*zetanem*q)./x - sin(psi)*fc);
                dfdv(2,5)   =  pi*L*(2*I1+x.*sin(psi).*(tss-2*zetanem*q))./(x.*x.*x);
                dfdv(2,10)  =  pi*L*(2*zetanem*sin(psi)./x);
                dfdv(2,12)  = -pi*L/(x.*x);
                
                dfdv(3,2)   =  0.5*pi*L;
                dfdv(3,4)   = -pi*L*zetacnem*sin(psi).*q./x;
                dfdv(3,5)   = -pi*L*zetacnem*cos(psi).*q./(x.*x);
                dfdv(3,10)  =  pi*L*zetacnem*cos(psi)./x;
                dfdv(3,11)  =  0.5*pi*L*zetacnem;
                
                dfdv(4,3)   =  pi*L;
                dfdv(4,4)   = -pi*L*cos(psi)./x;
                dfdv(4,5)   =  pi*L*sin(psi)./(x.*x);
                
                dfdv(5,4)   = -pi*L*sin(psi);
                
                dfdv(6,4)   =  pi*L*cos(psi);
                
                dfdv(7,4)   =  pi*L*(3/4)*x.*x.*cos(psi);
                dfdv(7,5)   =  pi*L*(3/2)*x.*sin(psi);
                
                dfdv(8,5)   =  0.5*pi*L;
                
                dfdv(9,1)   = -L*2*K*x./((2*K+tss-zeta+zetanem*q).^2 *x0(s0));
                dfdv(9,5)   =  L*2*K./((2*K+tss-zeta+zetanem*q)*x0(s0));
                dfdv(9,9)   = -L*2*K*x.*ds0x0(s0)./((2*K+tss-zeta+zetanem*q)*x0(s0)^2);
                dfdv(9,10)  = -L*2*K*zetanem*x./((2*K+tss-zeta+zetanem*q)^2 *x0(s0));
                
                dfdv(10,11) =  pi*L;
                
                dfdv(11,4)  =  pi*L*sin(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x);
                dfdv(11,5)  =  pi*L*cos(psi).*(dsq.*x - 8*q.*cos(psi))./(x.*x.*x);
                dfdv(11,10) =  pi*L*((-1+3*q.*q)./(2*lc^2)+4*cos(psi).*cos(psi)./(x.*x));
                dfdv(11,11) = -pi*L.*cos(psi)./x;
                
                dfdv(12,4)  =  pi*L*x.*cos(psi)*fconst;
                dfdv(12,5)  =  pi*L*(fc+fconst*sin(psi));
                
                %% partial derivatives parameters

                % wrt L:
                dfdvar(:,5) = pi*[  2*cos(psi).*T./x - c2.*tns - fexts; ...
                                    2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                                    0.5*(tns + zetacnem*(dsq+2*cos(psi).*q./x));...
                                    c2;...
                                    cos(psi);...
                                    sin(psi);...
                                    (3/4)*x.*x.*sin(psi);...
                                    (1/2)*x;...
                                    (1/pi)*(x./x0(v(9,:))).*abs(onevec./(onevec+u));...
                                    dsq;...
                                    (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                                    x.*(fc*onevec + fextz)];
                
                % wrt to P (if applicable)
                if strcmp(FixedPar,'V')
                    dfdvar(2,1) = -pi*0.5*L;
                end
                
                % wrt fc:
                dfdvar(1,2)  = -pi*L*sin(psi);
                dfdvar(2,2)  =  pi*L*cos(psi);
                dfdvar(12,2) =  pi*L*x;
                
        %% new contributions to dF/dpar due to parametric version:
        
        switch ControlPar
            case 'Nematic'
                zetanemred_noL = guess(1)-dot([varpar(1) 0.]-guess(2:end),tangentvec(2:end))/tangentvec(1);
                % wrt P or V:
                dfdvar(1,1) = dfdvar(1,1) + pi*L.*q.*(-tangentvec(2)/tangentvec(1)).*(2*cos(psi)./x);
                dfdvar(2,1) = dfdvar(2,1) + pi*L.*q.*(-tangentvec(2)/tangentvec(1)).*(2*sin(psi)./x);
                dfdvar(9,1) = dfdvar(9,1) + (x./x0(s0))*(-L*(q.*(-tangentvec(2)/tangentvec(1)))/(2*K*(1+u)^2));
                % wrt L: (this is instead of the above)
                dfdvar(1,5) = pi*(2*cos(psi).*(zetanemred_noL*q)./x - c2.*tns - fexts) -4*pi*L.*q.*(tangentvec(3)/tangentvec(1)).*cos(psi)./x;
                dfdvar(2,5) = pi*(2*c1.*(zetanemred_noL*q) + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn) -4*pi*L.*q.*(tangentvec(3)/tangentvec(1)).*sin(psi)./x;
                dfdvar(9,5) = (x./x0(s0))*(1/(1+u) - L*(q.*(-tangentvec(3)/tangentvec(1)))/(2*K*(1+u)^2));
            case 'Bending'
            case 'Force'
            case 'Tension'
            case 'BendingNematic'
                zetacnemred_noL = guess(1)-dot([varpar(1) 0.]-guess(2:end),tangentvec(2:end))/tangentvec(1);
                % wrt P or V:
                dfdvar(3,1) = dfdvar(3,1) + pi*L*0.5*(dsq+2*cos(psi).*q./x).*(-tangentvec(2)/tangentvec(1));
                % wrt L:
                dfdvar(3,5) = pi*0.5*(tns + zetacnemred_noL*(dsq+2*cos(psi).*q./x)) + pi*L*(dsq+2*cos(psi).*q./x).*(-tangentvec(3)/tangentvec(1)); %this is instead of the above
        end
        
        
end


end