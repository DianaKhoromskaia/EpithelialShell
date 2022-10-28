function dvds = ode(s, v, k, varpar, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar)
% RHS of ode for step profile, non-parametric mode
%v = (   tss,  tns,    C,  psi,    x,    z, v(s), a(s),   s0,      q,   dsq,    I1)
%  = (  v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9),  v(10), v(11), v(12))

switch FixedPar
    case 'V'
        P = varpar(1);
    case 'P'
        P = fixedpar;
end
fc = varpar(2);
dsw_0 = varpar(3);
dsw_L = varpar(4);
L1 = varpar(5);
L2 = varpar(6);

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
     
        %% RHS in main interval:
        dvds = pi*L1.*[ 2*cos(psi).*T./x - c2.*tns - fexts; ...
                        2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                        0.5*(tns + zetacnem*(dsq+2*cos(psi).*q./x));...
                        c2;...
                        cos(psi);...
                        sin(psi);...
                        (3/4)*x.*x.*sin(psi);...
                        (1/2)*x;...
                        (1/pi)*(x./x0(s0)).*(onevec./(onevec+u));...
                        dsq;...
                        (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                        x.*(fc*onevec + fextz)];
        
        %% specify force balance at SP:
        indices = find(~s);
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        u0 = (tss(indices) - zeta*onevec)/(2*K);
        
        dvds(:,indices) =   pi*L1.*[-fconst*onevec;
                                    0.5*(C(indices).*tss(indices)-P+fc);
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
        
    case 2 %% blue region
        
        T = zetanem*q;
        t = tss+T;
        u = (t - zeta*onevec)/(2*K);
        
        fexts = fc*sin(psi) + fconst*onevec;
        fextn = -cos(psi)*fc;
        fextz = sin(psi)*fconst;
        
        %% RHS in main interval:
        dvds = pi*L2.*[ 2*cos(psi).*T./x - c2.*tns - fexts; ...
                        2*c1.*T + c2.*tss - (P/2)*onevec - I1./(x.*x) - fextn;...
                        0.5*(tns + zetacnem*(dsq+2*cos(psi).*q./x));...
                        c2;...
                        cos(psi);...
                        sin(psi);...
                        (3/4)*x.*x.*sin(psi);...
                        (1/2)*x;...
                        (1/pi)*(x./x0(s0)).*(onevec./(onevec+u));...
                        dsq;...
                        (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(psi)./x).*(4*(cos(psi)./x).*q - dsq);...
                        x.*(fc*onevec + fextz)];
        
        %% specify force balance at NP:
        indices = find(~(s-2));
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        uL = (tss(indices) - zeta*onevec)/(2*K);
        
        dvds(:,indices) =  pi*L2.*[-0.5*C(indices).*tns(indices)-fconst*onevec;
                                    0.5*(C(indices).*tss(indices)-P-fc);
                                    0.5*tns(indices);
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

end