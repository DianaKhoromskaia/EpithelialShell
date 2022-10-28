function y = yguessfun_init(s, k, la, zeta) %k - region
% intial guess for functions t, tns, C, psi, x, z, v, a, s0, L, q, dsq, I1
% in form of the solution on the sphere; all dimensionless
% s0 in [0,1],[1,2] and s in [0,1],[1,2]

switch k
    case 1
        s = s*la; %retrieve s0 coordinate 
        y = [zeta
            0.
            2
            s*pi
            sin(s*pi)
            1-cos(s*pi)
            (2+cos(s*pi)).*(sin(s*pi/2).^4)
            (1-cos(s*pi))
            s
            1;
            0;
            0]; 
    case 2
        s = la + (1-la)*(s-1); %retrieve s0 coordinate 
        y = [zeta
            0.
            2
            s*pi
            sin(s*pi)
            1-cos(s*pi)
            (2+cos(s*pi)).*(sin(s*pi/2).^4)
            (1-cos(s*pi))
            s
            1;
            0;
            0];
end

end