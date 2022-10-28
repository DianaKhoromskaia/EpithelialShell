function y = yguessfun_init_fullint(s, L0, zeta) %k - region
% intial guess for functions t, tns, C, psi, x, z, v, a, s0, L, q, dsq, I1
% in form of the solution on the sphere; all dimensionless
% for the case of homogeneous profile (la=L0)
% s0 in [0,1],[1,2] and s in [0,1],[1,2]
s = s*L0; %retrieve s0 coordinate

Qinit = dlmread('Qinit.dat','\t');
q = griddedInterpolant(Qinit(2,:),Qinit(3,:));
dsq = griddedInterpolant(Qinit(2,:),Qinit(4,:));

y = [zeta
    0.
    2*2 %mss=2*Ckk
    s*pi
    sin(s*pi)
    1-cos(s*pi)
    (2+cos(s*pi)).*(sin(s*pi/2).^4)
    (1-cos(s*pi))
    s
    q(s);
    dsq(s);
    0];

end