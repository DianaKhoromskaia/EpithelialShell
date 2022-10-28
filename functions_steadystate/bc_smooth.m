function res = bc_smooth(yleft, yright, varpar, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar, Profile, width)
% residual function for solution with smooth profile, in non-parametric mode

%v = (t,   tns,  C  , psi,    x,    z,     v,     a,  s0,   q,  dsq,  I1,  L)
%  = (1,     2,  3  ,   4,    5,    6,     7,     8,   9,  10,  11,   12, 13)

switch FixedPar
    case 'V'
        V = fixedpar;
    case 'P'
        V = varpar(1);
end

res = [yleft(2)-0.;   %tns=0
    yleft(4)-0.;      %psi=0
    yleft(5)-0.;      %x=0
    yleft(6)-0.;      %z=0 
    yleft(7)-0.;      %v=0
    yleft(8)-0.;      %a=0
    yleft(9)-0.;      %s0=0
    yleft(10)-0.;     %q=0
    yleft(12)-0.;     %I1=0
    yright(4)-pi;     %psi=pi
    yright(5)-0.;     %x=0
    yright(7)-V;      %v=1
    yright(9)-1;      %s0=1
    yright(10)-0;     %q=0
    yright(12)-0.;    %I1=0
    yleft(11)-0.;     %dsq=0 on SP
    yright(11)-0];    %dsq=0 on NP
       
end
