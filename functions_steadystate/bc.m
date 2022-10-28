function res = bc(yleft, yright, varpar, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar)
% residual function for solution with step-profile, in non-parametric mode

%v = (t,   tns,  C  , psi,    x,    z,     v,     a,  s0,   q,  dsq,  I1)
%  = (1,     2,  3  ,   4,    5,    6,     7,     8,   9,  10,  11,   12)
switch FixedPar
    case 'V'
        V = fixedpar;
    case 'P'
        V = varpar(1);
end

res = [yleft(2,1)-0.;   %tns=0
    yleft(4,1)-0.;      %psi=0
    yleft(5,1)-0.;      %x=0
    yleft(6,1)-0.;      %z=0 
    yleft(7,1)-0.;      %v=0
    yleft(8,1)-0.;      %a=0
    yleft(9,1)-0.;      %s0=0
    yleft(10,1)-0.;     %q=0
    yleft(12,1)-0.;     %I1=0
    yleft(1:2,2)-yright(1:2,1); %matching tss and tns
    yleft(3,2)-yright(3,1)-(aM-zetacnemred*0.5*(yright(10,1)+yleft(10,2)))/2; %curvature-jump on C^k_k
    yleft(4:12,2)-yright(4:12,1); %matching psi,x,z,v,a,s0,q,dsq,I1
    yright(4,2)-pi;     %psi=pi
    yright(5,2)-0.;     %x=0
    yright(7,2)-V;      %v=1
    yright(9,2)-1;      %s0=1
    yright(10,2)-0;     %q=0
    yright(12,2)-0.;    %I1=0
    yright(9,1)-la;     %s0 at boundary
    yleft(11,1)-0.;     %dsq=0 at SP
    yright(11,2)-0];     %dsq=0 at NP  

end
