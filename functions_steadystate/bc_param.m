function res = bc_param(yleft, yright, varpar, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar)
% residual function for solution with step-profile, in parametric mode

%v = (t,   tns,  C  , psi,    x,    z,     v,     a,  s0,   q,  dsq,  I1,  L)
%  = (1,     2,  3  ,   4,    5,    6,     7,     8,   9,  10,  11,   12, 13)

switch FixedPar
    case 'V'
        V = fixedpar;
    case 'P'
        V = varpar(1);
end

switch ControlPar
    case 'Bending'
        aM = guess(1)-dot([varpar(1) varpar(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
    case 'BendingNematic'
        zetacnemred = guess(1)-dot([varpar(1) varpar(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
end

res = [yleft(2,1)-0.;
    yleft(4,1)-0.;
    yleft(5,1)-0.;
    yleft(6,1)-0.;
    yleft(7,1)-0.;
    yleft(8,1)-0.;
    yleft(9,1)-0.;
    yleft(10,1)-0.;
    yleft(12,1)-0.;
    yleft(1:2,2)-yright(1:2,1); %matching
    yleft(3,2)-yright(3,1)-(aM-zetacnemred*0.5*(yright(10,1)+yleft(10,2)))/2; %curvature-jump on C_k^k
    yleft(4:12,2)-yright(4:12,1); %matching
    yright(4,2)-pi;
    yright(5,2)-0.;
    yright(7,2)-V; %for fixed or changing volume
    yright(9,2)-1;
    yright(10,2)-0;
    yright(12,2)-0.;
    yright(9,1)-la; %boundary on s0
    yleft(11,1)-0.; %dsq=0 at SP
    yright(11,2)-0]; %dsq=0 at NP

end
