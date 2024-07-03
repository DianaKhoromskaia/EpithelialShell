function [C1new, dsC1new, C2new, Cnew, dsCnew, Xnew, Psinew, Znew, Unew, dsUnew, kappaNew, dskappaNew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew, xintegral, solnem, s0new, s0inv, Q] = evolvefunctions_t0(svec, zeta_controls, zeta_profiles, zeta_implementation_types, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac, dszetac, zetacnem, dszetacnem, vs, dsvs, vkk, vn, dsvn, mss, tns, U, C1, dsC1, C2, C, C0, dsC, kappa, dskappa, X, Psi, Z, snewfun, sfun, t, dt, tsigma, L, Lnew, L0, eps1, eps2, etacb, npoints, dsw0, lc, Q, dsQ, s0, optode_nem, N_regions, zetasrect, kappa0, write94)
%Lagrangian update of everything:
    
    svec3 = [0. svec((svec > eps1)&(svec < (L-eps2))) L];%svec;%
    %svec3 = svec;
    svec2 = svec3(2:end); % vector of svec3 without endpoints
    
    snew = [0. snewfun(svec3(2:end-1)) Lnew];%snewvec;%
    
    sgrid = snew;
    sgridnem = [0. snew((snew > eps1)&(snew < (Lnew-eps2))) Lnew];%snew;
    svec1 = linspace(0,Lnew,npoints);%uniform grid for derivatives
    ds = Lnew/(npoints-1);

    Cnew0 = C(0) + dt*(mss(0)-2*kappa(0)*(C(0)-C0)-zetac(0)+zetacnem(0))/etacb;

    %% Update shape on Lagrange grid: 
     Xnew = griddedInterpolant( snew, [0. (X(svec2) + dt*(sin(Psi(svec2)).*vn(svec2) + cos(Psi(svec2)).*vs(svec2)))] , 'spline');
     Psinew = griddedInterpolant( snew, [0. (Psi(svec2) + dt*( -dsvn(svec2) + C2(svec2).*vs(svec2)))], 'spline');
     Znew = griddedInterpolant(snew, Z(svec3)+dt*(-cos(Psi(svec3)).*vn(svec3) + sin(Psi(svec3)).*vs(svec3)), 'spline');

    %% Update curvatures on Lagrangian grid:
    %jac = sqrt(ones(size(svec2)) + 2*dt*(vn(svec2).*C2(svec2)+dsvs(svec2)) + (dt)^2*((vn(svec2).*C2(svec2)+dsvs(svec2)).^2+(dsvn(svec2)-C2(svec2).*vs(svec2)).^2));%abs(ones(size(svec2)) + dt*(vn(svec2).*C2(svec2)+dsvs(svec2))); 
    jac = ones(size(svec2)) + dt*(vn(svec2).*C2(svec2)+dsvs(svec2)); 
    Cnew = griddedInterpolant( snew, [Cnew0 (C(svec2) + dt*( (mss(svec2)-2*kappa(svec2).*(C(svec2)-C0)-zetac(svec2)+zetacnem(svec2))/etacb) )], 'spline');
    dsCnew = griddedInterpolant( snew, [0. (dsC(svec2) + dt*( (tns(svec2)+2*cos(Psi(svec2)).*zetacnem(svec2)./X(svec2)-2*kappa(svec2).*dsC(svec2)-2*dskappa(svec2).*(C(svec2)-C0)-dszetac(svec2)+dszetacnem(svec2))/etacb))./jac], 'spline');
    
    C1new = griddedInterpolant(snew, [Cnew0/2 C1(svec2)+dt*(-dsvn(svec2).*cos(Psi(svec2))./X(svec2) -vn(svec2).*C1(svec2).*C1(svec2) + vs(svec2).*dsC1(svec2) )], 'spline');
    C2new = griddedInterpolant(snew, Cnew(snew)-C1new(snew), 'spline');
    
    dsC1new = griddedInterpolant( snew, [0. cos(Psinew(snewfun(svec2))).*(C2new(snewfun(svec2))-C1new(snewfun(svec2)))./Xnew(snewfun(svec2))], 'spline');
    %% update U and stretches in Lagrangian view:
    Unew = griddedInterpolant(snew, U(svec3) + dt*((ones(size(svec3))+U(svec3)).*vkk(svec3)), 'spline');
    gradU = gradient(Unew(svec1),ds);
    dsUnew = griddedInterpolant(svec1, [0 gradU(2:end)], 'spline');
    %fsnew = griddedInterpolant(snew, fs(svec3) + dt*(dsvs(svec3)+C2(svec3).*vn(svec3))./fs(svec3), 'spline');

    %% s0 on new shape, if la<1:
    s0new = griddedInterpolant(snew, s0(svec3),'spline');
    s0inv = griddedInterpolant(s0(svec3), snew,'spline');
       
    %% solve nematic Euler-Lagrange equations
    parguess = dsw0;
    solinit_nem = bvpinit(sgridnem, @yguessfun_nematic_t0, parguess, sfun, Q, dsQ);
    
    solnem = bvp4c( @ode_nematic, ...
        @bc_nematic, ...
        solinit_nem, optode_nem, Psinew, Xnew, Lnew, lc);
    solnem_grid = deval(solnem,sgridnem);
    Qnew = solnem_grid(1,:);
    
    Q = griddedInterpolant(sgridnem, Qnew, 'spline');

    [kappaNew, dskappaNew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew] = actualiseprofiles(zeta_controls, zeta_profiles, zeta_implementation_types, zeta_consts, zeta_las, zeta_facs, zeta_sigmas, zeta_thalfs, Qnew, svec1, sgrid, sgridnem, s0new, N_regions, L0, zetasrect, t, dt, tsigma, ds, kappa0, write94);
      
    xintegral = integral(@(s) Xnew(s), 0., Lnew);
end

