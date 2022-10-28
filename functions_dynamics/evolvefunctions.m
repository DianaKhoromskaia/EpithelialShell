function [C1new, dsC1new, C2new, Cnew, dsCnew, Xnew, Psinew, Znew, Unew, dsUnew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew, xintegral, solnem, s0new, s0inv, Q] = evolvefunctions(svec, zeta, dszeta, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac, dszetac, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun, snewvec, sfun, t, dt, thalf, tsigma, L, Lnew, L0, eps1, eps2, etacb, npoints, solnem_old, lc, s0, optode_nem)
%% Lagrangian time evolution of the surface descriptors and fields, for times t>dt:    
    
    svec3 = [0. svec((svec > eps1)&(svec < (L-eps2))) L]; % grid on surface at time t
    svec2 = svec3(2:end-1); % svec3 without endpoints
    
    snew = [0. snewfun(svec2) Lnew]; % grid on displaced surface
    
    sgrid = snew;
    svec1 = linspace(0,Lnew,npoints);%uniform grid for derivatives
    h = Lnew/(npoints-1);
    
    Cnew0 = C(0) + dt*(mss(0)-2*kappa*C(0)-zetac(0)+zetacnem(0))/etacb;
    CnewL = C(L) + dt*(mss(L)-2*kappa*C(L)-zetac(L)+zetacnem(L))/etacb;
   
    %% Update shape on Lagrangian grid:
    Xnew = griddedInterpolant( snew, [0. (X(svec2) + dt*(sin(Psi(svec2)).*vn(svec2) + cos(Psi(svec2)).*vs(svec2))) 0.] , 'spline');
    Psinew = griddedInterpolant( snew, [0. (Psi(svec2) + dt*( -dsvn(svec2) + C2(svec2).*vs(svec2))) pi], 'spline');
    Znew = griddedInterpolant(snew, Z(svec3)+dt*(-cos(Psi(svec3)).*vn(svec3) + sin(Psi(svec3)).*vs(svec3)), 'spline');

    %% Update curvatures and their derivatives on Lagrangian grid:
    jac = ones(size(svec2)) + dt*(vn(svec2).*C2(svec2)+dsvs(svec2)); 
    Cnew = griddedInterpolant( snew, [Cnew0 (C(svec2) + dt*( (mss(svec2)-2*kappa*C(svec2)-zetac(svec2)+zetacnem(svec2))/etacb) ) CnewL], 'spline');
    dsCnew = griddedInterpolant( snew, [0. (dsC(svec2) + dt*( (tns(svec2)+2*cos(Psi(svec2)).*zetacnem(svec2)./X(svec2)-2*kappa*dsC(svec2)-dszetac(svec2)+dszetacnem(svec2))/etacb))./jac 0.], 'spline');

    C1new = griddedInterpolant(snew, [Cnew0/2 C1(svec2)+dt*(-dsvn(svec2).*cos(Psi(svec2))./X(svec2) - vn(svec2).*C1(svec2).*C1(svec2) + vs(svec2).*dsC1(svec2) ) CnewL/2], 'spline');  %
    C2new = griddedInterpolant(snew, Cnew(snew)-C1new(snew), 'spline');
    
    dsC1new = griddedInterpolant( snew, [0. cos(Psinew(snewfun(svec2))).*(C2new(snewfun(svec2))-C1new(snewfun(svec2)))./Xnew(snewfun(svec2)) 0.], 'spline');

    %% Update area strain on Lagrangian grid:
    Unew = griddedInterpolant(snew, U(svec3) + dt*((ones(size(svec3))+U(svec3)).*vkk(svec3)), 'spline');
    gradU = gradient(Unew(svec1),h);
    dsUnew = griddedInterpolant(svec1, [0 gradU(2:end-1) 0], 'spline');

    %% Calculate s0 on new shape:
    s0new = griddedInterpolant(snew, s0(svec3),'spline');
    s0inv = griddedInterpolant(s0(svec3), snew,'spline');
       
    %% Solve nematic Euler-Lagrange equations, with solution from previous time step as initial guess:
    sgrid_nem = snewfun(solnem_old.x);
    sgrid_nem = [0 sgrid_nem(2:(end-1)) Lnew];
    parguess = solnem_old.parameters;
    solinit_nem = bvpinit(sgrid_nem, @yguessfun_nematic, parguess, sfun, solnem_old);

    solnem = bvp4c( @ode_nematic, ...
        @bc_nematic, ...
        solinit_nem, optode_nem, Psinew, Xnew, Lnew, lc);
    solnem_sgrid = deval(solnem,sgrid);
    solnem_svec1 = deval(solnem,svec1);
    Qnew = solnem_sgrid(1,:);
    dsQnew = solnem_svec1(2,:);

    Q = griddedInterpolant(snew, Qnew, 'spline');

    %% Update active profiles:
    zetanemvec = zeros(size(sgrid));
    dszetanemvec = zeros(size(svec1));
    zetacnemvec = zeros(size(sgrid));
    dszetacnemvec = zeros(size(svec1));
    
    onevec = ones(size(sgrid));
    zetaNew =  griddedInterpolant(sgrid, (1-sigmoidal(t+dt,thalf,tsigma))*(zeta_const*onevec+zeta_par*sigmoidal(s0new(sgrid), zeta_la*L0, zeta_sigma*L0)), 'spline');
    zetacNew = griddedInterpolant(sgrid, (1-sigmoidal(t+dt,thalf,tsigma))*(zetac_const*onevec+zetac_par*sigmoidal(s0new(sgrid), zetac_la*L0, zetac_sigma*L0)), 'spline');
    
    % Define new zetanem(s) profile using recalculated q(s):
    switch zetanem_Profile
        case 'Gaussian'
            zetanemvec = Qnew.*superGaussian(s0new(sgrid), zetanem_la*L0, zetanem_par, zetanem_const, zetanem_sigma*L0, 1);
        case 'Sigmoidal'
            if zetanem_la==1
                zetanemvec = zetanem_const*Qnew;
            else
                zetanemvec = zetanem_const*Qnew+zetanem_par*Qnew.*sigmoidal(s0new(sgrid), zetanem_la*L0, zetanem_sigma*L0);
            end
    end
    
    % Define new zetacnem(s) profile using recalculated q(s):
    switch zetacnem_Profile
        case 'Gaussian'
            zetacnemvec = Qnew.*superGaussian(s0new(sgrid), zetacnem_la*L0, zetacnem_par, zetacnem_const, zetacnem_sigma*L0, 1);
        case 'Sigmoidal'
            if zetacnem_la==1
                zetacnemvec = zetacnem_const*Qnew;
            else
                zetacnemvec = zetacnem_const*Qnew+zetacnem_par*Qnew.*sigmoidal(s0new(sgrid), zetacnem_la*L0, zetacnem_sigma*L0);
            end
    end
    
    zetanemvec = (1-sigmoidal(t+dt,thalf,tsigma))*zetanemvec;
    zetacnemvec = (1-sigmoidal(t+dt,thalf,tsigma))*zetacnemvec;
    
    zetanemNew = griddedInterpolant(sgrid, zetanemvec, 'spline');       
    zetacnemNew = griddedInterpolant(sgrid, zetacnemvec, 'spline');
    
    % Profile gradients:
    dszetavec = gradient(zetaNew(svec1),h); 
    dszetacvec = gradient(zetacNew(svec1),h); 
    
    % Define dszetanem(s) profile using recalculated \partial_s q:
    if (zetanem_la==1)&&strcmp(zetanem_Profile,'Sigmoidal')
        dszetanemvec = zetanem_const*dsQnew;
        dszetanemvec = (1-sigmoidal(t,thalf,tsigma))*dszetanemvec;
    else
        dszetanemvec = gradient(zetanemNew(svec1),h);
    end
    
    % Define dszetacnem(s) profile using recalculated \partial_s q:
    if (zetacnem_la==1)&&strcmp(zetacnem_Profile,'Sigmoidal')
        dszetacnemvec = zetacnem_const*dsQnew;
        dszetacnemvec = (1-sigmoidal(t,thalf,tsigma))*dszetacnemvec;
    else
        dszetacnemvec = gradient(zetacnemNew(svec1),h);
    end
    
    dszetaNew = griddedInterpolant(svec1, [0 dszetavec(2:end-1) 0], 'spline');
    dszetacNew = griddedInterpolant(svec1, [0 dszetacvec(2:end-1) 0], 'spline'); 
    dszetanemNew = griddedInterpolant(svec1, [0 dszetanemvec(2:end-1) 0], 'spline'); 
    dszetacnemNew = griddedInterpolant(svec1, [0 dszetacnemvec(2:end-1) 0], 'spline');
       
    %% Update integral of the shape:
    xintegral = integral(@(s) Xnew(s), 0., Lnew);

end

