function [C1new, dsC1new, C2new, Cnew, dsCnew, Xnew, Psinew, Znew, Unew, dsUnew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew, xintegral, s0, s0inv] = evolvefunctions_continued(svec, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun, snewvec, sfun, t, dt, L, Lnew, L0, eps1, eps2, etacb, npoints, s0)
%Lagrangian update of everything:
    
    svec3 = [0. svec((svec > eps1)&(svec < (L-eps2))) L];%svec;%
    %svec3 = svec;
    svec2 = svec3(2:end-1); % vector of svec3 without endpoints
    
    snew = [0. snewfun(svec2) Lnew];%snewvec;%
    
    sgrid = snew;
    sgridnem = [0. snew((snew > eps1)&(snew < (Lnew-eps2))) Lnew];%snew;
    svec1 = linspace(0,Lnew,npoints);%uniform grid for derivatives
    h = Lnew/(npoints-1);
    
    Cnew0 = C(0) + dt*(mss(0)-2*kappa*C(0)-zetac(0)+zetacnem(0))/etacb;
    CnewL = C(L) + dt*(mss(L)-2*kappa*C(L)-zetac(L)+zetacnem(L))/etacb;

    %% Update shape on Lagrange grid: 
     Xnew = griddedInterpolant( snew, [0. (X(svec2) + dt*(sin(Psi(svec2)).*vn(svec2) + cos(Psi(svec2)).*vs(svec2))) 0.] , 'spline');
     Psinew = griddedInterpolant( snew, [0. (Psi(svec2) + dt*( -dsvn(svec2) + C2(svec2).*vs(svec2))) pi], 'spline');
     Znew = griddedInterpolant(snew, Z(svec3)+dt*(-cos(Psi(svec3)).*vn(svec3) + sin(Psi(svec3)).*vs(svec3)), 'spline');

    %% Update curvatures on Lagrangian grid:
    jac = sqrt(ones(size(svec2)) + 2*dt*(vn(svec2).*C2(svec2)+dsvs(svec2)) + (dt)^2*((vn(svec2).*C2(svec2)+dsvs(svec2)).^2+(dsvn(svec2)-C2(svec2).*vs(svec2)).^2));%abs(ones(size(svec2)) + dt*(vn(svec2).*C2(svec2)+dsvs(svec2))); 
    %jac = ones(size(svec2)) + dt*(vn(svec2).*C2(svec2)+dsvs(svec2)); 
    Cnew = griddedInterpolant( snew, [Cnew0 (C(svec2) + dt*( (mss(svec2)-2*kappa*C(svec2)-zetac(svec2)+zetacnem(svec2))/etacb) ) CnewL], 'spline');
    dsCnew = griddedInterpolant( snew, [0. (dsC(svec2) + dt*( (tns(svec2)+2*cos(Psi(svec2)).*zetacnem(svec2)./X(svec2)-2*kappa*dsC(svec2)-dszetac(svec2)+dszetacnem(svec2))/etacb))./jac 0.], 'spline');
    %dsCnewvec = gradient(Cnew(svec1),h);
    %dsCnew = griddedInterpolant( svec1, [0. dsCnewvec(2:end-1) 0.], 'spline');
    
    C1new = griddedInterpolant(snew, [Cnew0/2 C1(svec2)+dt*(-dsvn(svec2).*cos(Psi(svec2))./X(svec2) -vn(svec2).*C1(svec2).*C1(svec2) + vs(svec2).*dsC1(svec2) ) CnewL/2], 'spline');
    C2new = griddedInterpolant(snew, Cnew(snew)-C1new(snew), 'spline');
    
    dsC1new = griddedInterpolant( snew, [0. cos(Psinew(snewfun(svec2))).*(C2new(snewfun(svec2))-C1new(snewfun(svec2)))./Xnew(snewfun(svec2)) 0.], 'spline');
    %% update U and stretches in Lagrangian view:
    Unew = griddedInterpolant(snew, U(svec3) + dt*((ones(size(svec3))+U(svec3)).*vkk(svec3)), 'spline');
    gradU = gradient(Unew(svec1),h);
    dsUnew = griddedInterpolant(svec1, [0 gradU(2:end-1) 0], 'spline');
    %fsnew = griddedInterpolant(snew, fs(svec3) + dt*(dsvs(svec3)+C2(svec3).*vn(svec3))./fs(svec3), 'spline');

    %% s0 on new shape, if la<1:
    s0 = griddedInterpolant(snew, s0(svec3),'spline');
    s0inv = griddedInterpolant(s0(svec3), snew,'spline');
     
    %% Checking geometric relations:
%     int1 = trapz(svec1,gradient(Xnew(svec1),h)-cos(Psinew(svec1)));
%     
%     figure(50)
%     subplot(1,2,1)
%     plot(svec1, gradient(Xnew(svec1),h)-cos(Psinew(svec1)), 'b-');
%     legend('\partial_s x-\cos(\psi)');
%     subplot(1,2,2)
%     hold on;
%     plot(t, int1,'r.')
%     
%     
%     %% solve nematic Euler-Lagrange equations
%     %if strcmp(ControlPar,'Nematic')||strcmp(ControlPar,'BendingNematic')
%         parguess = [dsw0 dswL];
%         solinit_nem = bvpinit(sgridnem, @yguessfun_nematic_t0, parguess, sfun, Q, dsQ);
%         
%         solnem = bvp4c( @ode_nematic, ...
%             @bc_nematic, ...
%             solinit_nem, optode_nem, Psinew, Xnew, Lnew, lc);
%         solnem_grid = deval(solnem,sgridnem);
%         solnem_svec1 = deval(solnem,svec1);
%         Qnew = solnem_grid(1,:);
%         dsQnew = solnem_svec1(2,:);
%         
%         Q = griddedInterpolant(sgridnem, Qnew, 'spline');
%      %else
%      %    solnem = 0;
%      %end
    
 % plot profiles to check:
%  figure(2)
%  subplot(2,2,1)
%  plot(svec3, zeta(svec3), svec3, dszeta(svec3));
%  legend('\zeta','\partial_s\zeta');
%  subplot(2,2,2)
%  plot(svec3, zetac(svec3), svec3, dszetac(svec3));
%  legend('\zeta_c','\partial_s\zeta_c');
%  subplot(2,2,3)
%  plot(svec3, zetanem(svec3), svec3, dszetanem(svec3));
%  legend('\zeta_{n}','\partial_s\zeta_{n}');
%  subplot(2,2,4)
%  plot(svec3, zetacnem(svec3), svec3, dszetacnem(svec3));
%  legend('\zeta_{cn}','\partial_s\zeta_{cn}');

    %% update profiles:
    %zetavec = zeros(size(sgrid));
    %dszetavec = zeros(size(svec1));
    %zetacvec = zeros(size(sgrid));
    %dszetacvec = zeros(size(svec1));
    zetanemvec = zeros(size(sgrid));
    dszetanemvec = zeros(size(svec1));
    zetacnemvec = zeros(size(sgrid));
    dszetacnemvec = zeros(size(svec1));
    
    zetaNew = griddedInterpolant(sgrid, zeta(svec3), 'spline');
    %zetacNew = griddedInterpolant(sgrid, zetac(svec3)-dt*vs(svec3).*dszeta(svec3), 'spline');
    zetacNew = griddedInterpolant(sgrid, ((1-sigmoidal(t+dt,0.01,0.002))/(1-sigmoidal(t,0.01,0.002)))*zetac(svec3), 'spline');
    
%     % define new zetanem(s) profile: (nematic tension)
%     switch zetanem_Profile
%         case 'Gaussian'
%             zetanemvec = Qnew.*superGaussian(s0(sgridnem), zetanem_la*L0, zetanem_par, zetanem_const, zetanem_sigma*L0, 1);
%         case 'Sigmoidal'
%             if zetanem_la==1
%                 zetanemvec = zetanem_const*Qnew;
%             else
%                 zetanemvec = zetanem_const*Qnew+zetanem_par*Qnew.*sigmoidal(s0(sgridnem), zetanem_la*L0, zetanem_sigma*L0);
%             end
%     end
%     
%     % define new zetacnem(s) profile: (nematic bending moment)
%     switch zetacnem_Profile
%         case 'Gaussian'
%             zetacnemvec = Qnew.*superGaussian(s0(sgridnem), zetacnem_la*L0, zetacnem_par, zetacnem_const, zetacnem_sigma*L0, 1);
%         case 'Sigmoidal'
%             if zetacnem_la==1
%                 zetacnemvec = zetacnem_const*Qnew;
%             else
%                 zetacnemvec = zetacnem_const*Qnew+zetacnem_par*Qnew.*sigmoidal(s0(sgridnem), zetacnem_la*L0, zetacnem_sigma*L0);
%             end
%     end
    
    zetanemNew = griddedInterpolant(sgridnem, zetanemvec, 'spline');       
    zetacnemNew = griddedInterpolant(sgridnem, zetacnemvec, 'spline');

    % profile gradients:
    dszetavec = gradient(zetaNew(svec1),h); 
    dszetacvec = gradient(zetacNew(svec1),h); 
    
%     % define dszetanem(s) profile: (nematic tension)
%     if (zetanem_la==1)&&strcmp(zetanem_Profile,'Sigmoidal')
%         dszetanemvec = zetanem_const*dsQnew;
%     else
%         dszetanemvec = gradient(zetanemNew(svec1),h);
%     end
%     
%     % define dszetacnem(s) profile: (nematic bending)
%     if (zetacnem_la==1)&&strcmp(zetacnem_Profile,'Sigmoidal')
%         dszetacnemvec = zetacnem_const*dsQnew;
%     else
%         dszetacnemvec = gradient(zetacnemNew(svec1),h);
%     end
    
    dszetaNew = griddedInterpolant(svec1, [0 dszetavec(2:end-1) 0], 'spline');
    dszetacNew = griddedInterpolant(svec1, [0 dszetacvec(2:end-1) 0], 'spline'); 
    dszetanemNew = griddedInterpolant(svec1, [0 dszetanemvec(2:end-1) 0], 'spline'); 
    dszetacnemNew = griddedInterpolant(svec1, [0 dszetacnemvec(2:end-1) 0], 'spline');    
    
     % plot profiles to check:
%  figure(2)
%  subplot(2,2,1)
%  plot(svec1, zetaNew(svec1), svec1, dszetaNew(svec1));
%  legend('\zeta','\partial_s\zeta');
%  subplot(2,2,2)
%  plot(svec1, zetacNew(svec1), svec1, dszetacNew(svec1));
%  legend('\zeta_c','\partial_s\zeta_c');
%  subplot(2,2,3)
%  plot(svec1, zetanemNew(svec1), svec1, dszetanemNew(svec1));
%  legend('\zeta_{n}','\partial_s\zeta_{n}');
%  subplot(2,2,4)
%  plot(svec1, zetacnemNew(svec1), svec1, dszetacnemNew(svec1));
%  legend('\zeta_{cn}','\partial_s\zeta_{cn}');
      
    xintegral = integral(@(s) Xnew(s), 0., Lnew);
end

