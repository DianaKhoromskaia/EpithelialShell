function JacobianCheck(sol, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt)
%% Checking Jacobian wrt functions:

    %% centre of interval:
    pos=0.25;
    vec0 = deval(sol, pos);
    Jacpert = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    
    epsjac = 0.1;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters,U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,1)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,2)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,3)
    contourf(Jacpert-Fpert);
    colorbar
    
    %% SP:
    pos=0.;
    vec0 = deval(sol, pos);
    Jacpert = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    epsjac = 0.1;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,4)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,5)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,6)
    contourf(Jacpert-Fpert);
    colorbar
    
    %% NP:
    pos=1;
    vec0 = deval(sol, pos);
    Jacpert = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    epsjac = 0.1;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(10)
    subplot(3,3,7)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    subplot(3,3,8)
    contourf(Jacpert-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10);
    for k = 1:10
        pert = zeros(size(vec0));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0+pert, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    subplot(3,3,9)
    contourf(Jacpert-Fpert);
    colorbar
    
%% Checking Jacobian wrt parameters:
    pos=0.25;
    vec0 = deval(sol, pos);
    [Jacpert, Jacpertvar] = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    epsjac = 0.1;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    figure(11)
    subplot(3,3,1)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    subplot(3,3,2)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
        %(Fpert1-Fpert2)./pert' - Jacpert
        %plot(Jacpert*pert-(Fpert1-Fpert2),'o'); hold on;
    end
    subplot(3,3,3)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    %% SP:
    pos=0.;
    vec0 = deval(sol, pos);
    [Jacpert, Jacpertvar] = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    epsjac = 0.1;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,4)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,5)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,6)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    %% NP:
    pos=1;
    vec0 = deval(sol, pos);
    [Jacpert, Jacpertvar] = fjac(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    F0 = ode(pos, vec0, sol.parameters, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
    
    epsjac = 0.1;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,7)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.01;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,8)
    contourf(Jacpertvar-Fpert);
    colorbar
    
    epsjac = 0.001;
    Fpert = zeros(10,length(sol.parameters));
    for k = 1:length(sol.parameters)
        pert = zeros(size(sol.parameters));
        pert(k,1) = epsjac;
        F1 = ode(pos, vec0, sol.parameters+pert, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, dt);
        Fpert(:,k)=((F1-F0)./abs(F0))/epsjac;
    end
    subplot(3,3,9)
    contourf(Jacpertvar-Fpert);
    colorbar


    
    
end