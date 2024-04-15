function [sint, Vint, Par, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, dsmss, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, SolFound] = forcebalance_dthalf(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, delt, FixedPar, P0, thalf_P, tsigma, varargin)
SolFound = true;

if t==0
    sgridinit = varargin{1};
    tss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
    tns = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
    mss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
    
    sgrid = sgridinit; 

    %with fc as second parameter:
    if strcmp(FixedPar,'V')
        parguess = [0 0]; 
    elseif strcmp(FixedPar,'P')
        parguess = [0 0];
    end
    solinit = bvpinit(sgrid, @yguessfun_steadystate, parguess, mss, tss, tns);
    
elseif t>0
    sol = varargin{1};
    sfun = varargin{2};
    snewfun = varargin{3};
    
    solold = sol;
    
    sgrid = snewfun(solold.x);
    sgrid = [0. sgrid(2:(end-1)) L];
    
    fint = integral(@(s) X(s).*fext(s), 0, L);
    
    parguess = solold.parameters;
    
    solinit = bvpinit(sgrid, @yguessfun, parguess, solold, sfun);
end

try
    sol = bvp4c( @ode, ...
                 @bc, ...
                 solinit,optode, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, t, P0, thalf_P, tsigma);
catch ME
    ME
    disp(strcat('could not find solution at time t=',num2str(t)));
    SolFound = false;
    sol = solold; %this works only if t>0
end
    
sint = sol.x;
Vint = sol.y;
derivatives = sol.yp;
Par = sol.parameters;
dV = sol.y(6,end);
dX0 = sol.y(7,end);

% figure(30)
% plot(sint, Vint,'.-');
% axis tight

%% saving new arc length with dt=delt
Lnew_dt = L + delt*Vint(9,end);
eps1new_dt = eps1abs*Lnew_dt;

sfun_dt = griddedInterpolant( [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt], [0. sint(2:(end-1)) L], 'spline'); % interpolant of s(snew)
snewfun_dt = griddedInterpolant( [0. sint(2:(end-1)) L], [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt], 'spline'); % interpolant of snew(s)
snewvec_dt = [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt]; % vector snew

% 
% Lnew_dt = Vint(9,end);
% eps1new_dt = eps1abs*Lnew_dt;
% 
% sfun_dt = griddedInterpolant( Vint(9,:), sint, 'spline'); % interpolant of s(snew)
% snewfun_dt = griddedInterpolant( sint, Vint(9,:), 'spline'); % interpolant of snew(s)
% snewvec_dt = Vint(9,:); % vector snew

%% interpolants of functions and derivatives
%tss = griddedInterpolant(sint, Vint(1,:), 'spline');
dsvs = griddedInterpolant(sint, Vint(1,:), 'spline');
vn = griddedInterpolant(sint, Vint(2,:) , 'spline');
dsvn = griddedInterpolant(sint, Vint(3,:) , 'spline');
mss = griddedInterpolant(sint, Vint(4,:), 'spline');
dsmss = griddedInterpolant(sint, Vint(5,:), 'spline');
vs = griddedInterpolant(sint, Vint(8,:), 'spline');
%vkk = griddedInterpolant(sint, (Vint(1,:)+zetanem(sint)-2*K*U(sint)-zeta(sint))/etab, 'spline');
vkk = griddedInterpolant(sint, Vint(1,:)+C(sint).*Vint(2,:)+[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)], 'spline');
tss = griddedInterpolant(sint, 2*K*U(sint)+zeta(sint)-zetanem(sint)-2*kappa*C(sint).*(C2(sint)-0.5*C(sint))+(etab+eta)*Vint(1,:)+(etab*C(sint)+eta*(C2(sint)-C1(sint))).*Vint(2,:)+(etab-eta)*[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)], 'spline');
%dsvs = griddedInterpolant(sint, [((Vint(1,1)+zetanem(0)-2*K*U(0)-zeta(0))/etab-C(0)*Vint(2,1))/2 (Vint(1,2:end-1)+zetanem(sint(2:end-1))-2*K*U(sint(2:end-1))-zeta(sint(2:end-1)))/etab-cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1))-C(sint(2:end-1)).*Vint(2,2:end-1) ((Vint(1,end)+zetanem(L)-2*K*U(L)-zeta(L))/etab-C(L)*Vint(2,end))/2], 'spline');

ds2vn = griddedInterpolant(sint, derivatives(3,:), 'spline');

end

