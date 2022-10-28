% Dynamical deformations of epithelial shells with active tensions and bending moments
% (isotropic and nematic).
% Code for the publication "Active morphogenesis of patterned epithelial
% shells".

% Authors: Diana Khoromskaia and Guillaume Salbreux
% Date last edited: 27/10/2022

function Epithelia_Dynamics(dt, tmax, FixedPar, ControlPar, param_const, param_la, param_fac, dtmax, Adaptive, xi, tol)
% Input arguments:
% dt - initial time step
% tmax - end time
% FixedPar - set to 'V' or 'P' for fixed or free volume, respectively
% ControlPar - control parameter, i.e. choice of active effect to be considered: 'Tension', 'Bending',
% 'Nematic, 'BendingNematic'
% param_const, param_la, param_fac - parameters defining the active
% profile, e.g. for 'Bending' those correspond to \zeta_c^0, l_a/L_0, and
% \delta\zeta_c in Eq. (19)
% dtmax - maximal time step
% Adaptive - adaptive time stepping 'On' or 'Off'
% xi - normal friction coefficient; set to 0, not required here. 
% tol - relative tolerance on adaptive time step

addpath('./functions_dynamics');

%% set physical parameters:
R0 = 1; % radius of initial sphere
h = 0.1; % thickness of epithelium in units of R0
eta = 1; % shear viscosity (in units of 2d shear viscosity of epithelium)
etab = 1; % bulk viscosity (in units of 2d shear viscosity of epithelium)
etacb = (h^2)*etab; % bulk bending viscosity (in units of eta/R_0^2)
etap = 1e-4; % pressure viscosity for the case FixedPar='P' (in units of eta/R_0^4)
kappa = 1; % bending modulus (eqns are normalised by kappa)
%xi = 0; %normal friction coefficient
K = 1000; % bulk elastic modulus (in units of kappa/R_0^2)
lc = 0.1*R0; % nematic length scale in (units of R0)

Plotting = 'On'; % set to 'On' or 'Off' to save intermediate results as plots/movies (to enable saving as movies adjust plotting.m)

%% set numerical parameters:
% set cut-off for nematic equations at SP and at NP of the surface in units of interval length:
eps1abs = 1e-4;
eps2abs = 1e-4;
% relative and absolute error bounds for bvp-solver:
epsoderel = 1e-4;
epsodeabs = 1e-6;
% options for bvp-solvers: -- set 'stats' to 'off' if details of solver output not required
optode = bvpset('RelTol', epsoderel, 'AbsTol', epsodeabs, 'stats','on', 'Vectorized', 'on','FJacobian',@fjac,'BCJacobian',@bcjac,'NMax', 1e+5); 
optode_nem = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'stats','on','Vectorized', 'on','FJacobian',@fjac_nematic,'BCJacobian',@bcjac_nematic,'NMax', 1e+5);

% set remaining time stepping parameters:
t = 0;
% parameters for t-dependent sigmoid: 
thalf = 0.01;  % mu in Eq (18)
tsigma = 0.002; % sigma in Eq (18)

%% sizes of different grids
nplot = 1e+3;       % for plotting and saving to file
npoints = 1e+6;     % for initialising functions on sphere and for calculating derivatives

%% initialise spherical shape:
z0 = 0.; %offset in z-direction
L = pi*R0; % current perimeter length
L0 = pi*R0; % initial perimeter length
[C1, C2, C, dsC1, dsC2, dsC, Psi, X, Z, X0, xintegral, svec1, U, dsU, Q, dsQ, dsw0, dswL, fs, s0, Qgrid] = initialisesphere(L0, R0, z0, npoints);
eps1 = eps1abs*L;
eps2 = eps2abs*L;
V0 = (4/3)*pi*R0^3; % initial volume
V = V0; % current volume

%if strcmp(ControlPar,'Nematic')||strcmp(ControlPar,'BendingNematic')
%    sgrid = Qgrid;
%else
    sgrid = linspace(0,L0,200);
%end

% save initial shape for plotting:
sinit = svec1;
Xinit = X;
Zinit = Z;
Uinit = U;

%% initilise 'active' profiles:
% ..._Profile - set to 'Sigmoidal' or 'Gaussian'
% ..._const - constant offset
% ..._la - position of half-maximum or peak, resp.
% ..._sigma - controls sharpness of sigmoid or width of Gaussian
% ..._par - prefactor for profile
% if 'Sigmoidal' and ..._la=1 then profile is flat: profile= ..._const
zeta_Profile = 'Sigmoidal';
zeta_const = 0;
zeta_la = 0;
zeta_sigma = 0.005;
zeta_par = 0;
zetac_Profile = 'Sigmoidal';
zetac_const = 0;
zetac_la = 0.;
zetac_sigma = 0.005;
zetac_par = 0;
zetanem_Profile = 'Sigmoidal';
zetanem_const = 0;
zetanem_la = 0;
zetanem_sigma = 0.005;
zetanem_par = 0;
zetacnem_Profile = 'Sigmoidal';
zetacnem_const = 0;
zetacnem_la = 0;
zetacnem_sigma = 0.005;
zetacnem_par = 0;

switch ControlPar 
    case 'Tension'
        zeta_const = param_const;
        zeta_la = param_la;
        zeta_par = param_fac;
    case 'Bending'
        zetac_const = param_const;
        zetac_la = param_la;
        zetac_par = param_fac;
    case 'Nematic'
        zetanem_const = param_const;
        zetanem_la = param_la;
        zetanem_par = param_fac;
    case 'BendingNematic'
        zetacnem_const = param_const;
        zetacnem_la = param_la;
        zetacnem_par = param_fac;
end
        
[zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, dszetacnem] = initialiseprofiles(zeta_Profile, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac_Profile, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem_Profile, zetanem_const, zetanem_la, zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, sgrid, svec1, L, npoints, L0, s0, Q, dsQ, t, thalf, tsigma);

%% initialise pressure and output directory, depending on FixedPar
if strcmp(FixedPar,'V')
    P0 = 2*zeta_const/R0; % Laplace pressure, initial guess for the bvp-solver
    dir1 = strcat('Dynamics_',ControlPar,'_V_la=',num2str(param_la),'_const=',num2str(param_const),'_fac=',num2str(param_fac),'_xi=',num2str(xi),'_tol=',num2str(tol));%datestr(now,30));
elseif strcmp(FixedPar,'P')
    P0 = 2*zeta_const/R0; % pressure is initialised to this value, and should return to it in steady state
    dir1 = strcat('Dynamics_',ControlPar,'_P_la=',num2str(param_la),'_const=',num2str(param_const),'_fac=',num2str(param_fac),'_xi=',num2str(xi),'_tol=',num2str(tol));%,datestr(now,30));
end
P = P0;

%% output files:
mkdir(dir1);
cd(dir1);

% define special point to evaluate functions at for saving, e.g. equator or
% pole.
if param_la==1
    seval=0.5*L;
else
    seval=0;
end

filename1 = 'observables.dat';
fileID = fopen(filename1,'w');
fprintf(fileID, '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n', 't', 'dt', 'v_n(seval)', 'P', 'X(seval)', 'C2(seval)','pole-pole', 'area', 'tcomp [sec]', 'L', 'V','X0', 'nmesh','intersect', 'seval','la','f_c','I(L)');
formatSpec = '%e \t %e \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %i \t %i \t %2.1f \t %3.2f \t %6.5f \t %9.8f \n';
fprintf(fileID, formatSpec, t, dt, 0., P0, X(seval), C2(seval), Z(L0)-Z(0.), 2*pi*xintegral, 0., L0, V0, X0, 0., 0., seval, param_la, 0, 0);

svecuni = linspace(0., L, nplot);%uniform grid for saving

% save functions of initial shape:
filename2 = 'curvatures.dat';
dlmwrite(filename2, svecuni, 'precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C1(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C2(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, dsC1(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, dsC(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');

filename3 = 'x.dat';
dlmwrite(filename3, X(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');
filename4 = 'z.dat';
dlmwrite(filename4, Z(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');
filename41 = 'psi.dat';
dlmwrite(filename41, Psi(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');

% save output of bvp-solver:
filename5 = 'bvpsolution.dat';
filename6 = 'vs.dat';
filename7 = 'vn.dat';
filename9 = 'profiles.dat'; 
filename10 = 'u.dat';
filename11 = 'q.dat';
filename12 = 'tss.dat';
filename13 = 's0.dat';

vidObj1 = VideoWriter('shape_sequence.avi');%,'MPEG-4');
vidObj1.Quality = 25;
vidObj1.FrameRate = 5;
open(vidObj1);

% time step for saving into data file:
delt_data = 10; %every 10-th step
n = 1; %number of saved time steps
i = 0; %number of time steps
error = 0; % initialise relative error for time stepping

%external force magnitude along z-direction: DO NOT USE THIS FOR NOW
fextmag = 0;
fext = griddedInterpolant(svec1, fextmag*svec1, 'spline');

%% save parameters into file
T1 = table({'FixedPar';'ControlPar';'param_const';'param_la'; 'param_fac';'zeta_sigma';'R0';'h';'eta';'etab';'etacb';'etap';'kappa';'xi';'K';'lc'},...
    [{FixedPar};{ControlPar};param_const;param_la;param_fac;zeta_sigma;R0;h;eta;etab;etacb;etap;kappa;xi;K;lc],...
    {'fixed volume (V) or free volume (P)';'control parameter';'const. offset';'size of active region at t=0, in units of L';'magnitude of active region';'sigma for Sigmoidal profile';'radius';'thickness';'shear viscosity';'bulk viscosity';'bending viscosity';'volume viscosity';'bending modulus';'friction';'elastic modulus';'nematic length, in units of R0'},'VariableNames',{'Parameter','Value','Description'});
writetable(T1,'parameters_physical.csv');

T2 = table({'Adaptive';'dt';'dtmax';'tmax';'tol';'epsoderel';'epsodeabs';'eps1abs';'nplot';'npoints';'delt_data'},...
    [{Adaptive};dt;dtmax;tmax;tol;epsoderel;epsodeabs;eps1abs;nplot;npoints;delt_data],...
    {'adaptive stepping';'initial time step';'max time step';'max total time';'rel. error tolerance for adaptive stepping';'rel. error for sover';'abs. error for solver';'cutoff at poles, used for nem. solver only';'size of uniform grid for saving initial shape';'size of uniform grid for profiles and gradients';'every n-th time step saved'},'VariableNames',{'Parameter','Value','Description'});
writetable(T2,'parameters_numerical.csv');

while t < tmax
    tnew = t
    
    %% solve force balance at time t (with step dt and dt/2):
    if t==0
        tic
        [svec, v, P1, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] ...
            = forcebalance(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, dt, Adaptive, FixedPar, sgrid);
        tcomp = toc;
    else
        tic
        [svec, v, P1, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] ...
            = forcebalance(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, dt, Adaptive, FixedPar, solold, sfun, snewfun);
        tcomp = toc;
    end
    
    if SolFound
        %% evolve system with time step dt:
        if t==0
            [C1_dt, dsC1_dt, C2_dt, C_dt, dsC_dt, X_dt, Psi_dt, Z_dt, U_dt, dsU_dt, zeta_dt, dszeta_dt, zetac_dt, dszetac_dt, zetanem_dt, dszetanem_dt, zetacnem_dt, dszetacnem_dt, xintegral_dt, solnem_dt, s0_dt, s0inv_dt, Q_dt] = ...
                evolvefunctions_t0(svec, zeta, dszeta, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac, dszetac, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, ...
                zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, ...
                U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dt, snewvec_dt, sfun_dt, t, dt, thalf, tsigma, L, Lnew_dt, L0, eps1, eps2, etacb, npoints, dsw0, dswL, lc, Q, dsQ, s0, optode_nem);
        else
            [C1_dt, dsC1_dt, C2_dt, C_dt, dsC_dt, X_dt, Psi_dt, Z_dt, U_dt, dsU_dt, zeta_dt, dszeta_dt, zetac_dt, dszetac_dt, zetanem_dt, dszetanem_dt, zetacnem_dt, dszetacnem_dt, xintegral_dt, solnem_dt, s0_dt, s0inv_dt, Q_dt] = ...
                evolvefunctions(svec, zeta, dszeta, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac, dszetac, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, ...
                zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, ...
                U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dt, snewvec_dt, sfun_dt, t, dt, thalf, tsigma, L, Lnew_dt, L0, eps1, eps2, etacb, npoints, solnem, lc, s0, optode_nem);
        end
        
        if strcmp(Adaptive,'On')
            %% evolve system with time step dt/2:
            if t==0
                [C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, X_dthalf, Psi_dthalf, Z_dthalf, U_dthalf, dsU_dthalf, zeta_dthalf, dszeta_dthalf, zetac_dthalf, dszetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, dszetacnem_dthalf, xintegral_dthalf, solnem_dthalf, s0_dthalf, s0inv_dthalf, Q_dthalf] = ...
                    evolvefunctions_t0(svec, zeta, dszeta, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac, dszetac, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, ...
                    zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, ...
                    U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dthalf, snewvec_dthalf, sfun_dthalf, t, dt/2, thalf, tsigma, L, Lnew_dthalf, L0, eps1, eps2, etacb, npoints, dsw0, dswL, lc, Q, dsQ, s0, optode_nem);
            else
                [C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, X_dthalf, Psi_dthalf, Z_dthalf, U_dthalf, dsU_dthalf, zeta_dthalf, dszeta_dthalf, zetac_dthalf, dszetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, dszetacnem_dthalf, xintegral_dthalf, solnem_dthalf, s0_dthalf, s0inv_dthalf, Q_dthalf] = ...
                    evolvefunctions(svec, zeta, dszeta, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac, dszetac, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem, dszetanem, zetacnem, dszetacnem, zetanem_Profile, zetanem_const, zetanem_la, ...
                    zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, ...
                    U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dthalf, snewvec_dthalf, sfun_dthalf, t, dt/2, thalf, tsigma, L, Lnew_dthalf, L0, eps1, eps2,etacb, npoints, solnem, lc, s0, optode_nem);
            end
            
            if strcmp(FixedPar,'V')
                P_dthalf = P1(1);
                V_dthalf = V0;
            elseif strcmp(FixedPar,'P')
                P_dthalf = -etap*P1(1);
                V_dthalf = V + 0.5*dt*dV;
            end
            L_dthalf = Lnew_dthalf;
            eps1_dthalf = eps1new_dthalf;
            X0_dthalf = X0 + 0.5*dt*dX0;
            %% solve force balance at time t+dt/2 with step dt/2:
            [svec_dt2half, v_dt2half, P1_dt2half, sol_dt2half, vs_dt2half, dsvs_dt2half, vkk_dt2half, tss_dt2half, vn_dt2half, dsvn_dt2half, mss_dt2half, tns_dt2half, ds2vn_dt2half, dV_dt2half, dX0_dt2half, Lnew_dt2half, eps1new_dt2half, sfun_dt2half, snewfun_dt2half, snewvec_dt2half, SolFound_dt2half] ...
                = forcebalance_dthalf(U_dthalf, dsU_dthalf, C1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, Psi_dthalf, X_dthalf, Z_dthalf, X0_dthalf, L_dthalf, L0, zeta_dthalf, dszeta_dthalf, zetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, eta, etab, etacb, etap, eps1abs, xintegral_dthalf, fext, kappa, K, xi, optode, t+dt/2, dt/2, FixedPar, sol, sfun_dthalf, snewfun_dthalf);
            
            %% evolve system with time step dt/2:
            [C1_dt2half, dsC1_dt2half, C2_dt2half, C_dt2half, dsC_dt2half, X_dt2half, Psi_dt2half, Z_dt2half, U_dt2half, dsU_dt2half, zeta_dt2half, dszeta_dt2half, zetac_dt2half, dszetac_dt2half, zetanem_dt2half, dszetanem_dt2half, zetacnem_dt2half, dszetacnem_dt2half, xintegral_dt2half, solnem_dt2half, s0_dt2half, s0inv_dt2half, Q_dt2half] = ...
                evolvefunctions(svec_dt2half, zeta_dthalf, dszeta_dthalf, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac_dthalf, dszetac_dthalf, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, dszetacnem_dthalf, zetanem_Profile, zetanem_const, zetanem_la, ...
                zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, vs_dt2half, dsvs_dt2half, vkk_dt2half, tss_dt2half, vn_dt2half, dsvn_dt2half, mss_dt2half, tns_dt2half, ds2vn_dt2half, ...
                U_dthalf, C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, kappa, X_dthalf, Psi_dthalf, Z_dthalf, snewfun_dt2half, snewvec_dt2half, sfun_dt2half, t+dt/2, dt/2, thalf, tsigma, L_dthalf, Lnew_dt2half, L0, eps1, eps2, etacb, npoints, solnem_dthalf, lc, s0_dthalf, optode_nem);
            
            %% compare results for adaptive time step after dt and 2*(dt/2):
            scompare = linspace(0,L,1000);
            scompare_dt = snewfun_dt(scompare);
            scompare_dt2half = snewfun_dt2half(snewfun_dthalf(scompare));

            % exclude values in Z and dsC from relative error calculation which are too close to zero, and exclude poles from X:
            cutoff = 1e-3;
            Z_dt_vec = Z_dt(scompare_dt);
            Z_dt2half_vec = Z_dt2half(scompare_dt2half);
            dsC_dt_vec = dsC_dt(scompare_dt);
            dsC_dt2half_vec = dsC_dt2half(scompare_dt2half);
            
            errorvec = [max(abs(X_dt(scompare_dt(2:end-1))-X_dt2half(scompare_dt2half(2:end-1)))./abs(X_dt(scompare_dt(2:end-1))));...
                max(abs((Z_dt_vec(abs(Z_dt_vec)>cutoff)-Z_dt2half_vec(abs(Z_dt_vec)>cutoff))./Z_dt_vec(abs(Z_dt_vec)>cutoff)));...
                max(abs((dsC_dt_vec(abs(dsC_dt_vec)>cutoff)-dsC_dt2half_vec(abs(dsC_dt_vec)>cutoff))./dsC_dt_vec(abs(dsC_dt_vec)>cutoff)))];

            error = max(errorvec);

            %% determine size of next time step:
            dtnext = min(0.9*dt*min([max([tol/error 0.5]) 1.5]), dtmax)
            
            if dtnext<1e-8 % increase ode-solver accuracy if time step goes below 1e-8
                epsoderel=1e-5;
            else
                epsoderel=1e-4;
            end
            optode = bvpset('RelTol', epsoderel, 'AbsTol', epsodeabs, 'stats','on','Vectorized', 'on','FJacobian',@fjac,'BCJacobian',@bcjac,'NMax', 1e+5); %'FJacobian',@fjac,'BCJacobian',@bcjac, 

            if dtnext < 1e-10 % stop, if time step goes below 1e-10
                disp('Time step decreased below 1e-10.')
                break;
            end
        else 
            dtnext = dt; 
        end
        
       if (error<tol)||(dtnext==1e-10)||strcmp(Adaptive,'Off') %% accept step only if tol satisfied, accept first step always;
            %% update remaining quantities with dt:
            C1 = C1_dt;
            dsC1 = dsC1_dt;
            C2 = C2_dt;
            C = C_dt;
            dsC = dsC_dt;
            X = X_dt;
            Psi = Psi_dt;
            Z = Z_dt;
            U = U_dt;
            dsU = dsU_dt;
            Q = Q_dt;
            zeta =zeta_dt;
            dszeta=dszeta_dt;
            zetac=zetac_dt;
            dszetac=dszetac_dt;
            zetanem=zetanem_dt;
            dszetanem=dszetanem_dt;
            zetacnem=zetacnem_dt;
            dszetacnem=dszetacnem_dt;
            xintegral=xintegral_dt;
            solnem=solnem_dt;
            s0=s0_dt;
            s0inv=s0inv_dt;
            sfun=sfun_dt;
            snewfun=snewfun_dt;
            solold=sol;
            
            lacurrent = s0inv(param_la*L0);
            
            if strcmp(FixedPar,'V')
                P = P1(1);
                V = V0;
            elseif strcmp(FixedPar,'P')
                P = -etap*P1(1);
                V = V + dt*dV;
            end
            L = Lnew_dt;
            eps1 = eps1new_dt;
            X0 = X0 + dt*dX0;
            
            tnew = t + dt;
            
            t = tnew;
            told = t-dt;
            i = i+1;
            
            %% save to file and plot (if applicable)
            if (mod(i,delt_data)==0)||(i==1)

                if  strcmp(Plotting,'On')
                    plotting(snewvec_dt, C1, C2, C, dsC1, dsC, sol, t, told, P1, P, P0, V, V0, R0, X0, U, zetanem, X, Z, Uinit, Xinit, Zinit, sinit, vidObj1);
                end

                % where to evaluate functions for saving:
                if param_la==1 %if hom. profile, evaluate at equator
                    seval=0.5*L;
                else
                    seval=0; %else, evaluate at SP
                end

                savetofile(X, Z, Psi, svec, snewvec_dt, seval, lacurrent, fileID, formatSpec, t, dt, P1, P, C1, C2, C, dsC1, dsC, xintegral, tcomp, L, dX0, V, X0, sol.stats.nmeshpoints, v, vs, vn, tss, U, Q, s0, zeta, zetac, zetanem, zetacnem, filename2, filename3, filename4, filename41, filename5, filename6, filename7, filename9, filename10, filename11, filename12, filename13, n, ControlPar);

                n = n + 1;
            end
       else
           told=t;
       end
       dt = dtnext;

    else
        break;
    end

    %% additional stopping criterion for successfull relaxation to steady state:
    if (i>1000)&&(max(abs(vn(svec)))<1e-4)
        if  strcmp(Plotting,'On')
            plotting(snewvec_dt, C1, C2, C, dsC1, dsC, sol, t, told, P1, P, P0, V, V0, R0, X0, U, zetanem, X, Z, Uinit, Xinit, Zinit, sinit, vidObj1);
        end
        
        % where to evaluate functions for saving
        if param_la==1 %if hom. profile, evaluate at equator
            seval=0.5*L;
        else
            seval=0; %else, evaluate at SP
        end
        
        savetofile(X, Z, Psi, svec, snewvec_dt, seval, lacurrent, fileID, formatSpec, t, dt, P1, P, C1, C2, C, dsC1, dsC, xintegral, tcomp, L, dX0, V, X0, sol.stats.nmeshpoints, v, vs, vn, tss, U, Q, s0, zeta, zetac, zetanem, zetacnem, filename2, filename3, filename4, filename41, filename5, filename6, filename7, filename9, filename10, filename11, filename12, filename13, n, ControlPar);
        
        disp('Steady state reached.')
        break;
    end
end
close(vidObj1);
close all;

end
