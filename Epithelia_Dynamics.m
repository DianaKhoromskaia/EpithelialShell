
% Dynamics for epithelial shells with active tensions and bending moments (isotropic and nematic)
% Code for the Publication "Active morphogenesis of patterned epithelial shells"

% authors: Diana Khoromskaia, Nicolas Cuny and Guillaume Salbreux
% date last edited: 14/02/2023

function Epithelia_Dynamics(l_open,dt, tmax, FixedPar, ProfileFile, dtmax, Adaptive, xi, tol, P0, thalf_P)
% arguments are:
% l_open - value of la at which is situated the free boundary (has to be between 0 and 1)
% dt - initial time step
% tmax - end time
% FixedPar - set to 'V' or 'P' for fixed or free volume, respectively
% ProfileFile - a string with the name of a file containing N lines of the
% form [ControlPar, ProfileType, param_const, param_la, param_fac, param_sigma, thalf],
% each line corresponds to one active region with:
% ControlPar : control parameter/ active effect: 'Tension', 'Bending', 'Nematic, 'BendingNematic'
% ProfileType : 'Sigmoidal', 'Gaussian' or 'Rectangle'
% param_const, param_la, param_fac, param_sigma, thalf : parameters defining the active profile in the region
% dtmax - maximal time step
% Adaptive - 'On' or 'Off'
% xi - normal friction coefficient; set to 0, not required. 
% tol - relative tolerance on adaptive time step

addpath('./functions_dynamics');

%% set physical parameters:
R0 = 1;
h = 1/3; % thickness of epithelium in units of R0
eta = 1; % shear viscosity
etab = 1; % bulk viscosity
etacb = (h^2)*etab; % bulk bending viscosity
etap = 1e-4; % pressure viscosity (for FixedPar='P')
kappa = 1; % bending modulus
%xi = 0; %normal friction coefficient
K = 110; % bulk elastic modulus
lc = 0.1*R0; % nematic length scale in units of R0

Plotting = 'Off'; % set 'On' or 'Off' to save results as movies

%% set numerical parameters:
%cut-off for nematic equations at SP and at NP of the surface in units of interval length:
eps1abs = 1e-4;% 1e-4;
eps2abs = 1e-4;%1e-2;
%relative and absolute error bounds for bvp-solver:
epsoderel = 1e-4;
epsodeabs = 1e-6;
%options for bvp-solvers: -- set 'stats' to 'off' if details of solver output not required
optode = bvpset('RelTol', epsoderel, 'AbsTol', epsodeabs, 'stats','on','Vectorized', 'on','FJacobian',@fjac,'BCJacobian',@bcjac,'NMax', 1e+5); %'FJacobian',@fjac,'BCJacobian',@bcjac, 
optode_nem = bvpset('RelTol', 1e-4, 'AbsTol', 1e-6, 'stats','on','Vectorized', 'on','FJacobian',@fjac_nematic,'BCJacobian',@bcjac_nematic,'NMax', 1e+5);

%tol = 1e-5; %relative tolerance for time step

% set time parameters
t = 0;
%dtmax = 1e-5;
%tmax = timesteps*dt;
%thalf = 0.01;
tsigma = 0.002;

%% sizes of different grids
nplot = 1e+3;       % for plotting and saving to file
npoints = 1e+6;     % for initialising functions on sphere, and for derivatives

%% initialise spherical shape on half of s-interval
z0 = 0.; %offset in z-direction
L = pi*R0;
L0 = l_open*pi*R0;
[C1, C2, C, dsC1, dsC, Psi, X, Z, X0, xintegral, svec1, U, dsU, Q, dsQ, dsw0, dswL, s0] = initialisesphere(L0, R0, z0, npoints);
eps1 = eps1abs*L;
eps2 = eps2abs*L;
V0 = (4/3)*pi*R0^3;
V = V0;

sgrid = linspace(0,L0,300);

% save for plotting:
sinit = svec1;
Xinit = X;
Zinit = Z;

%% initialise pressure and output directory, depending on FixedPar

if strcmp(FixedPar,'V')
    dir1 = strcat('Dynamics_FixedV_K=',num2str(K));
elseif strcmp(FixedPar,'P')
    if P0 ~= 0
        dir1 = strcat('Dynamics_FixedP_P0=',num2str(P0),'_thalfP=',num2str(thalf_P),'_K=',num2str(K));
    else
        dir1 = strcat('Dynamics_FixedP_K=',num2str(K));
    end
end
P = P0*(1-sigmoidal(0,thalf_P,tsigma));

%% initialise 'active' profiles:
% ..._Profile - set to 'Sigmoidal', 'Gaussian' or 'Rectangle'
% ..._const - constant offset
% ..._la - position of half-maximum or peak, resp.
% ..._sigma - controls sharpness of sigmoid or width of Gaussian/Rectangle
% ..._par - prefactor for profile
% ..._thalf - time at whcich the profile start appearing (the active profiles are modulated by a sigmoid in time)
% if 'Sigmoidal' and ..._la=1 then profile is flat: profile= ..._const

zetasrect=0.002; %width of the sigmoidal defining the 'Rectangle' active profile
        
[zeta, dszeta, zetac, dszetac, zetanem, dszetanem , zetacnem, dszetacnem, dir2, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_facs, zeta_sigmas, zeta_thalfs, N_regions, write9, write91, write92, write93] = initialiseprofiles(ProfileFile, sgrid, svec1, L, npoints, L0, s0, Q, t, tsigma, dir1, zetasrect);

%% output files:
mkdir(dir2);
copyfile(ProfileFile,strcat(dir2,'/profiles_file.dat'));
cd(dir2);

% plot profiles to check:
figure(1)
subplot(2,2,1)
plot(svec1, zeta(svec1), svec1, dszeta(svec1));
legend('\zeta','\partial_s\zeta');
subplot(2,2,2)
plot(svec1, zetac(svec1), svec1, dszetac(svec1));
legend('\zeta_c','\partial_s\zeta_c');
subplot(2,2,3)
plot(svec1, zetanem(svec1), svec1, dszetanem(svec1));
legend('\zeta_{n}','\partial_s\zeta_{n}');
subplot(2,2,4)
plot(svec1, zetacnem(svec1), svec1, dszetacnem(svec1));
legend('\zeta_{cn}','\partial_s\zeta_{cn}');

saveas(figure(1),'activeprofiles.fig');
%saveas(figure(1),'activeprofiles.png');

seval=0;

filename1 = 'observables.dat';
fileID = fopen(filename1,'w');
fprintf(fileID, '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n', 't', 'dt', 'v_n(seval)', 'P', 'X(seval)', 'C2(seval)','pole-pole', 'area', 'tcomp [sec]', 'L', 'V','X0', 'nmesh','intersect', 'seval', 'f_c','I(L)');
formatSpec = '%e \t %e \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %i \t %i \t %2.1f \t %6.5f \t %9.8f \n';
fprintf(fileID, formatSpec, t, dt, 0., P0, X(seval), C2(seval), Z(L0)-Z(0.), 2*pi*xintegral, 0., L0, V0, X0, 0., 0., seval, 0, 0);

svecuni = linspace(0., L, nplot);%uniform grid for saving

% save functions of initial shape:
filename2 = 'curvatures.dat';
dlmwrite(filename2, svecuni, 'precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C1(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C2(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, C(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
dlmwrite(filename2, dsC1(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');
%dlmwrite(filename2, dsC2(svecuni), '-append','precision', '%6.5f' ,'delimiter', '\t');
dlmwrite(filename2, dsC(svecuni), '-append','precision', '%10.9f' ,'delimiter', '\t');

filename3 = 'x.dat';
dlmwrite(filename3, X(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');
filename4 = 'z.dat';
dlmwrite(filename4, Z(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');
filename41 = 'psi.dat';
dlmwrite(filename41, Psi(svecuni), 'precision', '%10.9f' ,'delimiter', '\t');

%save output of bvp-solver:
filename5 = 'bvpsolution.dat';
filename6 = 'vs.dat';
filename7 = 'vn.dat';
filename9 = 'profile.dat';
filename91 = 'profilec.dat';
filename92 = 'profilenem.dat';
filename93 = 'profilecnem.dat';
filename10 = 'u.dat';
filename11 = 'q.dat';
filename12 = 'tss.dat';
filename13 = 's0.dat';
filename14 = 'svecnew.dat';

%time step for saving into data file:
delt_data = 10; %every 10-th step
n = 1.; %number saved time steps
nsteps = 0;
i = 0; %number time steps
error = 0;

%external force magnitude along z-direction: DO NOT USE THIS FOR NOW
fextmag = 0;
fext = griddedInterpolant(svec1, fextmag*svec1, 'spline');

%% save parameters into file

T1 = table({'FixedPar';'zetasrect';'R0';'h';'eta';'etab';'etacb';'etap';'kappa';'xi';'K';'lc'},...
        [{FixedPar};zetasrect;R0;h;eta;etab;etacb;etap;kappa;xi;K;lc],...
        {'fixed volume (V) or free volume (P)';'width of the transition for Rectangle profiles';'radius';'thickness';'shear viscosity';'bulk viscosity';'bending viscosity';'volume viscosity';'bending modulus';'friction';'elastic modulus';'nematic length, in units of R0'},'VariableNames',{'Parameter','Value','Description'});
writetable(T1,'parameters_physical.csv');

T2 = table({'Adaptive';'dt';'dtmax';'tmax';'tol';'epsoderel';'epsodeabs';'eps1abs';'nplot';'npoints';'delt_data'},...
    [{Adaptive};dt;dtmax;tmax;tol;epsoderel;epsodeabs;eps1abs;nplot;npoints;delt_data],...
    {'adaptive stepping';'initial time step';'max time step';'max total time';'rel. error tolerance for adaptive stepping';'rel. error for sover';'abs. error for solver';'cutoff at poles, used for nem. solver only';'size of uniform grid for saving initial shape';'size of uniform grid for profiles and gradients';'every n-th time step saved'},'VariableNames',{'Parameter','Value','Description'});
writetable(T2,'parameters_numerical.csv');

while t < tmax
    tnew = t
    
    %% solve force balance at time t (with step dt and dt/2 for arc length):
    if t==0
        tic
        [svec, v, P1, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] ...
            = forcebalance(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, dt, Adaptive, FixedPar, P0, thalf_P, tsigma, sgrid);
        tcomp = toc;
    else
        tic
        [svec, v, P1, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] ...
            = forcebalance(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, dt, Adaptive, FixedPar, P0, thalf_P, tsigma, solold, sfun, snewfun);
        tcomp = toc;
    end
    
    if SolFound
        %% evolve system with time step dt:
        if t==0
            [C1_dt, dsC1_dt, C2_dt, C_dt, dsC_dt, X_dt, Psi_dt, Z_dt, U_dt, dsU_dt, zeta_dt, dszeta_dt, zetac_dt, dszetac_dt, zetanem_dt, dszetanem_dt, zetacnem_dt, dszetacnem_dt, xintegral_dt, solnem_dt, s0_dt, s0inv_dt, Q_dt] = ...
                evolvefunctions_t0(svec, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac, dszetac, zetacnem, dszetacnem, ...
                vs, dsvs, vkk, vn, dsvn, mss, tns, ...
                U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dt, sfun_dt, t, dt, tsigma, L, Lnew_dt, L0, eps1, eps2, etacb, npoints, dsw0, dswL, lc, Q, dsQ, s0, optode_nem, N_regions, zetasrect);
        else
            [C1_dt, dsC1_dt, C2_dt, C_dt, dsC_dt, X_dt, Psi_dt, Z_dt, U_dt, dsU_dt, zeta_dt, dszeta_dt, zetac_dt, dszetac_dt, zetanem_dt, dszetanem_dt, zetacnem_dt, dszetacnem_dt, xintegral_dt, solnem_dt, s0_dt, s0inv_dt, Q_dt] = ...
                evolvefunctions(svec, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac, dszetac, zetacnem, dszetacnem, ...
                vs, dsvs, vkk, vn, dsvn, mss, tns, ...
                U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dt, sfun_dt, t, dt, tsigma, L, Lnew_dt, L0, eps1, eps2, etacb, npoints, solnem, lc, s0, optode_nem, N_regions, zetasrect);
        end
        
        if strcmp(Adaptive,'On')
            %% evolve system with time step dt/2:
            if t==0
                [C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, X_dthalf, Psi_dthalf, Z_dthalf, U_dthalf, dsU_dthalf, zeta_dthalf, dszeta_dthalf, zetac_dthalf, dszetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, dszetacnem_dthalf, xintegral_dthalf, solnem_dthalf, s0_dthalf, s0inv_dthalf, Q_dthalf] = ...
                    evolvefunctions_t0(svec, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac, dszetac, zetacnem, dszetacnem, ...
                    vs, dsvs, vkk, vn, dsvn, mss, tns, ...
                    U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dthalf, sfun_dthalf, t, dt/2, tsigma, L, Lnew_dthalf, L0, eps1, eps2, etacb, npoints, dsw0, dswL, lc, Q, dsQ, s0, optode_nem, N_regions, zetasrect);
            else
                [C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, X_dthalf, Psi_dthalf, Z_dthalf, U_dthalf, dsU_dthalf, zeta_dthalf, dszeta_dthalf, zetac_dthalf, dszetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, dszetacnem_dthalf, xintegral_dthalf, solnem_dthalf, s0_dthalf, s0inv_dthalf, Q_dthalf] = ...
                    evolvefunctions(svec, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac, dszetac, zetacnem, dszetacnem, ...
                    vs, dsvs, vkk, vn, dsvn, mss, tns, ...
                    U, C1, dsC1, C2, C, dsC, kappa, X, Psi, Z, snewfun_dthalf, sfun_dthalf, t, dt/2, tsigma, L, Lnew_dthalf, L0, eps1, eps2,etacb, npoints, solnem, lc, s0, optode_nem, N_regions, zetasrect);
            end
            
            L_dthalf = Lnew_dthalf;
            eps1_dthalf = eps1new_dthalf;
            X0_dthalf = X0 + 0.5*dt*dX0;
            %% solve force balance at time t+dt/2 with step dt/2:
            [svec_dt2half, v_dt2half, P1_dt2half, sol_dt2half, vs_dt2half, dsvs_dt2half, vkk_dt2half, tss_dt2half, vn_dt2half, dsvn_dt2half, mss_dt2half, tns_dt2half, ds2vn_dt2half, dV_dt2half, dX0_dt2half, Lnew_dt2half, eps1new_dt2half, sfun_dt2half, snewfun_dt2half, snewvec_dt2half, SolFound_dt2half] ...
                = forcebalance_dthalf(U_dthalf, dsU_dthalf, C1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, Psi_dthalf, X_dthalf, Z_dthalf, X0_dthalf, L_dthalf, L0, zeta_dthalf, dszeta_dthalf, zetac_dthalf, zetanem_dthalf, dszetanem_dthalf, zetacnem_dthalf, eta, etab, etacb, etap, eps1abs, xintegral_dthalf, fext, kappa, K, xi, optode, t+dt/2, dt/2, FixedPar, P0, thalf_P, tsigma, sol, sfun_dthalf, snewfun_dthalf);
            
            %% evolve system with time step dt/2:
            [C1_dt2half, dsC1_dt2half, C2_dt2half, C_dt2half, dsC_dt2half, X_dt2half, Psi_dt2half, Z_dt2half, U_dt2half, dsU_dt2half, zeta_dt2half, dszeta_dt2half, zetac_dt2half, dszetac_dt2half, zetanem_dt2half, dszetanem_dt2half, zetacnem_dt2half, dszetacnem_dt2half, xintegral_dt2half, solnem_dt2half, s0_dt2half, s0inv_dt2half, Q_dt2half] = ...
                evolvefunctions(svec_dt2half, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_sigmas, zeta_facs, zeta_thalfs, zetac_dthalf, dszetac_dthalf, zetacnem_dthalf, dszetacnem_dthalf, ...
                vs_dt2half, dsvs_dt2half, vkk_dt2half, vn_dt2half, dsvn_dt2half, mss_dt2half, tns_dt2half, ...
                U_dthalf, C1_dthalf, dsC1_dthalf, C2_dthalf, C_dthalf, dsC_dthalf, kappa, X_dthalf, Psi_dthalf, Z_dthalf, snewfun_dt2half, sfun_dt2half, t+dt/2, dt/2, tsigma, L_dthalf, Lnew_dt2half, L0, eps1, eps2, etacb, npoints, solnem_dthalf, lc, s0_dthalf, optode_nem, N_regions, zetasrect);
            
            %% compare results for adaptive time step after dt and 2*(dt/2):
            scompare = linspace(0,L,1000);
            scompare_dt = snewfun_dt(scompare);
            scompare_dt2half = snewfun_dt2half(snewfun_dthalf(scompare));
            
            % exclude values from relative error calculation which are too close to zero:
            cutoff = 1e-3;
            Z_dt_vec = Z_dt(scompare_dt);
            Z_dt2half_vec = Z_dt2half(scompare_dt2half);
            dsC_dt_vec = dsC_dt(scompare_dt);
            dsC_dt2half_vec = dsC_dt2half(scompare_dt2half);
            
            format long 
            
            errorvec = [max(abs(X_dt(scompare_dt(2:end-1))-X_dt2half(scompare_dt2half(2:end-1)))./abs(X_dt(scompare_dt(2:end-1))));...
                max(abs((Z_dt_vec(abs(Z_dt_vec)>cutoff)-Z_dt2half_vec(abs(Z_dt_vec)>cutoff))./Z_dt_vec(abs(Z_dt_vec)>cutoff)));...
                max(abs((dsC_dt_vec(abs(dsC_dt_vec)>cutoff)-dsC_dt2half_vec(abs(dsC_dt_vec)>cutoff))./dsC_dt_vec(abs(dsC_dt_vec)>cutoff)))]
            
            error = max(errorvec);
            dtnext = min(0.9*dt*min([max([tol/error 0.5]) 1.5]), dtmax)
            if dtnext < 1e-10
                dtnext = 1e-10
            end
        else 
            dtnext = dt; 
        end

       if (error<tol)||(dtnext==1e-10)||strcmp(Adaptive,'Off') %% accept step only if tol satisfied, accept first step always; instead: first 1000 steps; ||(i<1000)
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
                        
            if strcmp(FixedPar,'V')
                P = P1(1);
                V = V0;
            elseif strcmp(FixedPar,'P')
                P = P0*(1-sigmoidal(t,thalf_P,tsigma))-etap*P1(1);
                V = V + dt*dV;
            end
            L = Lnew_dt;
            eps1 = eps1new_dt;
            X0 = X0 + dt*dX0;
            
            tnew = t + dt;
            nsteps = nsteps + 1;
            
            t = tnew;
            told = t-dt;
            i = i+1;
            
            %% save to file and plot (if applicable)
            if (mod(i,delt_data)==0)||(i==1)  %t >= n*delt_data
                
                if  strcmp(Plotting,'On')
                    plotting(snewvec_dt, sol, t, told, P, P0, V, V0, R0, X0, X, Z, Xinit, Zinit, sinit);
                end
                
                %% where to evaluate
                                
                savetofile(X, Z, Psi, svec, snewvec_dt, seval, fileID, formatSpec, t, dt, P1, P, C1, C2, C, dsC1, dsC, xintegral, tcomp, L, dX0, V, X0, sol.stats.nmeshpoints, v, vs, vn, tss, U, Q, s0, zeta, zetac, zetanem, zetacnem, filename2, filename3, filename4, filename41, filename5, filename6, filename7, filename9, filename91, filename92, filename93, write9, write91, write92, write93, filename10, filename11, filename12, filename13, filename14, n);
                
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
            plotting(snewvec_dt, sol, t, told, P, P0, V, V0, R0, X0, X, Z, Xinit, Zinit, sinit);
        end
        
        %% where to evaluate
                
        savetofile(X, Z, Psi, svec, snewvec_dt, seval, fileID, formatSpec, t, dt, P1, P, C1, C2, C, dsC1, dsC, xintegral, tcomp, L, dX0, V, X0, sol.stats.nmeshpoints, v, vs, vn, tss, U, Q, s0, zeta, zetac, zetanem, zetacnem, filename2, filename3, filename4, filename41, filename5, filename6, filename7, filename9, filename91, filename92, filename93, write9, write91, write92, write93, filename10, filename11, filename12, filename13, filename14, n);
        
        disp('Steady state reached.')
        break;
    end
end
close all;
cd ..
end
