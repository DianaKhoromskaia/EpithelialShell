% This code calculates the steady state solution branches for epithelial shells with active tensions and bending moments
% (isotropic and nematic), using step-profiles or smooth profiles (Sigmoidal or Gaussian). 
% Code for the publication "Active morphogenesis of patterned epithelial
% shells".

% Authors: Diana Khoromskaia and Guillaume Salbreux
% Date last edited: 27/10/2022

function Epithelia_SteadyState(la, K, reltol, ControlPar, FixedPar, BC, Sign, Profile, Running)

% Computational parameters (see below, read in from file):
% la = size of red region; special case for la=1
% K = area elastic modulus of the surface
% nsteps = total number of iterations, e.g. 40'000
% NonParseq = number of steps (out of nsteps), which are done via explicit stepping 
% lstep1 = step size in control parameter (explicit stepping), e.g. 0.1
% lstep2 = step size in control parameter (parametric curve stepping), e.g. 0.1*lstep1
% stepsave1 = saving every ... iteration to file in explicit stepping
% stepssave2 = saving every ... iteration to file in parametric stepping
% allsave1 = saving every ... iteration current shape and functions to file
% (explicit stepping) as input for dynamics simulation
% allsave2 = saving every ... iteration current shape and functions to file
% (parametric stepping) as input for dynamics simulation
% 
% Input parameters:
% ControlPar = the physical parameter varied in current simulation (control
% parameter): 'Bending', 'Nematic', 'Tension', 'Force', 'BendingNematic'
% Running = specifying whether running on 'Cluster' or on 'Laptop'
% FixedPar = 'V' or 'P' (volume or pressure fixed)
% BC = 'Sphere', 'DiskFree', 'DiskFixed' (closed sphere, disk with free edge or with edge clamped at initial position)
% Sign = 'pos', 'neg' referring to sign of controlpar sequence

addpath('./functions_steadystate');
warning('off','all'); % this turns off all warnings, in particular the one of bvp4c.m in the case of split interval and discontinuities across that boundary

%% read in remaining model parameters:
% R0=1;
physparam = table2array(readtable('physical_par.csv'));
C0 = physparam(1);
aM = physparam(2); %jump in curvature C_k^k
zeta = physparam(3)*K;
zetared = physparam(4)*K;
zetanem = physparam(5);
zetanemred = physparam(6);
fconst = physparam(7); %f = fconst*\vec{e}_s   traction force along the boundary
fconstred = physparam(8);
zetacnem = physparam(9);
zetacnemred = physparam(10);
L0 = 1; %dim-less, in units of pi*R0
lc = physparam(11); %0.1; %characteristic length of the nematic in units of R0

%% read in numerical parameters:
filenameNumPar = strcat('numerical_par_',ControlPar,'_fixed',FixedPar,'_K=',num2str(K),'_',Sign,'.csv');
numparam = table2array(readtable(filenameNumPar));
nsteps = numparam(1);
NonParseq = numparam(2); 
lstep1 = numparam(3);
lstep2 = numparam(4);
stepsave1 = numparam(5);
stepsave2 = numparam(6);
allsave1 = numparam(7);
allsave2 = numparam(8);
gradientswitchlower = numparam(9);
gradientswitchupper = numparam(10);
profilewidth = 200;%numparam(11);

%% solver parameters:
%reltol = 1e-4;
abstol = 1e-6;
Stats = 'on';
ninit = 500; %#points in initial grid for first iteration
optode = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized', 'on','FJacobian',@fjac, 'BCJacobian', @bcjac,'NMax', 1e+4); 
optode_param = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized','on','FJacobian',@fjac_param,'BCJacobian', @bcjac_param,'NMax', 1e+7); 

optode_full = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized', 'on','FJacobian',@fjac_full,'BCJacobian', @bcjac_full,'NMax', 1e+7); 
optode_full_param = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized', 'on','FJacobian',@fjac_full_param,'BCJacobian', @bcjac_full_param,'NMax', 1e+7); 

optode_smooth = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized', 'on','FJacobian', @fjac_smooth, 'BCJacobian', @bcjac_smooth, 'NMax', 1e+7); %,
optode_smooth_param = bvpset('RelTol', reltol, 'AbsTol', abstol, 'stats',Stats, 'Vectorized', 'on','FJacobian',@fjac_smooth_param,'BCJacobian', @bcjac_smooth_param,'NMax', 1e+7);

%% create folder for given control parameter, containing date and time:
switch Running
    case 'Laptop'
        dirnameControlPar = strcat(ControlPar,'_fixed',FixedPar,'_',Profile,'_',BC,'_',Sign,'_',datestr(now,30)); % date and time format: 'yyyymmddTHHMMSS'
    case 'Cluster'
        dirnameControlPar = strcat(ControlPar,'_fixed',FixedPar,'_',Profile,'_',BC,'_',Sign,'_',datestr(now,29)); % date format: 'yyyy-mm-dd'
end

mkdir(dirnameControlPar);
cd(dirnameControlPar) % go to control parameter folder

%% set plotting options for saveobservables.m : 
plotopt = 'movie';
if strcmp(Running, 'Cluster')
    plotopt = 'off';
end

%% initialise remaining simulation parameters, e.g. parameter array = (CPar, P, L1) for parameteric curve stepping
controlparstep = lstep1; %value for i=2
controlparinit = 0; %value for i=1
controlpar = controlparinit;
P0 = 2*zeta;
FixedParOffset = 1; %for first step along control parameter, since with P=0 the solution is usually harder to find
V0 = 1;

% set initial grid for the solver:
Qinit = dlmread('Qinit.dat','\t'); %solution for nematic on sphere
if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
    switch FixedPar 
        case 'V'
            Param = [controlparinit P0 L0;...
                    0. 0. 0.];
            paraminit = [P0 0. Qinit(1,1) Qinit(1,1) 1]; % = (P, fc, dsw_0, dsw_L, L) 
            fixedpar = V0; %numerical value of fixed parameter (volume)
        case 'P'
            Param = [controlparinit V0 L0;...
                    0. 0. 0.];
            paraminit = [V0 0. Qinit(1,1) Qinit(1,1) 1]; % = (V, fc, dsw_0, dsw_L, L) 
            fixedpar = P0; %numerical value of fixed parameter (pressure)
        otherwise
            disp('FixedPar value invalid. Choose V(olume) or P(ressure).')
    end
    sgridinit = linspace(0,1,ninit);
    solinit = bvpinit(sgridinit, @(s) yguessfun_init_fullint(s,L0,zeta), paraminit); 
elseif strcmp(Profile,'Step')
    switch FixedPar 
        case 'V'
            Param = [controlparinit P0 la*L0;...
                    0. 0. 0.];
            paraminit = [P0 0. Qinit(1,1) Qinit(1,1) la 1-la]; % = (P, fc, dsw_0, dsw_L, L1, L2) 
            fixedpar = V0;
        case 'P'
            Param = [controlparinit V0 la*L0;...
                    0. 0. 0.];
            paraminit = [V0 0. Qinit(1,1) Qinit(1,1) la 1-la]; % = (P, fc, dsw_0, dsw_L, L1, L2) 
            fixedpar = P0;
        otherwise
            disp('FixedPar value invalid. Choose V(olume) or P(ressure).')
    end
    sgridinit = [linspace(0,1,ninit) linspace(1,2,ninit)];
    solinit = bvpinit(sgridinit, @(s,k) yguessfun_init(s,k,la,zeta), paraminit);
end
paramvec = zeros(4,length(paraminit));
paramextrap = zeros(length(paraminit),1);

% sol-evaluation and saving parameters:
nsample = 1000; % points in uniform grid for function evaluation
epsstep = 1e-6; % for evaluating at the interval boundary

%% create folder and data files:
dirname = strcat('Results_la=',num2str(la),'_K=',num2str(K));
mkdir(dirname);
cd(dirname)

filename1 = 'observables.dat';
fileID = fopen(filename1, 'w');
fprintf(fileID, 'ControlPar \t P \t L \t L1 \t L2 \t V \t A \t fc \t C(0) \t C(L) \t max(|c2|) \t s_c2max \t max(|c1|) \t s_c1max \t z(L) \t tcomp \t mesh \t intersect \t dsw0 \t dswL \t index \t seq \t Fhalf \t x(half)');
formatSpec = '\n %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f \t %6.5f  \t %6.5f \t %6.5f \t %6.5f  \t %6.5f \t %6.5f \t %d \t %d \t %6.5f \t %6.5f \t %d \t %d \t %6.5f \t %6.5f';

% write all parameters in one file:
dlmwrite('parameters.dat', [la K reltol abstol nsteps NonParseq lstep1 lstep2 stepsave1 stepsave2 allsave1 allsave2 gradientswitchlower gradientswitchupper ninit nsample FixedParOffset C0 aM zeta zetared zetanem zetanemred fconst fconstred zetacnem zetacnemred L0 lc], 'precision', '%6.5f' ,'delimiter', '\t');

% create video file if necessary
if strcmp(plotopt,'movie')
    vidObj1 = VideoWriter('shapes.avi');%, 'MPEG-4');%,'MPEG-4');
    vidObj1.Quality = 25;
    vidObj1.FrameRate = 25;
    open(vidObj1);
    movie01 = vidObj1;
else 
    movie01 = '';
end

%% do several steps in ControlPar-sequence, then follow the parametric curve for remaining of nsteps: 
for i=1:(NonParseq-1)

    if strcmp(Running, 'Laptop') 
        %i 
    end
    if i>1 %saving for case when new solution can't be found 
        solold = sol;
        controlparold = controlpar;
    end
    SolFound = 'True';
    if i==1 % first step: recover the spherical solution
        switch ControlPar
            case 'Bending'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    C0 = controlpar;
                else
                    aM = controlpar;
                end
            case 'Nematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetanem = controlpar;
                else
                    zetanemred = controlpar;
                end
            case 'Tension'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zeta = controlpar;
                else
                    zetared = controlpar;
                end
            case 'Force'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    fconst = controlpar;
                else
                    fconstred = controlpar;
                end
            case 'BendingNematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetacnem = controlpar;
                else
                    zetacnemred = controlpar;
                end                   
            otherwise
                disp('ControlPar value invalid.')
        end
        
        tic
        try
            if la==1
                sol = bvp4c( @ode_full, @bc_full, solinit, optode_full, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&strcmp(Profile,'Step')
                sol = bvp4c( @ode, @bc, solinit, optode, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&~strcmp(Profile,'Step')
                sol = bvp4c( @ode_smooth, @bc_smooth, solinit, optode_smooth, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar, Profile, profilewidth);
            end
            paramvec(1,:) = sol.parameters;
        catch ME
            ME
            disp(strcat('Could not find solution at i=',num2str(i),' and controlpar = ',num2str(controlpar),'. Exiting.'));
            SolFound = 'False';
            %break; % go out of main loop
        end
        tcomp = toc;
        
    elseif i==2 % second step: first deformation
        controlpar = controlparstep;
        switch ControlPar
            case 'Bending'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    C0 = controlpar;
                else
                    aM = controlpar;
                end
            case 'Nematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetanem = controlpar;
                else
                    zetanemred = controlpar;
                end
            case 'Tension'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zeta = controlpar;
                else
                    zetared = controlpar;
                end
            case 'Force'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    fconst = controlpar;
                else
                    fconstred = controlpar;
                end
            case 'BendingNematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetacnem = controlpar;
                else
                    zetacnemred = controlpar;
                end
        end
        
        tic
        try
            if la==1
                paroffset = zeros(size(sol.parameters));
                paroffset(1) = FixedParOffset;
                solinit = bvpinit(sgridinit, @(s) yguessfun_init_fullint(s,L0,zeta), sol.parameters + paroffset);
                sol = bvp4c( @ode_full, @bc_full, solinit, optode_full, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&strcmp(Profile,'Step')
                paroffset = zeros(size(sol.parameters));
                paroffset(1) = FixedParOffset;
                solinit = bvpinit(sgridinit, @(s,k) yguessfun_init(s,k,la,zeta), sol.parameters + paroffset);
                sol = bvp4c( @ode, @bc, solinit, optode, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&~strcmp(Profile,'Step')
                paroffset = zeros(size(sol.parameters));
                paroffset(1) = FixedParOffset;
                solinit = bvpinit(sgridinit, @(s) yguessfun_init_fullint(s,L0,zeta), sol.parameters + paroffset);
                sol = bvp4c( @ode_smooth, @bc_smooth, solinit, optode_smooth, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar, Profile, profilewidth);
            end
            Param(2,:) = [controlpar sol.parameters(1) sol.parameters(5)];
            paramvec(2,:) = sol.parameters;
        catch ME
            ME
            disp(strcat('Could not find solution at i=',num2str(i),' and controlpar = ',num2str(controlpar),'. Exiting.'))
            saveobservables(solold, la, K, controlparold, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 0, movie01, ControlPar, Profile, profilewidth);
            SolFound = 'False';
            break;
        end
        tcomp = toc;
        
    elseif (i>2)&&(i<NonParseq) % follow solution branch with non-parametric sequence of steps along control parameter
        controlpar = controlpar + lstep1;
        paramextrap = sol.parameters; 
        if i>=5 % extrapolate parameter values for the solver guess
            paramextrap(1) = interp1([1 2 3 4], paramvec(:,1), 5, 'cubic','extrap');
            paramextrap(5) = interp1([1 2 3 4], paramvec(:,5), 5, 'cubic','extrap');
            if (la<1)&&strcmp(Profile,'Step')
                paramextrap(6) = interp1([1 2 3 4], paramvec(:,6), 5, 'cubic','extrap');   
            end
        end
        if (la<1)&&strcmp(Profile,'Step')
            solinit = bvpinit(sol.x, @(s,k) yguessfun_extrap(s,k,sol), paramextrap);
        elseif (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
            solinit = bvpinit(sol.x, @(s) yguessfun_extrap(s,sol), paramextrap);
        end
        switch ControlPar
            case 'Bending'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    C0 = controlpar;
                else
                    aM = controlpar;
                end
            case 'Nematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetanem = controlpar;
                else
                    zetanemred = controlpar;
                end
            case 'Tension'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zeta = controlpar;
                else
                    zetared = controlpar;
                end
            case 'Force'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    fconst = controlpar;
                else
                    fconstred = controlpar;
                end
            case 'BendingNematic'
                if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                    zetacnem = controlpar;
                else
                    zetacnemred = controlpar;
                end        
        end
        
        tic
        try
            if la==1  
                sol = bvp4c( @ode_full, @bc_full, solinit, optode_full, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&strcmp(Profile,'Step')
                sol = bvp4c( @ode, @bc, solinit, optode, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, FixedPar, fixedpar);
            elseif (la<1)&&~strcmp(Profile,'Step')
                sol = bvp4c( @ode_smooth, @bc_smooth, solinit, optode_smooth, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, FixedPar, fixedpar, Profile, profilewidth);   
            end            
            Param(1,:) = Param(2,:);
            Param(2,:) = [controlpar sol.parameters(1) sol.parameters(5)];
            if i==3
                paramvec(3,:) = sol.parameters;
            elseif i==4
                paramvec(4,:) = sol.parameters;
            else
                paramvec(1:3,:) = paramvec(2:4,:);
                paramvec(4,:) = sol.parameters;
            end
        catch ME
            ME
            disp(strcat('Could not find solution at i=',num2str(i),' and controlpar = ',num2str(controlpar),' in explicit stepping. Starting the parametric stepping.'))
            saveobservables(solold, la, K, controlparold, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 0, movie01, ControlPar, Profile, profilewidth);
            SolFound = 'False';
            ifinal = i; %go to end of explicit stepping sequence
            break;
        end
        tcomp = toc;
    end
    
    if strcmp(SolFound,'True')
        %% Calculate rmin:
        [c1max, ind_c1max] = max(abs(sin(sol.y(4,2:end-1))./sol.y(5,2:end-1)));
        rmin = 1/c1max;
        
        %% Test shape for intersection
        if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
            s_sample = linspace(0., 1, nsample);
        else
            s_sample_left = linspace(0., 1-epsstep, nsample*la);
            s_sample_right = linspace(1+epsstep, 2, nsample*(1-la));
            s_sample = [s_sample_left s_sample_right];
        end
        [sol_sample, solp_sample] = deval(sol, s_sample);
        Xvec = sol_sample(5,:);
        Zvec = sol_sample(6,:);
        [X0, Y0, I, J] = intersections(Xvec, Zvec);
        intersec = int8(~isempty(X0));
    else 
        intersec=0;
    end
    
    %% saving shape and functions to file for dynamics: on true length interval
    if ((i==1)||((i<NonParseq)&&(mod(i,allsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,allsave2)==0)))&&(strcmp(SolFound,'True'))
        savestates(sol, la, controlpar, ControlPar, NonParseq, epsstep, allsave1, allsave2, i, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, FixedPar, Profile, profilewidth);
    end

    %% saving observables and functions to file: on unit interval
    if ((i==1)||((i<NonParseq)&&(mod(i,stepsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,stepsave2)==0))||((rmin<0.01)&&strcmp(ControlPar,'Bending')))&&(strcmp(SolFound,'True'))  %||(rmin<0.1)
        saveobservables_extrap(sol, la, K, controlpar, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 0, movie01, paramextrap, ControlPar, Profile,profilewidth);
    end
    
    %% stop simulation if shape is self-intersecting
    if intersec==1
        disp(strcat('Shape has self-intersected. Stopping simulation at controlpar = ',num2str(controlpar)));
        saveobservables(sol, la, K, controlpar, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 0, movie01, ControlPar, Profile,profilewidth);
        break;
    end
    
    %% switch to parameteric, if constriction approaching
%     if ~strcmp(ControlPar,'BendingNematic')
%         if rmin<0.1
%             disp(strcat('Constriction approaching, rmin<0.1. Switching to parameteric sequence at controlpar = ',num2str(controlpar)));
%             break;
%         end
%     end
    
    %% criterion to switch to parametric sequence:
    tangentnormalised = norm((Param(2,2)-Param(1,2))*Param(2,1)/((Param(2,1)-Param(1,1))*Param(2,2))); %(dP/dcontrolpar)/P   *sol.parameters(1)
    if ((tangentnormalised<gradientswitchlower)||(tangentnormalised>gradientswitchupper))&&strcmp(ControlPar,'Nematic')
        disp(strcat('Relative gradient in P is too large/big. Switching to parameteric sequence at controlpar=',num2str(controlpar)));
        break;
    end

end

disp('Starting parametric stepping.');

if intersec~=1
for i=(NonParseq):nsteps  
        solold = sol;
        controlparold = controlpar;
        SolFound = 'True';
%         if strcmp(Running, 'Laptop')
%             %i
%         end
        % calculate tangent vector to parametric curve at point i-1:
        tangentvec = (Param(2,:)-Param(1,:))/norm((Param(2,:)-Param(1,:)));
        guess = Param(2,:) + lstep2*tangentvec;
        paramextrap = sol.parameters; 
        if i>=5
            paramextrap(1) = interp1([1 2 3 4], paramvec(:,1), 5, 'cubic','extrap');
            paramextrap(5) = interp1([1 2 3 4], paramvec(:,5), 5, 'cubic','extrap');
            if (la<1)&&strcmp(Profile,'Step')
                paramextrap(6) = interp1([1 2 3 4], paramvec(:,6), 5, 'cubic','extrap');   
            end
        end
        if (la<1)&&strcmp(Profile,'Step')
            solinit = bvpinit(sol.x, @(s,k) yguessfun_extrap(s,k,sol), paramextrap);
        elseif (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
            solinit = bvpinit(sol.x, @(s) yguessfun_extrap(s,sol), paramextrap);
        end
    
        tic
        try
            switch ControlPar
                case 'Bending'
                    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                        C0 = controlpar;
                    else
                        aM = controlpar;
                    end
                case 'Nematic'
                    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                        zetanem = controlpar;
                    else
                        zetanemred = controlpar;
                    end
                case 'Tension'
                    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                        zeta = controlpar;
                    else
                        zetared = controlpar;
                    end
                case 'Force'
                    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                        fconst = controlpar;
                    else
                        fconstred = controlpar;
                    end
                case 'BendingNematic'
                    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
                        zetacnem = controlpar;
                    else
                        zetacnemred = controlpar;
                    end 
            end
            if la==1
                sol = bvp4c( @ode_full_param, @bc_full_param, solinit, optode_full_param, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar);
            elseif (la<1)&&strcmp(Profile,'Step')
                sol = bvp4c( @ode_param, @bc_param, solinit, optode_param, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar);
            elseif (la<1)&&~strcmp(Profile,'Step')
                sol = bvp4c(@ode_smooth_param, @bc_smooth_param, solinit, optode_smooth_param, K, zeta, zetanem, zetacnem, C0, fconst, lc, la, guess, tangentvec, ControlPar, FixedPar, fixedpar, Profile, profilewidth); 
            end           
            controlpar = guess(1)-dot([sol.parameters(1) sol.parameters(5)]-guess(2:end),tangentvec(2:end))/tangentvec(1);
            Param(1,:) = Param(2,:);
            Param(2,:) = [controlpar sol.parameters(1) sol.parameters(5)];
            paramvec(1:3,:) = paramvec(2:4,:);
            paramvec(4,:) = sol.parameters;
        catch ME
            ME
            disp(strcat('Could not find solution at i=',num2str(i),' and controlpar = ',num2str(controlpar),' in parametric stepping. Exiting.'))
            saveobservables(solold, la, K, controlparold, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 0, movie01, ControlPar, Profile, profilewidth);
            SolFound = 'False';
            break; %go out of main loop.
        end
        tcomp = toc;
    
    %% Test shape for intersection
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        s_sample = linspace(0., 1, nsample);
    else
        s_sample_left = linspace(0., 1-epsstep, nsample*la);
        s_sample_right = linspace(1+epsstep, 2, nsample*(1-la));
        s_sample = [s_sample_left s_sample_right];
    end
    [sol_sample, solp_sample] = deval(sol, s_sample);
    Xvec = sol_sample(5,:);
    Zvec = sol_sample(6,:);
    [X0, Y0, I, J] = intersections(Xvec, Zvec);
    intersec = int8(~isempty(X0));        
        
    %% Calculate rmin:
    [c1max, ind_c1max] = max(abs(sin(sol.y(4,2:end-1))./sol.y(5,2:end-1)));
    rmin = 1/c1max;
    
    %% saving shape and functions to file for dynamics: on true length interval
    if ((i==1)||((i<NonParseq)&&(mod(i,allsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,allsave2)==0)))&&(strcmp(SolFound,'True'))
        savestates(sol, la, controlpar, ControlPar, NonParseq, epsstep, allsave1, allsave2, i, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, FixedPar, Profile, profilewidth);
    end
        
        %% saving observables and functions to file: on unit interval and on true (rescaled interval)
    if ((i==1)||((i<NonParseq)&&(mod(i,stepsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,stepsave2)==0)))&&(strcmp(SolFound,'True'))
        saveobservables_extrap(sol, la, K, controlpar, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 1, movie01, paramextrap, ControlPar, Profile,profilewidth);
    end
   
    %% stop simulation if shape is self-intersecting
    if intersec==1
        disp(strcat('Shape has self-intersected. Stopping simulation at controlpar = ',num2str(controlpar)));
        saveobservables(sol, la, K, controlpar, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotopt, 1, movie01, ControlPar, Profile, profilewidth);
        break;
    end
    
end
end

fclose('all');
close all;
cd('../')

end
