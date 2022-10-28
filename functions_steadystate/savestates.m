function savestates(sol, la, controlpar, ControlPar, NonParseq, epsstep, allsave1, allsave2, i, K, zeta, zetared, zetanem, zetanemred, zetacnem, zetacnemred, C0, aM, fconst, fconstred, lc, FixedPar, Profile, width)
% this function takes the sol-structure from the ode-solver and saves
% measurements to file.

%% evaluate solution: parameters 
    P = sol.parameters(1);
    fc = sol.parameters(2);
    dsw0 = sol.parameters(3);
    dswL = sol.parameters(4);
    V = sol.y(7,end);
    A = sol.y(8,end);
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        L=sol.parameters(5);
        L1=L; % saving them redundantly, but ok..
        L2=L;
    else
        L1 = sol.parameters(5);
        L2 = sol.parameters(6);
        L = L1+L2;
    end
   
%% evaluate on solver grid:
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        solvergrid = sol.x;
    else
        solvergrid = [sol.x(sol.x<1) 1-epsstep 1+epsstep sol.x(sol.x>1)];

        indred = find(solvergrid<=1); 
        indblue = find(solvergrid>1);
    end
    [solution, derivative] = deval(sol, solvergrid);
    onevec = ones(size(solvergrid));
    %solvergrid_left = solvergrid(solvergrid<=la);
    %solvergrid_right = solvergrid(solvergrid<=la);

    x = solution(5,:);
    z = solution(6,:);
    tss = solution(1,:);
    %t = [sol_sample_left(1,:) sol_sample_right(1,:)];
    tns = solution(2,:);
    %tns = [sol_sample_left(2,:) sol_sample_right(2,:)];
    psi = solution(4,:);
    s0 = solution(9,:);
    q = solution(10,:);
    dsq = solution(11,:);
    I1 = solution(12,:);
    if la<1
        mss = solution(3,:);
        switch Profile
            case 'Sigmoidal'
                profile = intensity(s0, la, width);
                %dxprofile = dxintensity(s0, la, width);
            case 'Gaussian'
        end
        if strcmp(ControlPar,'Bending')
            if la<1
                C = 0.5*(mss - controlpar*profile);
            else 
                C = 0.5*mss; 
            end
        elseif strcmp(ControlPar,'BendingNematic')
            if la<1
                C = 0.5*(mss + controlpar*profile.*q);
            else
                C = 0.5*(mss + controlpar*q);
            end
        else
            C = 0.5*mss;
        end
    elseif la==1
        C = solution(3,:);
        switch Profile
            case 'Sigmoidal'
                profile = intensity(s0, la, width);
                %dxprofile = dxintensity(s0, la, width);
            case 'Gaussian'
        end
        if strcmp(ControlPar,'Bending')
            if la<1
                mss = 2*C + controlpar*profile;
            else 
                mss = 2*C; 
            end
        elseif strcmp(ControlPar,'BendingNematic')
            if la<1
                mss = 2*C - controlpar*profile.*q;
            else
                mss = 2*C - controlpar*q;
            end
        else
            mss = 2*C;
        end
    end
    
    c1 = [0.5*C(1) sin(solution(4,2:end-1))./solution(5,2:end-1) 0.5*C(end)];
    c2 = C-c1;
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        fs = L*onevec./derivative(9,:);
        fphi = [L/derivative(9,1) solution(5,2:end-1)./(x0(solution(9,2:end-1))) L/derivative(9,end)];
    else
        fs = [L1*onevec(indred)./derivative(9,indred) L2*onevec(indblue)./derivative(9,indblue)];
        fphi = [L1/derivative(9,1) solution(5,2:end-1)./(x0(solution(9,2:end-1))) L2/derivative(9,end)];
    end
    u = (fs.*fphi)-onevec;

%% save to file on (interval-wise) uniform grid
if strcmp(FixedPar,'V')
    dirname = strcat('SteadyState_',ControlPar,'_la=',num2str(la),'_j=',num2str(i),'_par=',num2str(controlpar),'_P=',num2str(P),'_L=',num2str(L)); %here L always the total length
elseif strcmp(FixedPar,'P')
    dirname = strcat('SteadyState_',ControlPar,'_la=',num2str(la),'_j=',num2str(i),'_par=',num2str(controlpar),'_V=',num2str(P),'_L=',num2str(L));
end
mkdir(dirname);
cd(dirname)

filename1 = 'functions.dat';
filename2 = 'parameters.dat';

    if (i==1)||((i<NonParseq)&&(mod(i,allsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,allsave2)==0))
        dlmwrite(filename1, [solvergrid; x; z; tss; tns; c1; c2; psi; s0; q; dsq; I1; u; mss; C; fs; fphi] , 'precision', '%12.11f' ,'delimiter', '\t');
        dlmwrite(filename2, [sol.parameters' V A L K zeta zetared zetanem zetanemred C0 aM fconst fconstred lc fc dsw0 dswL controlpar la i] ,'precision', '%12.11f' ,'delimiter', '\t');
    end
cd ..

end

