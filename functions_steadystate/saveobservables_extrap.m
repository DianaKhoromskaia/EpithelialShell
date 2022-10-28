function saveobservables_extrap(sol, la, K, controlpar, tcomp, intersec, epsstep, fileID, formatSpec, NonParseq, stepsave1, stepsave2, i, plotting, sequence, moviename, paramextrap, ControlPar, Profile, width)
% this function takes the sol-structure from the ode-solver and saves
% measurements to file.
% options for plotting:
% 'shape' - only shape is plotted
% 'shapemovie' - only shape is plotted and saved as a movie
% 'all' - shape, functions, and parametric curve are plotted
% 'allmovie' - shape, functions, and parametric curve are plotted and saved as a movie
% 'off'- no plotting

%% calculate jump across the boundary in all functions: (all should be zero)
    if (la<1)&&strcmp(Profile,'Step')
        delta_nonmatch_abs = (deval(sol,1+epsstep)-deval(sol,1-epsstep));
        delta_nonmatch_rel = (deval(sol,1+epsstep)-deval(sol,1-epsstep))./deval(sol, 1);
        boundary_val = deval(sol, 1); %values at red-blue boundary
    end

%% evaluate solution: parameters 
    P = sol.parameters(1);
    fc = sol.parameters(2);
    dsw0 = sol.parameters(3);
    dswL = sol.parameters(4);
    V = sol.y(7,end);
    A = sol.y(8,end);
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        L=sol.parameters(5);
        L1=L; 
        L2=L;
    else
        L1 = sol.parameters(5);
        L2 = sol.parameters(6);
        L = L1+L2;
    end
    
%% evaluate on (almost) solver grid - except for la-point \pm epsstep:
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        solvergrid = sol.x;
        solvergridtransformed = L*sol.x;
    else
        solvergrid = [sol.x(sol.x<1) 1-epsstep 1+epsstep sol.x(sol.x>1)];
        onevecgrid = ones(size([1+epsstep sol.x(sol.x>1)]));
        solvergridtransformed = [L1*([sol.x(sol.x<1) 1-epsstep]) (L1-L2)*onevecgrid+L2*([1+epsstep sol.x(sol.x>1)])];
    end
    indred = find(solvergrid<=1); 
    indla = indred(end);
    indblue = find(solvergrid>1);
    [solution, derivative] = deval(sol, solvergrid);
    onevec = ones(size(solvergrid));
    %solvergrid_left = solvergrid(solvergrid<=la);
    %solvergrid_right = solvergrid(solvergrid<=la);
       
    x = solution(5,:);
    z = solution(6,:);
    tss = solution(1,:);
    tns = solution(2,:);
    psi = solution(4,:);
    s0 = solution(9,:);
    q = solution(10,:);
    dsq = solution(11,:);
    I1 = solution(12,:);
    I1sol = griddedInterpolant(solvergrid, I1, 'spline');
    xsol = griddedInterpolant(solvergrid, x, 'spline');
    
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
    C_0 = C(1);
    C_L = C(end);
        
    %% extrema of absolut principal curvatures
    [c1max, ind_c1max] = max(abs(c1));
    scaledposition_c1max = sol.x(ind_c1max);
    [c2max, ind_c2max] = max(abs(c2));
    scaledposition_c2max = sol.x(ind_c2max);
    if (la==1)||strcmp(Profile,'Sigmoidal')||strcmp(Profile,'Gaussian')
        position_c1max = L*scaledposition_c1max;
        position_c2max = L*scaledposition_c2max;
    else
        if scaledposition_c1max<=1
            position_c1max = L1*scaledposition_c1max;
        else
            position_c1max = L1-L2+L2*scaledposition_c1max;
        end
        if scaledposition_c2max<=1
            position_c2max = L1*scaledposition_c2max;
        else
            position_c2max = L1-L2+L2*scaledposition_c2max;
        end
    end

%% test shape for intersection
%     Xvec = sol_sample(5,:);
%     Zvec = sol_sample(6,:);
%     [X0, Y0, I, J] = intersections(Xvec, Zvec);
%     intersec = int8(~isempty(X0));   
    
%% save to file on solver grid

    if i==1
        fprintf(fileID, formatSpec, controlpar, P, L, L1, L2, V, A, fc, C_0, C_L, c2max, position_c2max, c1max, position_c1max, z(end), tcomp, sol.stats.nmeshpoints, intersec, dsw0, dswL, indla, sequence, 2*pi*I1sol(0.5), xsol(0.5));
        dlmwrite('sgrid.dat', solvergrid, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('sgridtransf.dat', solvergridtransformed, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('x.dat', x, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('z.dat', z, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('tss.dat', tss, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('mss.dat', mss, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('tns.dat', tns, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('c1.dat', c1, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('c2.dat', c2, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('C.dat', C, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('psi.dat', psi, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('s0.dat', s0, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('fs.dat', fs, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('fphi.dat',fphi, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('q.dat', q, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('dsq.dat', dsq, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('I1.dat',I1, 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('u.dat',u, 'precision', '%6.5f' ,'delimiter', '\t');
        if (la<1)&&strcmp(Profile,'Step')
            dlmwrite('deltas_abs.dat',delta_nonmatch_abs', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite('deltas_rel.dat',delta_nonmatch_rel', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite('boundaryval.dat',boundary_val', 'precision', '%10.9f' ,'delimiter', '\t');
        end
        
    else % ((i<NonParseq)&&(mod(i,stepsave1)==0))||(i==NonParseq)||((i>NonParseq)&&(mod(i,stepsave2)==0))
        fprintf(fileID, formatSpec, controlpar, P, L, L1, L2, V, A, fc, C_0, C_L, c2max, position_c2max, c1max, position_c1max, z(end), tcomp, sol.stats.nmeshpoints, intersec, dsw0, dswL, indla, sequence, 2*pi*I1sol(0.5), xsol(0.5));
        dlmwrite('sgrid.dat', solvergrid, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('sgridtransf.dat', solvergridtransformed, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('x.dat', x, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('z.dat', z, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('tss.dat', tss, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('mss.dat', mss, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('tns.dat', tns, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('c1.dat', c1, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('c2.dat', c2, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('C.dat', C, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('psi.dat', psi, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('s0.dat', s0, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('fs.dat', fs, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('fphi.dat', fphi, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('q.dat', q, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('dsq.dat', dsq, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('I1.dat', I1, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        dlmwrite('u.dat', u, '-append', 'precision', '%6.5f' ,'delimiter', '\t');
        if (la<1)&&strcmp(Profile,'Step')
            dlmwrite('deltas_abs.dat',delta_nonmatch_abs', '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite('deltas_rel.dat',delta_nonmatch_rel', '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite('boundaryval.dat',boundary_val', '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        end
        
    end

%% produce plots if required
    switch plotting
        case 'shape'
            figure(1)
            subplot(1,2,1)
%             plot(Xvec,Zvec,'b-');hold on;
%             plot(-Xvec,Zvec,'b-');hold off; 
            plot(x(1:indla), z(1:indla), 'r-');hold on;
            plot(-x(1:indla), z(1:indla), 'r-');
            if (la<1)&&strcmp(Profile,'Step')
                plot(x(indla:end), z(indla:end), 'b-');
                plot(-x(indla:end), z(indla:end), 'b-');
            end
            hold off;
            axis([-3 3 0 6])
            axis square
            subplot(1,2,2)
%             plot(solvergrid, c2);
%             title('meridional curvature');
            plot(controlpar, sol.parameters(1), 'r.'); hold on;
            if i>=5
                plot(controlpar, paramextrap(1), 'ko');
            end
            axis square
            set(findobj(gcf,'Type','Line'), 'LineWidth',2);
            %plot(controlpar, P,'.'); hold on;
            %subplot(1,3,3)
            %plot(solvergrid, q, solvergrid, dsq);
            
%              figure(2)
%              subplot(2,1,1)
%              plot(solvergrid, [derivative(10,indred)/L1 derivative(10,indblue)/L2]/pi, solvergrid, solution(11,:))
%             plot(solvergrid,fs,solvergrid,fphi, solvergrid, derivative(9,:), 'o-', s_sample, solp_sample(9,:),'*') %[derivative(9,indred)/L1 derivative(9,indblue)/L2]
%             subplot(2,1,2)
%             plot(solvergrid,u,solvergrid,utest1,solvergrid,utest2)
        case 'movie'
            figure(1)
            subplot(1,2,1)
%             plot(Xvec,Zvec,'b-');hold on;
%             plot(-Xvec,Zvec,'b-');hold off; 
            plot(x(1:indla), z(1:indla), 'r-');hold on;
            plot(-x(1:indla), z(1:indla), 'r-');
            if (la<1)&&strcmp(Profile,'Step')
                plot(x(indla:end), z(indla:end), 'b-');
                plot(-x(indla:end), z(indla:end), 'b-');
            end
            hold off;
            axis([-5 5 -5 5])
            axis square
            subplot(1,2,2)
%             plot(solvergrid, c2);
%             title('meridional curvature');
            plot(controlpar, sol.parameters(1), 'r.'); hold on;
            if i>=5
                plot(controlpar, paramextrap(1), 'ko');
            end
            axis square
            set(findobj(gcf,'Type','Line'), 'LineWidth',2);
            writeVideo(moviename, getframe(gcf));
            
        case 'all'
        case 'allmovie'
        case 'off'
    end
    
end

