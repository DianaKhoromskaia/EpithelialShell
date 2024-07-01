function [kappaNew, dskappaNew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew] = actualiseprofiles(zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_facs, zeta_sigmas, zeta_thalfs, Qnew, svec1, sgrid, sgridnem, s0new, N_regions, L0, zetasrect, t, dt, tsigma, ds, kappa0, write94)
% update active profiles

    zetavec= zeros(size(sgrid));
    zetacvec= zeros(size(sgrid)); 
    zetanemvec= zeros(size(sgridnem));
    zetacnemvec= zeros(size(sgridnem));
    kappavec= zeros(size(sgridnem));

    for i=1:N_regions
        zeta_control=zeta_controls{i};
        zeta_profile=zeta_profiles{i};
        zeta_const=zeta_consts(i);
        zeta_la=zeta_las(i);
        zeta_fac=zeta_facs(i);
        zeta_sigma=zeta_sigmas(i);
        zeta_thalf=zeta_thalfs(i);

        % choose the good grid vector
        if strcmp(zeta_control,'Tension')||strcmp(zeta_control,'Bending')
            svec=sgrid;
        else
            svec=sgridnem;
        end
        
        % calculate zeta profile of the considered region
        switch zeta_profile
            case 'Gaussian'
                zeta = (1-sigmoidal(t+dt,zeta_thalf,tsigma))*superGaussian(s0new(svec), zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, 1);
            case 'Sigmoidal'
                if zeta_la==1
                    zeta = (1-sigmoidal(t+dt,zeta_thalf,tsigma))*zeta_const*ones(size(svec));
                else
                    zeta = (1-sigmoidal(t+dt,zeta_thalf,tsigma))*(zeta_const*ones(size(svec))+zeta_fac*sigmoidal(s0new(svec), zeta_la*L0, zeta_sigma*L0));
                end
            case 'Rectangle'
                zeta=(1-sigmoidal(t+dt,zeta_thalf,tsigma))*rect(s0new(svec),zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, zetasrect*L0);
            case 'Linear'
                zeta=(1-sigmoidal(t+dt,zeta_thalf,tsigma))*linear(s0new(svec),zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, zetasrect*L0);
            case 'Exponential'
                zeta=(1-sigmoidal(t+dt,zeta_thalf,tsigma))*(zeta_const*ones(size(svec))+zeta_fac*exponential(s0new(svec),zeta_la*L0,zeta_sigma*L0));
        end
        
        % add to the global profile
        switch zeta_control 
            case 'Tension'
                zetavec = zetavec+zeta;
            case 'Bending'
                zetacvec = zetacvec+zeta;
            case 'Nematic'
                zetanemvec = zetanemvec+zeta;
            case 'BendingNematic'
                zetacnemvec = zetacnemvec+zeta;
            case 'BendingModulus'
                kappavec = kappavec+zeta;
        end
    end

    zetaNew =  griddedInterpolant(sgrid, zetavec, 'spline');
    zetacNew = griddedInterpolant(sgrid, zetacvec, 'spline');
    zetanemNew = griddedInterpolant(sgrid, zetanemvec.*Qnew, 'spline');
    zetacnemNew = griddedInterpolant(sgrid, zetacnemvec.*Qnew, 'spline');
    if write94
        kappaNew = griddedInterpolant(sgrid, kappavec, 'spline'); 
    else
        kappaNew = griddedInterpolant(sgrid, kappa0*ones(size(sgrid)), 'spline');
    end
          
    % profile gradients:
    dszetavec = gradient(zetaNew(svec1),ds); 
    dszetacvec = gradient(zetacNew(svec1),ds);
    dszetanemvec = gradient(zetanemNew(svec1),ds);
    dszetacnemvec = gradient(zetacnemNew(svec1),ds);
    dskappavec = gradient(kappaNew(svec1),ds);
    
    dszetaNew = griddedInterpolant(svec1, [0 dszetavec(2:end)], 'spline');
    dszetacNew = griddedInterpolant(svec1, [0 dszetacvec(2:end)], 'spline'); 
    dszetanemNew = griddedInterpolant(svec1, [0 dszetanemvec(2:end)], 'spline'); 
    dszetacnemNew = griddedInterpolant(svec1, [0 dszetacnemvec(2:end)], 'spline');
    dskappaNew = griddedInterpolant(svec1, [0 dskappavec(2:end)], 'spline');
end