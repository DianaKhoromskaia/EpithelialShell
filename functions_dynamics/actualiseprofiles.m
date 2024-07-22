function [kappaNew, dskappaNew, zetaNew, dszetaNew, zetacNew, dszetacNew, zetanemNew, dszetanemNew, zetacnemNew, dszetacnemNew] = actualiseprofiles(zeta_controls, zeta_profiles, zeta_implementation_types, zeta_consts, zeta_las, zeta_facs, zeta_sigmas, zeta_thalfs, Qnew, svec1, sgrid, sgridnem, s0new, N_regions, L0, zetasrect, t, dt, tsigma, ds, kappa0, write94)
% update active profiles

    kappavec= zeros(size(sgridnem));
    zetavec= zeros(size(sgrid));
    zetacvec= zeros(size(sgrid)); 
    zetanemvec= zeros(size(sgridnem));
    zetacnemvec= zeros(size(sgridnem));

    for i=1:N_regions
        zeta_control=zeta_controls{i};
        zeta_profile=zeta_profiles{i};
        zeta_implementation_type=zeta_implementation_types{i};
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

        % choose the profile coordinates

        switch zeta_implementation_type
            case 'Lagrangian'
                sprofile=s0new(svec);
                L=L0;
            case 'Eulerian'
                sprofile=svec;
                L=svec(end);
        end

        if zeta_thalf == 0
            prefac_t=1;
        else
            prefac_t=(1-sigmoidal(t,zeta_thalf,tsigma));
        end
        
        % calculate zeta profile of the considered region
        switch zeta_profile
            case 'Gaussian'
                zeta = prefac_t*superGaussian(sprofile, zeta_la*L, zeta_fac, zeta_const, zeta_sigma*L, 1);
            case 'Sigmoidal'
                if zeta_la==1
                    zeta = prefac_t*zeta_const*ones(size(svec));
                else
                    zeta = prefac_t*(zeta_const*ones(size(svec))+zeta_fac*sigmoidal(sprofile, zeta_la*L, zeta_sigma*L));
                end
            case 'Rectangle'
                zeta=prefac_t*rect(sprofile,zeta_la*L, zeta_fac, zeta_const, zeta_sigma*L, zetasrect*L);
            case 'Linear'
                zeta=prefac_t*linear(sprofile,zeta_la*L, zeta_fac, zeta_const, zeta_sigma*L, zetasrect*L);
            case 'Exponential'
                zeta=prefac_t*(zeta_const*ones(size(svec))+zeta_fac*exponential(sprofile,zeta_la*L,zeta_sigma*L));
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
    
    dszetaNew = griddedInterpolant(svec1, dszetavec, 'spline');
    dszetacNew = griddedInterpolant(svec1, dszetacvec, 'spline'); 
    dszetanemNew = griddedInterpolant(svec1, dszetanemvec, 'spline'); 
    dszetacnemNew = griddedInterpolant(svec1, dszetacnemvec, 'spline'); 
    dskappaNew = griddedInterpolant(svec1, dskappavec, 'spline');
end