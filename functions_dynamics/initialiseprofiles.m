function [zeta, dszeta, zetac, dszetac, zetanem, dszetanem , zetacnem, dszetacnem] = initialiseprofiles(zeta_Profile, zeta_const, zeta_la, zeta_sigma, zeta_par, zetac_Profile, zetac_const, zetac_la, zetac_sigma, zetac_par, zetanem_Profile, zetanem_const, zetanem_la, zetanem_sigma, zetanem_par, zetacnem_Profile, zetacnem_const, zetacnem_la, zetacnem_sigma, zetacnem_par, sgrid, svec1, L, npoints, L0, s0, Q, dsQ, t, thalf, tsigma)
%% Initialising profiles of active tensions and moments and their spatial derivatives as spline interpolants.

zetavec = zeros(size(sgrid));
dszetavec = zeros(size(svec1));
zetacvec = zeros(size(sgrid));
dszetacvec = zeros(size(svec1));
zetanemvec = zeros(size(sgrid));
dszetanemvec = zeros(size(svec1));
zetacnemvec = zeros(size(sgrid));
dszetacnemvec = zeros(size(svec1));
h = L/(npoints-1);
onevec = ones(size(sgrid));

% define zeta(s) profile: (isotropic tension)
switch zeta_Profile
    case 'Gaussian'
            zetavec = superGaussian(s0(sgrid), zeta_la*L0, zeta_par, zeta_const, zeta_sigma*L0, 1);
    case 'Sigmoidal'
        if zeta_la==1
            zetavec = zeta_const*onevec;  
        else
            zetavec = zeta_const*onevec+zeta_par*sigmoidal(s0(sgrid), zeta_la*L0, zeta_sigma*L0); 
        end
end
  
% define zetac(s) profile: (isotropic bending moment)
switch zetac_Profile
    case 'Gaussian'
            zetacvec = superGaussian(s0(sgrid), zetac_la*L0, zetac_par, zetac_const, zetac_sigma*L0, 1);
    case 'Sigmoidal'
        if zetac_la==1
            zetacvec = zetac_const*onevec;  
        else
            zetacvec = zetac_const*onevec+zetac_par*sigmoidal(s0(sgrid), zetac_la*L0, zetac_sigma*L0); 
        end
end

% define zetanem(s) profile: (nematic tension)
switch zetanem_Profile
    case 'Gaussian'
        zetanemvec = Q(sgrid).*superGaussian(s0(sgrid), zetanem_la*L0, zetanem_par, zetanem_const, zetanem_sigma*L0, 1);
    case 'Sigmoidal'
        if zetanem_la==1
            zetanemvec = zetanem_const*Q(sgrid);
        else
            zetanemvec = zetanem_const*Q(sgrid)+zetanem_par*Q(sgrid).*sigmoidal(s0(sgrid), zetanem_la*L0, zetanem_sigma*L0);
        end
end

% define zetacnem(s) profile: (nematic bending moment)
switch zetacnem_Profile
    case 'Gaussian'
        zetacnemvec = Q(sgrid).*superGaussian(s0(sgrid), zetacnem_la*L0, zetacnem_par, zetacnem_const, zetacnem_sigma*L0, 1);
    case 'Sigmoidal'
        if zetacnem_la==1
            zetacnemvec = zetacnem_const*Q(sgrid);
        else
            zetacnemvec = zetacnem_const*Q(sgrid)+zetacnem_par*Q(sgrid).*sigmoidal(s0(sgrid), zetacnem_la*L0, zetacnem_sigma*L0);
        end
end

% time sigmoidal:
zetavec = (1-sigmoidal(t,thalf,tsigma))*zetavec;
zetacvec = (1-sigmoidal(t,thalf,tsigma))*zetacvec;
zetanemvec = (1-sigmoidal(t,thalf,tsigma))*zetanemvec;
zetacnemvec = (1-sigmoidal(t,thalf,tsigma))*zetacnemvec;

% save as interpolants:
zeta = griddedInterpolant(sgrid, zetavec, 'spline');
zetac = griddedInterpolant(sgrid, zetacvec, 'spline');
zetanem = griddedInterpolant(sgrid, zetanemvec, 'spline');       
zetacnem = griddedInterpolant(sgrid, zetacnemvec, 'spline');                   


% define dszeta(s) profile: (isotropic tension)
if (zeta_la<1)||strcmp(zeta_Profile,'Gaussian')
    dszetavec = gradient(zeta(svec1),h); 
end

% define dszetac(s) profile: (isotropic bending)
if (zetac_la<1)||strcmp(zetac_Profile,'Gaussian')
    dszetacvec = gradient(zetac(svec1),h); 
end

% define dszetanem(s) profile: (nematic tension)
if (zetanem_la==1)&&strcmp(zetanem_Profile,'Sigmoidal')
    dszetanemvec = zetanem_const*dsQ(svec1);
    dszetanemvec = (1-sigmoidal(t,thalf,tsigma))*dszetanemvec;
else    
    dszetanemvec = gradient(zetanem(svec1),h); 
end

% define dszetacnem(s) profile: (nematic bending)
if (zetacnem_la==1)&&strcmp(zetacnem_Profile,'Sigmoidal')
    dszetacnemvec = zetacnem_const*dsQ(svec1);
    dszetacnemvec = (1-sigmoidal(t,thalf,tsigma))*dszetacnemvec;
else    
    dszetacnemvec = gradient(zetacnem(svec1),h); 
end

dszeta = griddedInterpolant(svec1, [0 dszetavec(2:end-1) 0], 'spline'); 
dszetac = griddedInterpolant(svec1, [0 dszetacvec(2:end-1) 0], 'spline'); 
dszetanem = griddedInterpolant(svec1, [0 dszetanemvec(2:end-1) 0], 'spline'); 
dszetacnem = griddedInterpolant(svec1, [0 dszetacnemvec(2:end-1) 0], 'spline'); 

end

