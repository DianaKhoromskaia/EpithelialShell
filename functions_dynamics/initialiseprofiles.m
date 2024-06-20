function [zeta, dszeta, zetac, dszetac, zetanem, dszetanem , zetacnem, dszetacnem, dir2, zeta_controls, zeta_profiles, zeta_consts, zeta_las, zeta_facs, zeta_sigmas, zeta_thalfs, N_regions, write9, write91, write92, write93] = initialiseprofiles(ProfileFile, sgrid, svec1, npoints, L0, s0, Q, t, tsigma, dir1, zetasrect)
% initialises profiles of active contributions to tensions or moments

zetavec = zeros(size(sgrid));
zetacvec = zeros(size(sgrid));
zetanemvec = zeros(size(sgrid));
zetacnemvec = zeros(size(sgrid));
ds = L0/(npoints-1);
onevec = ones(size(sgrid));

write9=0;
write91=0;
write92=0;
write93=0;

zeta_controls=[];
zeta_profiles=[];
zeta_consts=[];
zeta_las=[];
zeta_facs=[];
zeta_sigmas=[];
zeta_thalfs=[];
N_regions=0;

fid=fopen(ProfileFile);
tline=fgetl(fid);

if tline==-1
    fprintf('\n%s\n','<strong>Warning:</strong> Profiles file is empty')
end

dir2=dir1;

while tline ~= -1
    N_regions=N_regions+1;
    tline=split(tline);
    
    % get the active profile parameters
    zeta_control=tline{1};
    zeta_profile=tline{2};
    zeta_const=str2double(tline{3});
    zeta_la=str2double(tline{4});
    zeta_fac=str2double(tline{5});
    zeta_sigma=str2double(tline{6});
    zeta_thalf=str2double(tline{7});

    zeta_controls=[zeta_controls, {zeta_control}];
    zeta_profiles=[zeta_profiles, {zeta_profile}];
    zeta_consts=[zeta_consts, zeta_const];
    zeta_las=[zeta_las, zeta_la];
    zeta_facs=[zeta_facs, zeta_fac];
    zeta_sigmas=[zeta_sigmas, zeta_sigma];
    zeta_thalfs=[zeta_thalfs, zeta_thalf];
    
    if strcmp(zeta_profile,'Sigmoidal') && zeta_la==1
        dir2=strcat(dir2,'_',zeta_control,'_',zeta_profile,'_la=',num2str(zeta_la),'_fac=',num2str(zeta_const),'_sigma=',num2str(zeta_sigma),'_thalf=',num2str(zeta_thalf));
    else
        dir2=strcat(dir2,'_',zeta_control,'_',zeta_profile,'_la=',num2str(zeta_la),'_fac=',num2str(zeta_fac),'_sigma=',num2str(zeta_sigma),'_thalf=',num2str(zeta_thalf));
    end

    % define zeta profile of the considered region
    switch zeta_profile
        case 'Gaussian'
                zeta = (1-sigmoidal(t,zeta_thalf,tsigma))*superGaussian(s0(sgrid), zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, 1);
        case 'Sigmoidal'
            if zeta_la==1
                zeta = (1-sigmoidal(t,zeta_thalf,tsigma))*zeta_const*onevec;  
            else
                zeta = (1-sigmoidal(t,zeta_thalf,tsigma))*(zeta_const*onevec+zeta_fac*sigmoidal(s0(sgrid), zeta_la*L0, zeta_sigma*L0)); 
            end
        case 'Rectangle'
                zeta=(1-sigmoidal(t,zeta_thalf,tsigma))*rect(s0(sgrid),zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, zetasrect*L0);
        case 'Linear'
                zeta=(1-sigmoidal(t,zeta_thalf,tsigma))*linear(s0(sgrid),zeta_la*L0, zeta_fac, zeta_const, zeta_sigma*L0, zetasrect*L0);
        case 'Exponential'
                zeta=(1-sigmoidal(t,zeta_thalf,tsigma))*(zeta_const*onevec+zeta_fac*exponential(s0(sgrid),zeta_la*L0,zeta_sigma*L0));
    end

    % add to the global profile
    switch zeta_control 
        case 'Tension'
            zetavec = zetavec+zeta;
            write9=1;
        case 'Bending'
            zetacvec = zetacvec+zeta;
            write91=1;
        case 'Nematic'
            zetanemvec = zetanemvec+zeta;
            write92=1;
        case 'BendingNematic'
            zetacnemvec = zetacnemvec+zeta;
            write93=1;
    end

    tline=fgetl(fid); % read next line
end

% Nematic profiles:
zetanemvec = zetanemvec.*Q(sgrid);
zetacnemvec = zetacnemvec.*Q(sgrid);

% save as interpolants:
zeta = griddedInterpolant(sgrid, zetavec, 'spline');
zetac = griddedInterpolant(sgrid, zetacvec, 'spline');
zetanem = griddedInterpolant(sgrid, zetanemvec, 'spline');       
zetacnem = griddedInterpolant(sgrid, zetacnemvec, 'spline');                   


% define dszeta(s) profile: (isotropic tension)
dszetavec = gradient(zeta(svec1),ds); 

% define dszetac(s) profile: (isotropic bending)
dszetacvec = gradient(zetac(svec1),ds); 

% define dszetanem(s) profile: (nematic tension)
dszetanemvec = gradient(zetanem(svec1),ds); 

% define dszetacnem(s) profile: (nematic bending)
dszetacnemvec = gradient(zetacnem(svec1),ds); 


dszeta = griddedInterpolant(svec1, [0 dszetavec(2:end)], 'spline'); 
dszetac = griddedInterpolant(svec1, [0 dszetacvec(2:end)], 'spline'); 
dszetanem = griddedInterpolant(svec1, [0 dszetanemvec(2:end)], 'spline'); 
dszetacnem = griddedInterpolant(svec1, [0 dszetacnemvec(2:end)], 'spline'); 

fclose(fid);

end

