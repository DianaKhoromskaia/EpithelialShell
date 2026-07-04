function [sint, Vint, Par, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] = forcebalance(U, dsU, C1, C2, C, C0, dsC, dsC1, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, dszetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, dskappa, K, xi, optode, t, delt, Adaptive, FixedPar, P0, thalf_P, tsigma, zc0_on, etacb_DC_on, varargin)
SolFound = true;

if t==0
    sgridinit = varargin{1};

    tss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
    tns = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
    mss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');

    sgrid = sgridinit;

    if strcmp(FixedPar,'V')
        parguess = [0 0];
    elseif strcmp(FixedPar,'P')
        parguess = [0 0];
    end

    % varargin handling at t==0:
    %   old behavior: varargin{1} = sgridinit
    %   restart behavior: varargin{1} = sgridinit, varargin{2} = restartPath
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        restartPath = varargin{2};
        bvpfile = fullfile(restartPath, 'bvpsolution.dat');
        if exist(bvpfile, 'file') ~= 2
            error('Missing restart BVP file: %s', bvpfile);
        end
        solinit = bvpinit_from_saved_bvp(bvpfile, parguess);
    else
        solinit = bvpinit(sgrid, @yguessfun_steadystate, parguess, mss, tss, tns);
    end

elseif t>0
    sol = varargin{1};
    sfun = varargin{2};
    snewfun = varargin{3};

    solold = sol;

    sgrid = snewfun(solold.x);
    sgrid = [0. sgrid(2:(end-1)) L];
    sgrid = make_strictly_increasing(sgrid);

    fint = integral(@(s) X(s).*fext(s), 0, L); %#ok<NASGU>

    parguess = solold.parameters;

    solinit = bvpinit(sgrid, @yguessfun, parguess, solold, sfun);
end

try
    sol = bvp4c( @ode, ...
                 @bc, ...
                 solinit, optode, U, dsU, C1, C2, C, C0, dsC, dsC1, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, dszetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, dskappa, K, xi, FixedPar, t, P0, thalf_P, tsigma, zc0_on, etacb_DC_on);
catch ME
    disp(ME)
    disp(strcat('could not find solution at time t=',num2str(t)));
    SolFound = false;

    if t > 0
        sol = solold;
    else
        sol = solinit;
    end
end

if ~SolFound
    sint = sgrid;
    Vint = NaN(10, numel(sint));
    Par  = NaN(1,2);
    vs = []; dsvs = []; vkk = []; tss = [];
    vn = []; dsvn = []; mss = []; tns = []; ds2vn = [];
    dV = NaN; dX0 = NaN;

    Lnew_dt = L; eps1new_dt = eps1abs*L;
    sfun_dt = []; snewfun_dt = []; snewvec_dt = [];
    Lnew_dthalf = L; eps1new_dthalf = eps1abs*L;
    sfun_dthalf = []; snewfun_dthalf = []; snewvec_dthalf = [];

    return
end

sint = sol.x;
Vint = sol.y;
derivatives = sol.yp;
Par = sol.parameters;
dV = sol.y(6,end);
dX0 = sol.y(7,end);

%% saving new arc length with dt=delt
Lnew_dt = L + delt*Vint(9,end);
eps1new_dt = eps1abs*Lnew_dt;

snew_dt_raw = [0. sint(2:(end-1)) + delt*Vint(9,2:(end-1)) Lnew_dt];
sold_dt_raw = [0. sint(2:(end-1)) L];

snew_dt = make_strictly_increasing(snew_dt_raw);
sold_dt = make_strictly_increasing(sold_dt_raw);

sfun_dt = griddedInterpolant(snew_dt(:), sold_dt(:), 'linear');
snewfun_dt = griddedInterpolant(sold_dt(:), snew_dt(:), 'linear');
snewvec_dt = snew_dt_raw;

%% saving new arc length with dt=delt/2
Lnew_dthalf = L + (delt/2)*Vint(9,end);
eps1new_dthalf = eps1abs*Lnew_dthalf;

snew_dthalf_raw = [0. sint(2:(end-1)) + (delt/2)*Vint(9,2:(end-1)) Lnew_dthalf];
sold_dthalf_raw = [0. sint(2:(end-1)) L];

snew_dthalf = make_strictly_increasing(snew_dthalf_raw);
sold_dthalf = make_strictly_increasing(sold_dthalf_raw);

sfun_dthalf = griddedInterpolant(snew_dthalf(:), sold_dthalf(:), 'linear');
snewfun_dthalf = griddedInterpolant(sold_dthalf(:), snew_dthalf(:), 'linear');
snewvec_dthalf = snew_dthalf_raw;

%% interpolants of functions and derivatives
dsvs = griddedInterpolant(sint, Vint(1,:), 'spline');
vn = griddedInterpolant(sint, Vint(2,:) , 'spline');
dsvn = griddedInterpolant(sint, Vint(3,:) , 'spline');
mss = griddedInterpolant(sint, Vint(4,:), 'spline');
tns = griddedInterpolant(sint, Vint(5,:), 'spline');
vs = griddedInterpolant(sint, Vint(8,:), 'spline');
vkk = griddedInterpolant(sint, Vint(1,:)+C(sint).*Vint(2,:)+[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)], 'spline');

mss_int = Vint(4,:);
DCkk_int = (mss_int - 2*kappa(sint).*(C(sint)-C0) - zetac(sint) + zetacnem(sint)) ./ etacb;

tss = griddedInterpolant(sint, 2*K*U(sint)+zeta(sint)-zetanem(sint)-(2*kappa(sint).*(C(sint)-C0)+zc0_on*zetac(sint)).*(C2(sint)-0.5*C(sint))+(-kappa(sint)*C0+zc0_on*0.5*zetac(sint)).*(C(sint)-C0+zc0_on*0.5*zetac(sint)/kappa(sint))+(etab+eta)*Vint(1,:)+(etab*C(sint)+eta*(C2(sint)-C1(sint))).*Vint(2,:)+(etab-eta)*[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)] + etacb_DC_on*etacb*DCkk_int.*C2(sint), 'spline');
ds2vn = griddedInterpolant(sint, derivatives(3,:), 'spline');

end