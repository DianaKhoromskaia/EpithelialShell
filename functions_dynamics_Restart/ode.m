function dvds = ode(s, v, varpar, U, dsU, C1, C2, C, C0, dsC, dsC1, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, dszetac, zetanem, dszetanem, zetacnem, dszetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, dskappa, K, xi, FixedPar, t, P0, thalf_P, tsigma, zc0_on, etacb_DC_on)
% rhs for ode for the vector v = (dsvs, vn, dsvn, mss, tns, dV(s), dX(s),
% vs, dsnew(s), I(s)), solved on interval [0,L]
%figure(60)
%plot(s,v)
%legend('dsvs', 'vn', 'dsvn', 'mss', 'tns', 'dV(s)', 'dX(s)','vs', 'dsnew(s)', 'I')


%DEBUG
persistent nCalls t0
if isempty(nCalls)
    nCalls = 0;
    t0 = tic;
end
nCalls = nCalls + 1;

if mod(nCalls,2000)==0
    fprintf('ode calls: %d, elapsed %.1f s, size(s)=[%d %d]\n', ...
        nCalls, toc(t0), size(s,1), size(s,2));
end
%EDEBUG

if strcmp(FixedPar,'V')
    P = varpar(1);
elseif strcmp(FixedPar,'P')
    P = P0*(1-sigmoidal(t,thalf_P,tsigma))-etap*varpar(1); %here varpar(1)=dV/dt
end

fc=varpar(2);

dsvs = v(1,:);
vn = v(2,:);
dsvn = v(3,:);
mss = v(4,:);
tns = v(5,:);
%v0 = v(6,:);
%x0 = v(7,:);
vs = v(8,:);
% dsnew = v(9,:);
% I1 = v(10,:);

T = zetanem(s);
V = eta*(cos(Psi(s)).*vs./X(s)-dsvs+(C1(s)-C2(s)).*vn);

% geometry pieces: Cs^s - Cphi^phi and d_s Cs^s
Cs_minus_Cphi = C2(s) - C1(s);
dsCs          = dsC(s) - dsC1(s);   % since C = C1 + C2  ->  dsC2 = dsC - dsC1

% DC_kk/Dt definition from m_s^s constitutive equation
DCkk = (mss - 2*kappa(s).*(C(s)-C0) - zetac(s) + zetacnem(s))/etacb;

% bulk bending viscosity term for d2svs 
ds2vs_DCkk_term = ...
    ( (mss - 2*kappa(s).*(C(s)-C0) - zetac(s) + zetacnem(s)) ) .* ...
    ( cos(Psi(s)).*(C2(s)-C1(s))./X(s) + (dsC(s)-dsC1(s)) ) ...
  + ( 2*cos(Psi(s)).*zetacnem(s)./X(s) ...
      + tns ...
      - 2*dskappa(s).*(C(s)-C0) ...
      - 2*kappa(s).*dsC(s) ...
      - dszetac(s) ...
      + dszetacnem(s) ) .* C2(s);


fextn = -xi*vn;

onevec = ones(size(s));

ds2vn = -cos(Psi(s)).*dsvn./X(s) - vn.*(C1(s).^2 + C2(s).^2) + vs.*dsC(s) - (mss - 2*kappa(s).*(C(s)-C0) - zetac(s) + zetacnem(s))/etacb;

tss = 2*K*U(s) + zeta(s) - zetanem(s) -(kappa(s).*(C(s)-C0)+zc0_on*0.5*zetac(s)).*(C2(s)-C1(s)) + (zc0_on*0.5*zetac(s)-C0*kappa(s)).*(C(s)-C0+zc0_on*0.5*zetac(s)./kappa(s)) + (eta+etab)*dsvs + (etab-eta)*cos(Psi(s)).*vs./X(s) + (etab*C(s)+eta*(C2(s)-C1(s))).*vn - etacb_DC_on*etacb*DCkk.*C2(s);
dstns = 2*C1(s).*(T+V+(kappa(s).*(C(s)-C0)+zc0_on*0.5*zetac(s)).*(C2(s)-C1(s))) + C(s).*tss - cos(Psi(s)).*tns./X(s) - P + cos(Psi(s)).*fc  - fextn + etacb_DC_on*etacb*DCkk .* C1(s) .* (C2(s) - C1(s));
dsmss = tns + 2*cos(Psi(s)).*zetacnem(s)./X(s);

ds2vs = -cos(Psi(s)).*(dsvs-cos(Psi(s)).*vs./X(s))./X(s) - (eta-etab)*C1(s).*C2(s).*vs/(eta+etab) - dsC(s).*vn - (etab*C(s)+eta*(C2(s)-C1(s))).*dsvn/(eta+etab) - (2*K*dsU(s)+dszeta(s))/(eta+etab)  + (dszetanem(s) + 2*cos(Psi(s)).*zetanem(s)./X(s))/(eta+etab) - C2(s).*tns./(eta+etab) - sin(Psi(s)).*fc/(eta+etab) + (dskappa(s).*C(s).*(C2(s)-C1(s)) + 2*C2(s).*kappa(s).*dsC(s)+ zc0_on*C1(s).*dszetac(s) + (zc0_on*0.5*zetac(s)/kappa(s)-C0).*(zc0_on*dszetac(s)+(C0-zc0_on*0.5*zetac(s)/kappa(s)).*dskappa(s)/kappa(s)))./(eta+etab) + etacb_DC_on*ds2vs_DCkk_term./(eta+etab);


% force balance within the interval:
dvds = [ds2vs;
    dsvn;
    ds2vn;
    dsmss;
    dstns;
    2*pi*X(s).*vn;
    X(s).*vn.*(C(s).*(Z(s)-X0) - cos(Psi(s)))/xintegral;
    dsvs;
    dsvs + C2(s).*vn;%sqrt(onevec+2*dt*(dsvs + C2(s).*vn)+(dt^2)*((dsvs + C2(s).*vn).^2+(dsvn-C2(s).*vs).^2));%dsvs + C2(s).*vn;
    X(s).*(fc.*onevec - cos(Psi(s)).*fextn)];

%% Pole overwrite (same pole BC/regularity structure as OLD code)
% Assumption:
% - pole regularity conditions are unchanged
% - the only new pole contribution from eta_cb enters through t_s^s at the pole:
%       t_s^s -> t_s^s - etacb_DC_on * etacb * DCkk * C_s^s
% - therefore only dvds(5) changes relative to OLD, via tssPole
%
% With etacb_DC_on = 0, this reproduces the OLD pole overwrite exactly.

tolPole = 1e-12;

% --- South pole (s = 0)
idxSP = find(abs(s) < tolPole);
if ~isempty(idxSP)

    % values at s = 0
    C10   = C1(0);
    C20   = C2(0);
    C0tot = C(0);
    U0    = U(0);
    z0    = zeta(0);
    zc0p  = zetac(0);
    k0    = kappa(0);

    dz0   = dszeta(0);
    dzc0  = dszetac(0);
    dk0   = dskappa(0);

    for jj = 1:numel(idxSP)
        ii = idxSP(jj);

        one = 1;
        zer = 0;

        % --- OLD pole expressions kept exactly ---
        dv1_SP = (2/3) * ...
            ( -dz0 ...
            + zc0_on*dzc0*C10 ...
            + (zc0_on*0.5*zc0p/k0 - C0) * ...
              ( zc0_on*dzc0 + (C0 - zc0_on*0.5*zc0p/k0)*dk0/k0 ) ) ...
            / (eta + etab);

        dv3_SP = 0.5 * ...
            ( -vn(ii)*(C10^2 + C20^2) ...
            - (mss(ii) - 2*k0*(C0tot - C0) - zc0p)/etacb );

        dv4_SP = tns(ii);

        % OLD pole t_s^s
        tssPole_old_SP = ...
            2*K*U0 + z0 ...
            + (zc0_on*0.5*zc0p - k0*C0) * ...
              (C0tot - C0 + zc0_on*0.5*zc0p/k0) ...
            + 2*etab*dsvs(ii) + etab*C0tot*vn(ii);

        % NEW eta_cb correction enters only through t_s^s at the pole
        % Use the same DCkk definition as in the bulk/new code.
        DCkk0 = (mss(ii) - 2*k0*(C0tot - C0) - zc0p + zetacnem(0)) / etacb;

        tssPole_SP = tssPole_old_SP ...
                   - etacb_DC_on * etacb * DCkk0 * C20;

        % Same pole regularity formula as OLD code, but with updated tssPole
        dv5_SP = -0.5*P*one + C20*tssPole_SP + 0.5*xi*vn(ii) + 0.5*fc*one;

        dvds(:,ii) = [ ...
            dv1_SP; ...
            zer; ...
            dv3_SP; ...
            dv4_SP; ...
            dv5_SP; ...
            zer; ...
            zer; ...
            dsvs(ii); ...
            dsvs(ii) + C20*vn(ii); ...
            zer ];
    end
end

% --- North pole (s = L)
idxNP = find(abs(s - L) < tolPole);
if ~isempty(idxNP)

    % values at s = L
    C1L   = C1(L);
    C2L   = C2(L);
    CLtot = C(L);
    UL    = U(L);
    zL    = zeta(L);
    zcL   = zetac(L);
    kL    = kappa(L);

    dzL   = dszeta(L);
    dzcL  = dszetac(L);
    dkL   = dskappa(L);

    for jj = 1:numel(idxNP)
        ii = idxNP(jj);

        one = 1;
        zer = 0;

        % --- OLD pole expressions kept exactly ---
        dv1_NP = (2/3) * ...
            ( -dzL ...
            + zc0_on*dzcL*C1L ...
            + (zc0_on*0.5*zcL/kL - C0) * ...
              ( zc0_on*dzcL + (C0 - zc0_on*0.5*zcL/kL)*dkL/kL ) ) ...
            / (eta + etab);

        dv3_NP = 0.5 * ...
            ( -vn(ii)*(C1L^2 + C2L^2) ...
            - (mss(ii) - 2*kL*(CLtot - C0) - zcL)/etacb );

        dv4_NP = tns(ii);

        % OLD pole t_s^s
        tssPole_old_NP = ...
            2*K*UL + zL ...
            + (zc0_on*0.5*zcL - kL*C0) * ...
              (CLtot - C0 + zc0_on*0.5*zcL/kL) ...
            + 2*etab*dsvs(ii) + etab*CLtot*vn(ii);

        % NEW eta_cb correction enters only through t_s^s at the pole
        DCkkL = (mss(ii) - 2*kL*(CLtot - C0) - zcL + zetacnem(L)) / etacb;

        tssPole_NP = tssPole_old_NP ...
                   - etacb_DC_on * etacb * DCkkL * C2L;

        % Same pole regularity formula as OLD code, but with updated tssPole
        dv5_NP = -0.5*P*one + C2L*tssPole_NP + 0.5*xi*vn(ii) - 0.5*fc*one;

        dvds(:,ii) = [ ...
            dv1_NP; ...
            zer; ...
            dv3_NP; ...
            dv4_NP; ...
            dv5_NP; ...
            zer; ...
            zer; ...
            dsvs(ii); ...
            dsvs(ii) + C2L*vn(ii); ...
            zer ];
    end
end


end