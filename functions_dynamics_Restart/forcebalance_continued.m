function [sint, Vint, Par, sol, vs, dsvs, vkk, tss, vn, dsvn, mss, tns, ds2vn, dV, dX0, Lnew_dt, eps1new_dt, sfun_dt, snewfun_dt, snewvec_dt, Lnew_dthalf, eps1new_dthalf, sfun_dthalf, snewfun_dthalf, snewvec_dthalf, SolFound] = forcebalance_continued(U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, eps1abs, xintegral, fext, kappa, K, xi, optode, t, delt, Adaptive, FixedPar, varargin)
SolFound = true;

solcurrent = varargin{1};
Pcurrent = varargin{2};
fccurrent = varargin{3};

parguess = [Pcurrent fccurrent];

sgridcurrent = solcurrent(end-10,:);
[sgrid, indnonzero, ic] = unique(sgridcurrent);
Lcurrent = sgrid(end);

sol_dsvs = griddedInterpolant(sgrid, solcurrent(end-9,indnonzero), 'spline');
sol_vn = griddedInterpolant(sgrid, solcurrent(end-8,indnonzero), 'spline');
sol_dsvn = griddedInterpolant(sgrid, solcurrent(end-7,indnonzero), 'spline');
sol_mss = griddedInterpolant(sgrid, solcurrent(end-6,indnonzero), 'spline');
sol_tns = griddedInterpolant(sgrid, solcurrent(end-5,indnonzero), 'spline');
sol_dV = griddedInterpolant(sgrid, solcurrent(end-4,indnonzero), 'spline');
sol_dX = griddedInterpolant(sgrid, solcurrent(end-3,indnonzero), 'spline');
sol_vs = griddedInterpolant(sgrid, solcurrent(end-2,indnonzero), 'spline');
sol_dsnew = griddedInterpolant(sgrid, solcurrent(end-1,indnonzero), 'spline');
sol_I = griddedInterpolant(sgrid, solcurrent(end,indnonzero), 'spline');

%since I haven't saved sfun(s), try to connect the to arc length by a small
%stretch: 
sgridnext = [0 (L/Lcurrent)*sgrid(2:end-1) L];

sfun = griddedInterpolant(sgridnext, sgrid, 'spline');

solinit = bvpinit(sgridnext, @yguessfun_continued, parguess, sfun, sol_dsvs, sol_vn, sol_dsvn, sol_mss, sol_tns, sol_dV, sol_dX, sol_vs, sol_dsnew, sol_I);

% if t==0
%     sgridinit = varargin{1};
%     tss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
%     tns = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
%     mss = griddedInterpolant(sgridinit, 0*sgridinit, 'spline');
%     
%     sgrid = sgridinit;
% 
%     %with fc as second parameter:
%     if strcmp(FixedPar,'V')
%         parguess = [0 0]; 
%     elseif strcmp(FixedPar,'P')
%         parguess = [0 0];
%     end
%     solinit = bvpinit(sgrid, @yguessfun_steadystate, parguess, mss, tss, tns);
%     
% elseif t>0
%     sol = varargin{1};
%     sfun = varargin{2};
%     snewfun = varargin{3};
%     
%     solold = sol;
%     
%     sgrid = snewfun(solold.x);
%     %eps1=1e-2*L;
%     %eps2=1e-2*L;
%     sgrid = [0. sgrid(2:(end-1)) L];
%     %sgrid = [0 sgrid((sgrid > eps1)&(sgrid < (L-eps2))) L];
%     
%     fint = integral(@(s) X(s).*fext(s), 0, L);
%     
%     parguess = solold.parameters;
%     
%     solinit = bvpinit(sgrid, @yguessfun, parguess, solold, sfun);
% end

try
    sol = bvp4c( @ode, ...
                 @bc, ...
                 solinit,optode, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, delt);
    if strcmp(Adaptive,'On')
        sol_half = bvp4c( @ode, ...
                 @bc, ...
                 solinit,optode, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, delt/2);
    end
catch ME
    ME
    disp(strcat('could not find solution at time t=',num2str(t)));
    SolFound = false;
    sol = solold; %this works only if t>0
end
    
%JacobianCheck(sol, U, dsU, C1, C2, C, dsC, Psi, X, Z, X0, L, L0, zeta, dszeta, zetac, zetanem, dszetanem, zetacnem, eta, etab, etacb, etap, xintegral, fext, kappa, K, xi, FixedPar, delt)

sint = sol.x;
Vint = sol.y;
derivatives = sol.yp;
Par = sol.parameters;
dV = sol.y(6,end);
dX0 = sol.y(7,end);

% figure(60)
% subplot(1,2,1)
% plot(sint, zetacnem(sint), '.-');%, sint, dszetacnem(sint));
% legend('\zeta_{cn}');
% subplot(1,2,2)
% plot(sint, dsU(sint));%, sint, dszetacnem(sint));
% legend('dsU');

% figure(31)
% plot(sint, Vint,'*');
% legend('dsvs', 'vn', 'dsvn', 'mss', 'dV(s)', 'dX(s)','vs', 'dsnew(s)', 'I')
% axis tight

%% saving new arc length with dt=delt
% Lnew_dt = L + delt*Vint(9,end);
% eps1new_dt = eps1abs*Lnew_dt;
% 
% sfun_dt = griddedInterpolant( [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt], [0. sint(2:(end-1)) L], 'spline'); % interpolant of s(snew)
% snewfun_dt = griddedInterpolant( [0. sint(2:(end-1)) L], [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt], 'spline'); % interpolant of snew(s)
% snewvec_dt = [0. sint(2:(end-1))+delt*Vint(9,2:(end-1)) Lnew_dt]; % vector snew

%% saving new arc length with dt=delt
Lnew_dt = Vint(9,end);
eps1new_dt = eps1abs*Lnew_dt;

sfun_dt = griddedInterpolant( Vint(9,:), sint, 'spline'); % interpolant of s(snew)
snewfun_dt = griddedInterpolant( sint, Vint(9,:), 'spline'); % interpolant of snew(s)
snewvec_dt = Vint(9,:); % vector snew

%% saving new arc length with dt=delt/2
if strcmp(Adaptive,'On')
    Lnew_dthalf = sol_half.y(9,end);
    eps1new_dthalf = eps1abs*Lnew_dt;

    sfun_dthalf = griddedInterpolant(sol_half.y(9,:), sol_half.x, 'spline'); % interpolant of s(snew)
    snewfun_dthalf = griddedInterpolant(sol_half.x, sol_half.y(9,:),'spline'); % interpolant of snew(s)
    snewvec_dthalf = sol_half.y(9,:); % vector snew
else
    Lnew_dthalf = Lnew_dt;
    eps1new_dthalf = eps1new_dt;

    sfun_dthalf = sfun_dt; % interpolant of s(snew)
    snewfun_dthalf = snewfun_dt; % interpolant of snew(s)
    snewvec_dthalf = snewvec_dt; % vector snew
end
%% interpolants of functions and derivatives
%tss = griddedInterpolant(sint, Vint(1,:), 'spline');
dsvs = griddedInterpolant(sint, Vint(1,:), 'spline');
vn = griddedInterpolant(sint, Vint(2,:) , 'spline');
dsvn = griddedInterpolant(sint, Vint(3,:) , 'spline');
mss = griddedInterpolant(sint, Vint(4,:), 'spline');
tns = griddedInterpolant(sint, Vint(5,:), 'spline');
vs = griddedInterpolant(sint, Vint(8,:), 'spline');
%vkk = griddedInterpolant(sint, (Vint(1,:)+zetanem(sint)-2*K*U(sint)-zeta(sint))/etab, 'spline');
%dsvs = griddedInterpolant(sint, [((Vint(1,1)+zetanem(0)-2*K*U(0)-zeta(0))/etab-C(0)*Vint(2,1))/2 (Vint(1,2:end-1)+zetanem(sint(2:end-1))-2*K*U(sint(2:end-1))-zeta(sint(2:end-1)))/etab-cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1))-C(sint(2:end-1)).*Vint(2,2:end-1) ((Vint(1,end)+zetanem(L)-2*K*U(L)-zeta(L))/etab-C(L)*Vint(2,end))/2], 'spline');
vkk = griddedInterpolant(sint, Vint(1,:)+C(sint).*Vint(2,:)+[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)], 'spline');
tss = griddedInterpolant(sint, 2*K*U(sint)+zeta(sint)-zetanem(sint)+(etab+eta)*Vint(1,:)+(etab*C(sint)+eta*(C2(sint)-C1(sint))).*Vint(2,:)+(etab-eta)*[Vint(1,1) cos(Psi(sint(2:end-1))).*Vint(8,2:end-1)./X(sint(2:end-1)) Vint(1,end)], 'spline');

ds2vn = griddedInterpolant(sint, derivatives(3,:), 'spline');

%tns = (Vint(1,:).*sin(Psi(sint))-0.5*Par(1)*X(sint)+Vint(10,:)./X(sint))./cos(Psi(sint));
%Iformula = X(sint).*cos(Psi(sint)).*tns(sint)-X(sint).*sin(Psi(sint)).*tss(sint)+0.5*sol.parameters(1)*X(sint).*X(sint);

%figure(20)
%plot(sint, Vint, 'o-');

% figure(20)
% subplot(1,7,1);%hold on;
% plot(sint,Vint(5,:),'.-');%hold on;
% %plot(t,Vint(5,end),'.-');
% %plot(tns,'o');
% title('t_n^s')
% %hold off;
% axis tight
% grid on
% subplot(1,7,2)
% plot(sint,Vint(4,:),'.-');
% title('m^s_s')
% axis tight
% grid on
% subplot(1,7,3)
% plot(sint,Vint(1,:),'.-');
% title('t^s_s')
% axis tight
% grid on
% subplot(1,7,4)
% plot(sint,Vint(2,:),'.-');
% title('v_n')
% axis tight
% grid on
% subplot(1,7,5)
% plot(sint,Vint(8,:),'.-');
% title('v_s')
% axis tight
% grid on
% subplot(1,7,6)
% plot(sint,Vint(10,:),'.-');hold on;
% plot(sint(2:end-1),Iformula(2:end-1),'r-');hold off;
% title('I')
% axis tight
% grid on
% subplot(1,7,7)
% hold on;
% plot(t,abs(sol.parameters(2)/sol.parameters(1)),'b.');%Vint(10,1)
% %plot(sint(2:end-1),Vint(10,2:end-1)./(X(sint(2:end-1)).*X(sint(2:end-1))),'.-');
% title('|f_c/P|')
% %axis tight
% grid on

%% Check that TFB and NFB are satisfied:
% fc=Par(2)
% P=Par(1)
% TFB_rhs = -C2(sint).*tns(sint)-fc*sin(Psi(sint));
% TFB_lhs = derivatives(1,:);
% NFB_rhs = (C1(sint)+C2(sint)).*tss(sint)-P-[derivatives(5,1) (cos(Psi(sint(2:end-1)))./X(sint(2:end-1))).*tns(sint(2:end-1))  derivatives(5,end)]+fc*cos(Psi(sint));
% NFB_lhs = derivatives(5,:);
% TTB_rhs = tns(sint);
% TTB_lhs = derivatives(4,:);
% vn_rhs = -[derivatives(3,1) (cos(Psi(sint(2:end-1)))./X(sint(2:end-1))).*dsvn(sint(2:end-1)) derivatives(3,end)]-vn(sint).*(C1(sint).^2+C2(sint).^2)+vs(sint).*dsC(sint)-(mss(sint)-2*kappa*C(sint)-zetac(sint))/etacb;
% vn_lhs = derivatives(3,:);
% vs_lhs = -[derivatives(8,1) (cos(Psi(sint(2:end-1)))./X(sint(2:end-1))).*vs(sint(2:end-1)) derivatives(8,end)]-C(sint).*vn(sint)+(tss(sint)-2*K*U(sint))/etab;
% vs_rhs = derivatives(8,:);
% Ngrid=1000;
% sgrid=linspace(0,L,Ngrid);
% hgrid=L/(Ngrid-1);
% TFB_lhs_2_vec=gradient(tss(sgrid),hgrid);
% NFB_lhs_2_vec=gradient(tns(sgrid),hgrid);
% TTB_lhs_2_vec=gradient(mss(sgrid),hgrid);
% vn_lhs_2_vec=gradient(dsvn(sgrid),hgrid);
% vs_lhs_2_vec=gradient(vs(sgrid),hgrid);
% TFB_lhs_2 = griddedInterpolant(sgrid, TFB_lhs_2_vec,'spline');
% NFB_lhs_2 = griddedInterpolant(sgrid, NFB_lhs_2_vec,'spline');
% TTB_lhs_2 = griddedInterpolant(sgrid, TTB_lhs_2_vec,'spline');
% vn_lhs_2 = griddedInterpolant(sgrid, vn_lhs_2_vec,'spline');
% vs_lhs_2 = griddedInterpolant(sgrid, vs_lhs_2_vec,'spline');
% 
% figure(21)
% subplot(1,5,1)
% plot(sint,TFB_rhs-TFB_lhs);hold on;
% plot(sint,TFB_rhs-TFB_lhs_2(sint));hold off;
% title('difference TFB')
% legend({'\partial_st_s^s from solver','\partial_st_s^s from gradient'});
% subplot(1,5,2)
% plot(sint,NFB_rhs-NFB_lhs);hold on;
% plot(sint,NFB_rhs-NFB_lhs_2(sint));hold off;
% legend({'\partial_st_n^s from solver','\partial_st_n^s from gradient'});
% title('difference NFB');
% subplot(1,5,3)
% plot(sint,TTB_rhs-TTB_lhs);hold on;
% plot(sint,TTB_rhs-TTB_lhs_2(sint));hold off;
% legend({'\partial_sm_s^s from solver','\partial_sm_s^s from gradient'});
% title('difference TTB');
% subplot(1,5,4)
% plot(sint,vn_rhs-vn_lhs);hold on;
% plot(sint,vn_rhs-vn_lhs_2(sint));hold off;
% title('difference v_n eqn');
% subplot(1,5,5)
% plot(sint,vs_rhs-vs_lhs);hold on;
% plot(sint,vs_rhs-vs_lhs_2(sint));hold off;
% title('difference v_s eqn');

%% Check global force integral: (derivatives from num gradient)
% integrand_vec = X(sint).*(-sin(Psi(sint)).*TFB_lhs_2(sint)-sin(Psi(sint)).*C2(sint).*tns(sint)+cos(Psi(sint)).*NFB_lhs_2(sint)-cos(Psi(sint)).*C(sint).*tss(sint)+cos(Psi(sint)).*P+cos(Psi(sint)).*[derivatives(5,1) (cos(Psi(sint(2:end-1)))./X(sint(2:end-1))).*tns(sint(2:end-1))  derivatives(5,end)]);
% integrand_vec2 = X(sint).*(-sin(Psi(sint)).*derivatives(1,:)-sin(Psi(sint)).*C2(sint).*tns(sint)+cos(Psi(sint)).*derivatives(5,:)-cos(Psi(sint)).*C(sint).*tss(sint)+cos(Psi(sint)).*P+cos(Psi(sint)).*[derivatives(5,1) (cos(Psi(sint(2:end-1)))./X(sint(2:end-1))).*tns(sint(2:end-1))  derivatives(5,end)]);
% 
% Vint(10,end)
% 
% line1=trapz(sint,integrand_vec)
% %trapz(sint,integrand_vec2)
% 
% integrand_vec3 = X(sgrid).*(-sin(Psi(sgrid)).*TFB_lhs_2(sgrid)+gradient(cos(Psi(sgrid)),hgrid).*tns(sgrid)+cos(Psi(sgrid)).*NFB_lhs_2(sgrid)-gradient(sin(Psi(sgrid)),hgrid).*tss(sgrid)-cos(Psi(sgrid)).*C1(sgrid).*tss(sgrid)+cos(Psi(sgrid)).*P+cos(Psi(sgrid)).*[derivatives(5,1) (cos(Psi(sgrid(2:end-1)))./X(sgrid(2:end-1))).*tns(sgrid(2:end-1))  derivatives(5,end)]);
% line2=trapz(sgrid,integrand_vec3)
% 
% integrand_vec4 = X(sgrid).*(-sin(Psi(sgrid)).*TFB_lhs_2(sgrid)+gradient(cos(Psi(sgrid)),hgrid).*tns(sgrid)+cos(Psi(sgrid)).*NFB_lhs_2(sgrid)-gradient(sin(Psi(sgrid)),hgrid).*tss(sgrid)-gradient(X(sgrid),hgrid).*C1(sgrid).*tss(sgrid)+gradient(X(sgrid),hgrid).*P+gradient(X(sgrid),hgrid).*[derivatives(5,1) (cos(Psi(sgrid(2:end-1)))./X(sgrid(2:end-1))).*tns(sgrid(2:end-1))  derivatives(5,end)]);
% line3=trapz(sgrid,integrand_vec4)
% 
% integrand_vec5 = P*X(sgrid).*gradient(X(sgrid),hgrid)-X(sgrid).*gradient(sin(Psi(sgrid)).*tss(sgrid),hgrid)-gradient(X(sgrid),hgrid).*sin(Psi(sgrid)).*tss(sgrid)+X(sgrid).*gradient(cos(Psi(sgrid)).*tns(sgrid),hgrid)+gradient(X(sgrid),hgrid).*cos(Psi(sgrid)).*tns(sgrid);
% line4=trapz(sgrid,integrand_vec5)

% figure(40);
% plot(sint,Vint,'.-');
% axis tight


end

