function [C1, C2, C, dsC1, dsC2, dsC, Psi, X, Z, X0, xintegral, sgrid, svec1, U, Q, dsQ, fs, s0] = initialisesteadystate(init_func, init_par, npoints)

L0 = pi;
L = pi*init_par(8);
sgrid = L*init_func(1,:);
%uniform grid for gradient of curvature:
svec1 = linspace(0, L, npoints);
sone = ones(1, npoints);
h = L/(npoints-1);

c1vec = init_func(6,:);   %C_{\phi}^{\phi} = sin(psi)/x
c2vec = init_func(7,:);          %C_s^s = \partial_s psi
psivec = init_func(8,:);
xvec = init_func(2,:);
zvec = init_func(3,:);
uvec = init_func(13,:);
Cvec = init_func(15,:);
fsvec = init_func(16,:);
s0vec = L0*init_func(9,:);

C1 = griddedInterpolant(sgrid, c1vec, 'spline');
C2 = griddedInterpolant(sgrid, c2vec, 'spline');
C =  griddedInterpolant(sgrid, Cvec, 'spline');
dsc1vec = gradient(C1(svec1),h);
dsc2vec = gradient(C2(svec1),h);
dsCvec = gradient(C(svec1),h);
dsC1 = griddedInterpolant(svec1, [0. dsc1vec(2:end-1) 0.],'spline');
dsC2 = griddedInterpolant(svec1, [0. dsc2vec(2:end-1) 0.],'spline');
dsC = griddedInterpolant(svec1, [0. dsCvec(2:end-1) 0.],'spline');
Psi = griddedInterpolant(sgrid, psivec, 'spline');
X = griddedInterpolant(sgrid, xvec, 'spline');
Z = griddedInterpolant(sgrid, zvec, 'spline');
U = griddedInterpolant(sgrid, uvec, 'spline');
Q = griddedInterpolant(sgrid, init_func(10,:), 'spline');
dsQ = griddedInterpolant(sgrid, init_func(11,:), 'spline');
fs = griddedInterpolant(sgrid, fsvec, 'spline');
s0 = griddedInterpolant(sgrid, s0vec,'spline');

%z-coordinate of undeformed shape:
%x_0 = [0. 0. z0];
%[s,x] = ode45(@(s,x) shape(s, x, C2), [0. L0], x_0, optode2);

 figure(5)
 subplot(1,2,1)
 plot(xvec, zvec)
 subplot(1,2,2)
 plot(sgrid, Q(sgrid),'.', svec1, dsQ(svec1),'.')
% clf
% hold on
% % plot(x(:,2), x(:,3)-x(length(x),3)/2, 'b')
% % plot(-x(:,2),x(:,3)-x(length(x),3)/2, 'b')
% plot(x(:,2), x(:,3), 'b')
% plot(-x(:,2),x(:,3), 'b')
% hold off
% axis square
% axis(L0*[-1. 1 -1 1])

%% calculate centre of shape:
%R0 = L/pi;
%xintegral = trapz(svec1, X(svec1));
xintegral = integral(@(s) X(s), 0., L);
%xintergral in units of 2\pi
%xintegral = 2*R0^2; %2*trapz(svec1, X(svec1)) 
%X0 = z0+R0; 
%X0 = trapz(svec1, X(svec1).*Z(svec1))/xintegral;
X0 = integral(@(s) X(s).*Z(s), 0., L)/xintegral;
%xintegral = trapz(s, x(:,2));
end

