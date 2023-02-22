function dvds = ode_nematic(s, v, varpar, Psi, X, L, lc)
% rhs for ode for the vector v = (vs, dsvs, vn, dsvn, mss, dsmss, dV(s), dX(s), s0(s), U), solved starting from south pole s=0
dsw0 = varpar(1);
%dswL = varpar(2);
onevec = ones(size(s));

q = v(1,:);
dsq = v(2,:);

dvds = [    dsq;
            (1/(2*lc^2))*q.*(q.*q-onevec) + (cos(Psi(s))./X(s)).*(4*(cos(Psi(s))./X(s)).*q - dsq)];
% % 
% figure(1)
% subplot(2,1,1)
% plot(s, v(12,:), 'o-')
% %legend('s0')
% legend('tss','vn','dsvn','mss','tns','v0','x0','ds(E)','q','dsq','vs','ds(L)')
% subplot(2,1,2)
% plot(s, dvds(12,:),'o-')

% specify force balance at SP:
        indices = find(~s);
        zervec = zeros(size(s(indices)));
        onevec = ones(size(s(indices)));
        
        dvds(:,indices) =   [   zervec;
            dsw0*onevec];

 % specify force balance at NP:       
%         indices = find(~(s-L));
%         zervec = zeros(size(s(indices)));
%         onevec = ones(size(s(indices)));
%         
%         dvds(:,indices) =   [   zervec;
%             dswL*onevec]; 
        
% figure(1)
% subplot(2,1,1)
% plot(s, v, 'o-')
% % %legend('s0')
% % legend('tss','vn','dsvn','mss','tns','v0','x0','ds(E)','q','dsq','vs','ds(L)')
% subplot(2,1,2)
% plot(s, dvds,'o-')


end

