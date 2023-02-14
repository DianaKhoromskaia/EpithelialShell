function plotting(svecuni, sol, t, told, P, P0, V, V0, R0, X0, X, Z, Xinit, Zinit, sinit)


set(0,'defaultfigurecolor',[1 1 1])

    fig6 = figure(6);
    fig6.Position = [0.100 100.00 800. 400.];
    subplot(2,4,1)
    plot(sol.x, sol.y(2,:))
    legend('v_n(s)')
    subplot(2,4,2)
    hold on
    plot(told, P,'r.', told, P0,'k.');
    title('P')
    subplot(2,4,5)
    hold on
    %plot(told, max(abs(sol.y(2,:))),'b.');
    %plot(told, sol.y(2,1),'b.');
    plot(told, V/V0,'r.');%, told, V0,'k.');
    title('V')
    subplot(2,4,6)
    hold on;
    %plot(sinit, Uinit(sinit), sol.x, U(sol.x))
    %plot(sol.x, U(sol.x))
    plot(told, sol.y(10,end),'k.');
    %plot(sol.x, sol.y(10,:),'r.')
    title('I(L)')
%     plot(told, abs(P1(2)/P),'r.')
%     title('|f_c/P|')
    %plot(sol.x, zetanem(sol.x))
    %legend('\zeta_n q')%('u_k^k')
    axis tight
    subplot(1,2,2)
    %clf
    % plot original shape:
    plot(Xinit(sinit), Zinit(sinit), 'b'); hold on
    plot(-Xinit(sinit), Zinit(sinit), 'b')
    
    % plot shape at time t:
%     plotshape3D_2(Xfull,Zfull,svecuni,zeta, gamma+A)
    plot(X(svecuni), Z(svecuni), 'r','LineWidth',2);
    plot(-X(svecuni), Z(svecuni), 'r','LineWidth',2); hold off;
    %plot(0,X0, 'ob'); 
    axis square
    axis([-2.5*R0 2.5*R0 X0-2.5*R0 X0+2.5*R0])
     title(strcat('t=',num2str(t)))
    %saveas(figure(6),'currentshape','png');
    saveas(figure(6),'currentshape','fig');
    drawnow;

end

