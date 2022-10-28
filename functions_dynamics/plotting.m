function plotting(svecuni, C1, C2, C, dsC1, dsC, sol, t, told, P1, P, P0, V, V0, R0, X0, U, zetanem, X, Z, Uinit, Xinit, Zinit, sinit, vidObj1, FixedPar)

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
    plot(told, V/V0,'r.');
    title('V')
    subplot(2,4,6)
    hold on;
    plot(told, sol.y(10,end),'k.');
    title('I(L)')
    axis tight
    subplot(1,2,2)
    %clf
    % plot original shape:
    plot(Xinit(sinit), Zinit(sinit), 'b'); hold on
    plot(-Xinit(sinit), Zinit(sinit), 'b')
    % plot shape at time t:
    plot(X(svecuni), Z(svecuni), 'r','LineWidth',2);
    plot(-X(svecuni), Z(svecuni), 'r','LineWidth',2); hold off;
    axis square
    axis([-2.5*R0 2.5*R0 X0-2.5*R0 X0+2.5*R0])
     title(strcat('t=',num2str(t)))
    %writeVideo(vidObj1, getframe(gcf));
    saveas(figure(6),'currentshape','png');
    saveas(figure(6),'currentshape','fig');
    drawnow;
    
end

