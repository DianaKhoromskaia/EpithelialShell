function savetofile(X, Z, Psi, svec, svecnew, seval, fileID, formatSpec, t, dt, P1, P, C1, C2, C, dsC1, dsC, xintegral, tcomp, L, dX0, V, X0, nmesh, v, vs, vn, tss, U, Q, s0, kappa, zeta, zetac, zetanem, zetacnem, filename2, filename3, filename4, filename41, filename5, filename6, filename7, filename9, filename91, filename92, filename93, filename94, write9, write91, write92, write93, write94, filename10, filename11, filename12, filename13, filename14, n)

        %% test for intersections:
        Xfull = X(svecnew);
        Zfull = Z(svecnew);
        inters = intersections(Xfull, Zfull);
        
        if ~isempty(inters)
            intersect = 1;
        else 
            intersect = 0;
        end
        
        fprintf(fileID, formatSpec, t, dt, vn(seval), P(1), X(seval), C2(seval), (Z(L)-Z(0.)), 2*pi*xintegral, tcomp, L, V, X0, nmesh, intersect, seval, P1(2), v(10,end));

        dlmwrite(filename2, svecnew, '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename2, C1(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename2, C2(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename2, C(svecnew), '-append','precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename2, dsC1(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename2, dsC(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename3, X(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename4, Z(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
        dlmwrite(filename41, Psi(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');

        if n==1
            dlmwrite(filename5, svec, 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename5, v, '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename6, vs(svec), 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename7, vn(svec), 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename12, tss(svec), 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename10, U(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename13, s0(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename14, svecnew, 'precision', '%10.9f' ,'delimiter', '\t');
            if write92||write93
                dlmwrite(filename11, Q(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end

            % save active profile to file:
            if write9
                dlmwrite(filename9, zeta(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write91
                dlmwrite(filename91, zetac(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write92
                dlmwrite(filename92, zetanem(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write93
                dlmwrite(filename93, zetacnem(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write94
                dlmwrite(filename93, kappa(svecnew), 'precision', '%10.9f' ,'delimiter', '\t');
            end
   
        else 
            dlmwrite(filename5, svec, '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename5, v, '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename6, vs(svec), '-append','precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename7, vn(svec), '-append','precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename12, tss(svec), '-append','precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename10, U(svecnew), '-append','precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename13, s0(svecnew), '-append','precision', '%10.9f' ,'delimiter', '\t');
            dlmwrite(filename14, svecnew, '-append','precision', '%10.9f' ,'delimiter', '\t');
            if write92||write93
                dlmwrite(filename11, Q(svecnew), '-append','precision', '%10.9f' ,'delimiter', '\t');
            end

            % save active profile to file:
            if write9
                dlmwrite(filename9, zeta(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write91
                dlmwrite(filename91, zetac(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write92
                dlmwrite(filename92, zetanem(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write93
                dlmwrite(filename93, zetacnem(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            end
            if write94
                dlmwrite(filename94, kappa(svecnew), '-append', 'precision', '%10.9f' ,'delimiter', '\t');
            end
        end
end

