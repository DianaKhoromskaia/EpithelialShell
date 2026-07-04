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

        %% ----------------------------------------------------------------
        % Dense uniform grid for saving smooth geometry/profile fields
        %% ----------------------------------------------------------------
        nsave = 400;
        ssave = linspace(0, L, nsave);

	% DEBUG
	disp(['UPDATED SAVEFILE ACTIVE, nsave = ', num2str(nsave)])
	
        %% save observables
        fprintf(fileID, formatSpec, ...
            t, dt, vn(seval), P(1), X(seval), C2(seval), ...
            (Z(L)-Z(0.)), 2*pi*xintegral, tcomp, ...
            L, V, X0, nmesh, intersect, seval, P1(2), v(10,end));

        %% ----------------------------------------------------------------
        % save curvature-related quantities on dense uniform grid
        %% ----------------------------------------------------------------
        dlmwrite(filename2, ssave, '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename2, C1(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename2, C2(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename2, C(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename2, dsC1(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename2, dsC(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        %% ----------------------------------------------------------------
        % save geometry on dense uniform grid
        %% ----------------------------------------------------------------
        dlmwrite(filename3, X(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename4, Z(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        dlmwrite(filename41, Psi(ssave), '-append', ...
            'precision', '%10.9f' ,'delimiter', '\t');

        %% ----------------------------------------------------------------
        % save BVP solution quantities
        % keep these on adaptive BVP mesh svec
        %% ----------------------------------------------------------------
        if n==1

            dlmwrite(filename5, svec, ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename5, v, '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename6, vs(svec), ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename7, vn(svec), ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename12, tss(svec), ...
                'precision', '%10.9f' ,'delimiter', '\t');

            %% ------------------------------------------------------------
            % save auxiliary fields on dense uniform grid
            %% ------------------------------------------------------------
            dlmwrite(filename10, U(ssave), ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename13, s0(ssave), ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename14, ssave, ...
                'precision', '%10.9f' ,'delimiter', '\t');

            if write92 || write93
                dlmwrite(filename11, Q(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            %% ------------------------------------------------------------
            % save active profiles
            %% ------------------------------------------------------------
            if write9
                dlmwrite(filename9, zeta(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write91
                dlmwrite(filename91, zetac(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write92
                dlmwrite(filename92, zetanem(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write93
                dlmwrite(filename93, zetacnem(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write94
                dlmwrite(filename94, kappa(ssave), ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

        else

            dlmwrite(filename5, svec, '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename5, v, '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename6, vs(svec), '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename7, vn(svec), '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename12, tss(svec), '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            %% ------------------------------------------------------------
            % save auxiliary fields on dense uniform grid
            %% ------------------------------------------------------------
            dlmwrite(filename10, U(ssave), '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename13, s0(ssave), '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            dlmwrite(filename14, ssave, '-append', ...
                'precision', '%10.9f' ,'delimiter', '\t');

            if write92 || write93
                dlmwrite(filename11, Q(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            %% ------------------------------------------------------------
            % save active profiles
            %% ------------------------------------------------------------
            if write9
                dlmwrite(filename9, zeta(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write91
                dlmwrite(filename91, zetac(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write92
                dlmwrite(filename92, zetanem(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write93
                dlmwrite(filename93, zetacnem(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

            if write94
                dlmwrite(filename94, kappa(ssave), '-append', ...
                    'precision', '%10.9f' ,'delimiter', '\t');
            end

        end
end
