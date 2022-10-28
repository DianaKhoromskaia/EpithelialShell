function res = bc_nematic(yleft,yright, varpar, Psi, X, L, lc) 
% residual function specifying the boundary conditions for the system of
% ODE's set in ode_nematic.m

res = [ yleft(1);
        yleft(2);
        yright(1);
        yright(2)];
end


