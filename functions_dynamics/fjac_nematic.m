function [dfdv,dfdvar] = fjac_nematic(s, v, varpar, Psi, X, L, lc)
%Jacobian of RHS for the bvp4c solver

dfdv = zeros(2,2);
dfdvar = zeros(2,2);

switch s
    case 0.
        %partial derivatives wrt v-components 
        
        %partial derivatives wrt to parameters
        dfdvar(2,1) = 1;
       
    
    case L
        %partial derivatives wrt v-components
        
        %partial derivatives wrt to parameters
        dfdvar(2,2) = 1;
        
    otherwise
        
        %partial derivatives wrt v-components
        
        dfdv(1,2) = 1;
        
        dfdv(2,1) = ((-1+3*v(1,:).*v(1,:))./(2*lc^2) + 4*cos(Psi(s)).*cos(Psi(s))./(X(s).*X(s)));
        dfdv(2,2) = -cos(Psi(s))./X(s);
       
end


end

