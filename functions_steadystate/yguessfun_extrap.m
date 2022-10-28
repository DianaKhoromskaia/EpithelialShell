function y = yguessfun_extrap(varargin)  %(s, k, sol) % k-region
s = varargin{1};%first argument
sol = varargin{nargin}; %last argument
    
y = deval(sol, s);
end