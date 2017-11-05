%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u] = make_u_exp(t,t_shift,c,n)
% Outputs: u - nx1 array of input decaying exponential current 
% Inputs:
%   t - time at which to evaluate function (scalar)
%   t_shift - time shift
%   c - compartment at which current is applied
%   n - number of compartments
%
% This function returns the input current at a given time t with a time
% shift t_shift at the compartment c. It is used in dvdt and in make_G.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_exp] = make_u_exp(t,t_shift,c,n)
u_exp = zeros(2*n,1);   % initialize assuming input 0
tpeak = 50;         % us
g = 0.001;          % actually made up shhh

u_exp(2*c-1) = (t-t_shift)/tpeak*exp(1-(t-t_shift)/tpeak)*g;

if any(u_exp<0)
    ii0 = u_exp<0;
    u_exp(ii0) = 0;
end

end