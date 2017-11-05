%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u] = make_u_step(t,t_shift,c,n)
% Outputs: u - nx1 array of input step current
% Inputs:
%   t - time at which to evaluate function (scalar)
%   t_shift - time shift
%   c - compartment at which current is applied
%   n - number of compartments
%
% This function returns the input current at a given time t with a time
% shift t_shift at the compartment c. It is used in dvdt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_step] = make_u_step(t,t_shift,c,n)
u_step = zeros(n,1);
u_step(c) = heaviside(t-t_shift);
end