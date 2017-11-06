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
function [u_step] = make_u_step(t,t_shift,step_dur,c,n)
u_step = zeros(n,1);

if t <= t_shift 
    u = 0;
elseif t >= t_shift+step_dur;
    u = 0;
else
    u = 1;
end

u_step(c) = u;

end