%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [G] = make_g(t,t_shift,c,n)
% Outputs: u - nx1 array of input current
% Inputs:
%   t - time at which to evaluate function (scalar)
%   t_shift - time shift
%   c - compartment at which current is applied
%   n - number of compartments
%
% This function returns the input current at a given time t with a time
% shift t_shift at the compartment c. It is used within a loop to 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G] = make_G(u,Cm)
% G = diag(-(ge+gi)./Cm);
G = diag(-(u(1:2:end)+u(2:2:end))./Cm);

end