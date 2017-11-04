%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% t - current timepoint
% t_shift - shift in time to the left
% c - compartment number receiving input
% n - number of compartments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = make_u(t,t_shift,c,n)
u = zeros(n,1);

% decaying exp
tpeak = 50; % us
g = 0.001;

u(c) = ((t-t_shift)/tpeak)*exp(1-(t-t_shift)/tpeak)*g;

end