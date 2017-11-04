%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u] = make_u(t,t_shift,c,n)
% % step function
u = zeros(n,1);
u(c) = heaviside(t+t_shift);


% decaying exp
% tpeak = 1; % ms
% u(c) = exp(-t/tpeak);
end