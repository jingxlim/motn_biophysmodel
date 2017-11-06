%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% t - current timepoint
% t_shift - shift in time to the left
% ce - compartment number receiving excitatory input
% ci - compartment number receiving inhibitory input
% n - number of compartments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = make_u(t,t_shift_excite,t_shift_inhibit,ce,ci,n)
u = zeros(2*n,1);

% decaying exp
tpeak = 50; % us
g = 0.001;

u(2*ce-1) = ((t-t_shift_excite)/tpeak)*exp(1-(t-t_shift_excite)/tpeak)*g;
if ~isempty(t_shift_excite) && ~isempty(ci)
    u(2*ci) = ((t-t_shift_inhibit)/tpeak)*exp(1-(t-t_shift_inhibit)/tpeak)*g;
end

if any(u<0)
   ii0 = u<0;
   u(ii0) = 0;
end

end