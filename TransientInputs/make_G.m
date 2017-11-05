%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [G] = make_g(u,Cm)
% Outputs: G - nx1 array of conductances
% Inputs:
%   u  - nx1 array of input currents
%   Cm - membrane capacitance
%
% This function returns the conductances of the compartments using the
% input current vector u. It is used in the function dvdt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G] = make_G(u,Cm)
G = diag(-(u(1:2:end)+u(2:2:end))./Cm);

end