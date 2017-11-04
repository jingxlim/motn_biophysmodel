%make_G.m
% Function that makes a time-varying version of the G matrix based on the u
% vector.
%
% Inputs:
% u - input vector
% c - vector of capacitances
% 
% Output: Matrix G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = make_G(u,c)

G = zeros(length(c));
for i = 1:length(c)
   G(i,i) = -(u(2*i-1)+u(2*i))/c(i);
end

end
