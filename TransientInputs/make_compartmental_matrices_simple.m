%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A,B] = make_compartmental_matrices(Ri,Rm,Cm,r,l,n)
% Outputs: Square matrix A, Diagonal matrix B, Vector of radii R, and
% Vector of lengths L
% Inputs: Ri - axial resistance, Rm - membrane resistance, Cm - membrane
% capacitance, r - compartment radius, l - compartment length, n -
% number of compartments, and parents - a vector specifying the parent of 
% each compartment.
% 
% This function constructs the matrices for the compartmental model of a
% neuron dvdt = Av + Bu where v is the state vector of voltage in each
% compartment, A is the system matrix relating the compartments to each
% other, B is a diagonal matrix contain the inverse of the capacitances of
% each compartment, and u is the input vector containing the current
% applied to each compartment.
% 
% Simon (ty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,L,C] = make_compartmental_matrices_simple(Ri,Rm,Cm,r,l,n,parents)
%Check inputs
if nargin < 7
    error('Not enough input arguments');
elseif nargin > 7
    error('Too many input arguments');
end

if isrow(l), l = l'; end;

% num_rows = sum(n);
lambda = sqrt(Rm/Ri*r/2);
L = l./lambda;
C = 2.*r.*pi.*L.*Cm;           % capacitance of membrane segments
gi = (1./Ri).*(pi.*r.^2)./L;    % axial conductance of segments
gm = (pi./Rm).*2.*r.*L;         % membrane conductance of segments

% Set up Matrices
A = zeros(n);
B = zeros(n);
for j = 1:n
    B(j,j) = 1/C(j);
    children = find(parents==j);
    if parents(j) <= 0
        gi(j) = 0;
    end
    for i = 1:n
        if ismember(i,children)
            A(j,i) = gi(i).*B(j,j);
        elseif i==parents(j)
            A(j,i) = gi(j).*B(j,j);
        end
    end
    if isempty(children)
        A(j,j) = -(gi(j)+gm(j)).*B(j,j);        
    else
        A(j,j) = -(sum(gi(children))+gi(j)+gm(j)).*B(j,j);
    end
end

end
