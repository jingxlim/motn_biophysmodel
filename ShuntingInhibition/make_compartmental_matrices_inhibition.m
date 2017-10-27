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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,R,L] = make_compartmental_matrices_inhibition(Ri,Rm,Cm,Er,Ee,Ei,r,l,n,parents)
%Check inputs
if nargin < 10
    error('Not enough input arguments');
elseif nargin > 10
    error('Too many input arguments');
end

num_rows = sum(n);
N = cumsum(n);

%Fix length and radius vectors based on n
R = zeros(num_rows,1);
L = zeros(num_rows,1);
ni = 1;
for k = 1:num_rows
    if k <= N(ni)
        R(k) = r(ni);
        L(k) = l(ni)/n(ni);
    else
        ni = ni+1;
        R(k) = r(ni);
        L(k) = l(ni)/n(ni);
    end
end

c = 2.*R.*pi.*L.*Cm; %capacitance of membrane segments
gi = (1./Ri).*(pi.*R.^2)./L; %axial conductance of segments
gm = (pi./Rm).*2.*R.*L; %membrane conductance of segments

% Set up Matrices
A = zeros(num_rows);
B = zeros(num_rows,2*num_rows);
for j = 1:num_rows
    B(j,2*j-1) = 1/c(j)*((Ee(j)-Er(j))/(Ee(j)-Ei(j)));
    B(j,2*j) = 1/c(j)*((Ei(j)-Er(j))/(Ee(j)-Ei(j)));
    children = find(parents==j);
    if parents(j) <= 0
        gi(j) = 0;
    end
    for i = 1:num_rows
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
