%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [len] = comp_len(dend)
% Outputs: len - 1xcount array of the length of each compartment. Soma
%           compartments calculated to have length of 0
% Inputs: dend - matrix of dendritic trees (std format of neuromorpho)
%
% This function returns the absolute length of each compartment in the
% dendritic tree. It calls the function comp_len_xyz to get the x,y,z
% components of length and uses the sqrt function to calculate absolute
% length from these components.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l] = comp_len(dend)
% get the length of each component
l_xyz = comp_len_xyz(dend);
% magnitude
l = sqrt(l_xyz(:,1).^2+l_xyz(:,2).^2+l_xyz(:,3).^2);
% soma with zero value
soma = dend(:,2)==1;
z = l==0;
l(soma & z) = sum(l(soma))/sum(soma);

l = l';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [l_xyz] = comp_len_xyz(dend)
% Outputs: l_xyz - 3xcount matrix of x,y,z components of length of each
%           compartment. Soma compartment calculated to have length of 0.
% Inputs: dend - matrix of dendritic trees (std format of neuromorpho)
%
% This function calculates the absolute length of the x,y,z components of
% the dendrites. It is mainly called by the comp_len function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l_xyz] = comp_len_xyz(dend)
% variables
total = size(dend,1);
l_xyz = zeros(total,3);
t = 2; x = 3; y = 4; z = 5; p = 7;

% calculate length of compartments using parents
for i = 1:total                 % for each compartment
    j = find(dend(:,p)==i);     % find its children
    for n = 1:length(j)         % calculate len of comp using parent (comp)
        tempx = [dend(i,x) dend(j(n),x)];
        tempy = [dend(i,y) dend(j(n),y)];
        tempz = [dend(i,z) dend(j(n),z)];
        l_xyz(j(n),:) = [tempx(2)-tempx(1) tempy(2)-tempy(1) tempz(2)-tempz(1)];
    end
end

% calculate length of somas based on distance from origin
somas = find(dend(:,t)==1);
for k = 1:length(somas)
    soma = somas(k);
    l_xyz(soma,:) = dend(soma,x:z);
end
l_xyz = abs(l_xyz);
end