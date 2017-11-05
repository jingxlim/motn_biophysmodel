%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_tree(file) 
% Outputs: none, displays figure
% Inputs: dend - matrix of dendritic trees (std format from neuromorpho)
% 
% This function plots the dendritic tree of a neuron based on the
% compartments provided, assuming there is only soma and dendritic
% components. This has been derived from HW4 Problem 3a.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_tree(dend)
% variables
n = size(dend,1);
t = 2; x = 3; y = 4; z = 5; p = 7;

figure('Name','Dendrite Tree')
hold on

% finding the soma & making it a star!
soma = find(dend(:,t)==1);
plot3(dend(soma,x),dend(soma,y),dend(soma,z),'Marker','*','Color','r')

% plotting
for i = 1:n                 % each compartment
    j = find(dend(:,p)==i);
    for n = 1:length(j)     % each child of compartment
        tempx = [dend(i,x) dend(j(n),x)];
        tempy = [dend(i,y) dend(j(n),y)];
        tempz = [dend(i,z) dend(j(n),z)];
        plot3(tempx,tempy,tempz,'Color','g')
    end
end

title('Dendridic Tree', 'fontsize', 26)
grid on; view(3); rotate3d on;
xlabel('x (cm)', 'fontsize', 14); ylabel('y (cm)', 'fontsize', 14);
zlabel('z (cm)', 'fontsize', 14);

end