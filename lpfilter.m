%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 1
% Topic 5: low pass filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Download file
url = 'http://neuromorpho.org/dableFiles/cameron/CNG%20version/71INTER.CNG.swc';
name_split = strsplit(url, '/');
origname = name_split(end);
filename = char(strrep(origname, 'swc', 'txt'));
outfilename = websave(filename,url);

%% Import and preprocess data
dat = table2array(readtable(filename));  % import file and convert to array
fdata = dat(:,all(~isnan(dat)));  % remove columns with nans
data = fdata(fdata(:,2)~=2,:);  % remove axons
datafile = strrep(filename,'txt', 'mat');
save(datafile,'fdata','data'); % save data to file
%% Plot dendritic tree
load(datafile) % variable name is data

cmprt = data(:,1);
type = data(:,2);  % 1-soma, 2-axon, 3-dendrite
x = data(:,3) * 1e-4;  % cm
y = data(:,4) * 1e-4;  % cm
z = data(:,5) * 1e-4;  % cm
radius = data(:,6) * 1e-4;  % cm
parent = data(:,7);

figure(1); clf; hold on
co = get(gca,'colororder');

for i=1:numel(cmprt)
    compartment = i;
    
    % coordinates of compartment
    x_cmprt = x(i);
    y_cmprt = y(i);
    z_cmprt = z(i);
    
    % text(x_cmprt, y_cmprt, z_cmprt, num2str(i), 'FontSize',1);
    
    parent_cmprt = parent(i);
    
    % ignore soma because it is not a daughter segment
    if parent_cmprt == -1
        plot3(x_cmprt,y_cmprt,z_cmprt, '*', 'Color', co(end,:));
        continue
    end
    
    % parent coordinates
    x_parent = x(parent_cmprt);
    y_parent = y(parent_cmprt);
    z_parent = z(parent_cmprt);
    
    % join the compartment to its parent
    plot3([x_cmprt,x_parent],...
          [y_cmprt,y_parent],...
          [z_cmprt,z_parent],...
          'Color', co(1,:))
end
xlabel('x [cm]');ylabel('y [cm]');zlabel('z [cm]');
title('Dendritic tree morphology');
view(3); rotate3d on;
saveas(gcf, 'morpho.png')

%% check if the compartmental size constraint is fulfilled

% define intrinsic parameters
Ri = 100;  % ohm cm
Rm = 10000;  % ohm cm^2
Cm = 1;  % uF/cm^2
Er = 0;  % mV

% find length of compartments
length = zeros(size(cmprt));

for i=1:numel(cmprt)
    
    % coordinates of compartment
    x_cmprt = x(i);
    y_cmprt = y(i);
    z_cmprt = z(i);
    
    parent_cmprt = parent(i);
    
    % ignore soma because it is not a daughter segment
    if parent_cmprt == -1
        length(i) = 0.0005;  % cm; soma length=0.0005cm
        continue
    end
    
    % parent coordinates
    x_parent = x(parent_cmprt);
    y_parent = y(parent_cmprt);
    z_parent = z(parent_cmprt);
    
    length(i) = pdist([x_cmprt,y_cmprt,z_cmprt;...
                       x_parent,y_parent,z_parent],...
                       'euclidean');
end

% find space constant of compartments
lambda = sqrt((Rm*radius)/(Ri*2));

% check if size constraint is satisfied for all dendritic segments
constraint = length./lambda > 0.1;

if find(constraint)
    disp('Constraint is not satisfied for all dendritic segments')
else
    disp('Constraint is satisfied for all dendritic segments!')
end