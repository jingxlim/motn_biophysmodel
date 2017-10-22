%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 1
% Topic 5: low pass filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Refresh workspace
clr;

%% Download file
new = 0;

url = 'http://neuromorpho.org/dableFiles/cameron/CNG%20version/HP72N6B.CNG.swc';
name_split = strsplit(url, '/');
origname = name_split(end);
filename = char(strrep(origname, 'swc', 'txt'));
if new    outfilename = websave(filename,url); end
%% Import and preprocess data
dat = table2array(readtable(filename));  % import file and convert to array
fdata = dat(:,all(~isnan(dat)));  % remove columns with nans
data = fdata(fdata(:,2)~=2,:);  % remove axons
datafile = strrep(filename,'txt', 'mat');
if new  save(datafile,'fdata','data'); end% save data to file
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
ratio = length./lambda;
constraint =  ratio > 0.1;

if find(constraint)
    disp('Constraint is NOT satisfied for all dendritic segments')
    disp('Compartment '+string(find(constraint))+' does not satisfy the constraint')
else
    disp('Constraint is satisfied for all dendritic segments!')
end

%% Initialize compartmental model: construct matrix A

N = size(cmprt,1);  % number of compartments

% extrinsic values for compartment j
C = @(r,dl) 2*pi*r*dl*Cm;
gi = @(r,dl) (pi*r^2)/(dl*Ri);
gm = @(r,dl) (2*pi*r*dl)/Rm;

% initialize and construct matrix A
A_ = zeros(N);

for i=1:numel(cmprt)
    % find daughter branches
    ccmprt = cmprt(i);  % only 1
    daughter{i} = find(parent == ccmprt);
    
    % find daughter and parent compartments
    dcmprts = daughter{i};  % might have more than 1
    pcmprt = parent(i);  % only 1
    
    % populate A matrix
    A_(ccmprt,ccmprt) = A_(ccmprt,ccmprt)...
                              -gm(radius(ccmprt),length(ccmprt));
    if pcmprt~=-1
        A_(ccmprt,ccmprt) = A_(ccmprt,ccmprt)...
                                  -gi(radius(ccmprt),length(ccmprt));
        A_(ccmprt,pcmprt) = A_(ccmprt,pcmprt)...
                                  +gi(radius(ccmprt),length(ccmprt));
    end
                              
    for i=1:numel(dcmprts)
        dcmprt = dcmprts(i);
        A_(ccmprt,ccmprt) = A_(ccmprt,ccmprt)...
                                  - gi(radius(dcmprt),length(dcmprt));
        A_(ccmprt,dcmprt) = A_(ccmprt,dcmprt)...
                                  + gi(radius(dcmprt),length(dcmprt));
    end
end

Cm_all = zeros(N,1);

% calculate the membrane capacitance for each compartment
for i=1:numel(cmprt)
    ccmprt = cmprt(i);
    Cm_all(i) = C(radius(ccmprt),length(ccmprt));
end

Cm_matrix = repmat(Cm_all,1,N);
A = A_ ./ Cm_matrix;

%% Initialize comparmental model: construct matrices B and U
B = eye(N).*(1./Cm_matrix);

% construct U vector
Iapp = 1e-9;  % mA
inj_cmprt = 348;
U = zeros(N,1);
U(inj_cmprt,1) = Iapp;

%% Find steady state voltages

V = - inv(A)*B*U;

lambda = sqrt((Rm*radius)/(Ri*2));
L = length./lambda;

V_ = transpose(V);  % column vector of ss voltages over compartments

%% Find location of compartments and plot steady state
% find branch points
u_parent = unique(parent);
hist_parent = histc(parent(:),u_parent);
branch_pts = u_parent(hist_parent>1);

% find daughter branches
for i=1:numel(branch_pts)
    pbranch = branch_pts(i);
    dbranch = find(parent==pbranch);
    dbranches{pbranch} = dbranch;
end

branches={};

% find branches
for i=1:numel(branch_pts)
    pbranch = branch_pts(i)
    for j=1:numel(dbranches{pbranch})
        dbranch = dbranches{pbranch}(j)
        branches{size(branches,2)+1} = [pbranch,dbranch];
        if ismember(dbranch,branch_pts)
            continue
        else
            while ~ismember(dbranch,branch_pts)
                dbranch = find(parent==dbranch)
                branches{size(branches,2)} = [branches{size(branches,2)},dbranch];
            end
        end
    end
end
 
position = {};
X = zeros(size(cmprt));
X(1) = 0;
for i=1:numel(branches)
    counter = X(branches{i}(1));
    position{size(position,2)+1} = [counter];
    for j=2:numel(branches{i})
        comp = branches{i}(j);
        counter = counter + L(comp);
        position{size(position,2)} = [position{size(position,2)},counter];
        X(comp)= counter;
    end
end

figure(2); clf; hold on;

for i=1:numel(branches)
    plot(position{i},V_(branches{i}),'-');
    % text(position{i},V_(branches{i})+0.001, cellstr(string(branches{i})),'FontSize',8);
end

xlabel('Electrontonic distance from soma'); ylabel('Steady state voltage [mV]');
title('Steady-state voltage along cables: Injection into branch '+string(inj_cmprt));
saveas(gcf, 'ss_voltage.png');