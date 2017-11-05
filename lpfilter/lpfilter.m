%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 1
% Topic 5: low pass filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Refresh workspace
clr;

%% Global settings
% Set up ODE waitbar
options = odeset('OutputFcn',@odewbar);  % ODE wait bar
% options = odeset('OutputFcn',@odeprog,'Events',@odeabort);  % ODE Progress Bar and Interrupt

% Set up simulation parameters
end_time = 5e3;
sim_time = 0:1e0:end_time;  % us

% Initialise simulation
outdate = datestr(datetime('now'),'yyyyMMdd_HHmmss');

%% Download file
new = 0;

url = 'http://neuromorpho.org/dableFiles/jacobs/CNG%20version/202-2-23nj.CNG.swc';
name_split = strsplit(url, '/');
origname = name_split(end);
filename = char(strrep(origname, 'swc', 'txt'));
if new    outfilename = websave(filename,url); end

%% Import and preprocess data
dat = table2array(readtable(filename));  % import file and convert to array
fdata = dat(:,all(~isnan(dat)));  % remove columns with nans
data = fdata(fdata(:,2)~=2,:);  % remove axons
% data = [data(1,:);data(data(2:end,2)~=1,:)];  % remove extra somas
datafile = strrep(filename,'txt', 'mat');
if new  save(datafile,'fdata','data'); end  % save data to file

%% Load data and plot dendritic tree
load(datafile) % variable name is data

cmprt = data(:,1);
type = data(:,2);  % 1-soma, 2-axon, 3-basal dendrite, 4-apical dendrite
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
    
    % text(x_cmprt, y_cmprt, z_cmprt, num2str(i), 'FontSize',5);
    
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
E_rest = -65;  % mV

% % measurements by Núñez-Abades et al. (1993)
% tau = 6.8e-3  % ms
% rin = 39.2e6  % Mohm
% areas = [16851.7,25188.4,23732.2,31959.9,25699.3,24629.3,31242.6]*1e-8;  % um^2 converted to cm^2
% area = mean(areas);
% lens = [5138.46,6877.29,6771.46,8142.56,6879.21,4311.81,6297.64]*1e-4;  % um converted to cm
% len = mean(lens);  % mean length across 7 cells
% dia = [0.9,1.04,1.11,1.1,1.07,1.62,1.39]*1e-4;  % um coverted to cm
% rad = mean(dia)/2;
% % assume Cm
% Cm = 1e-6;  % uF/cm^2
% % calculate Rm and Ri
% Rm = tau/Cm;  % ohm cm^2
% rm = Rm/area;  % ohm
% ri = (4*rin^2)/rm;  % ohm
% Ri = (ri*pi*rad^2)/len;  % ohm cm
% Cm = 1;  % revert back to uF/cm^2

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
C = @(r,dl) 2*pi*r*dl*Cm;  % uF
gi = @(r,dl) (pi*r^2)/(dl*Ri);  % ohm^-1
gm = @(r,dl) (2*pi*r*dl)/Rm;  % ohm^-1

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

%% Response to current step: construct matrix B
B_mat = eye(N).*(1./Cm_matrix);

%% Find location of compartments
lambda = sqrt((Rm*radius)/(Ri*2));
L = length./lambda;

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
    pbranch = branch_pts(i);
    for j=1:numel(dbranches{pbranch})
        dbranch = dbranches{pbranch}(j);
        branches{size(branches,2)+1} = [pbranch,dbranch];
        if ismember(dbranch,branch_pts)
            continue
        else
            while ~ismember(dbranch,branch_pts)
                dbranch = find(parent==dbranch);
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
    plot(position{i},branches{i},'.-');
end
xlabel('Electrontonic distance from soma'); ylabel('Compartment number');
saveas(gcf, 'branching.png');

%% Response to current step: construct matrix U
% create Iapp input current
Iapp = 1e-9;  % mA

mu = Iapp;  % mean
sigma = 0.8;  % standard deviation

seed = 0;
rng(seed)  % set random generator seed to 0
Iapp_rand = mu+sigma.*rand(1,size(sim_time,2));  % to plot

% check changes in conductance in U matrix
figure(3); clf; hold on;
plot(sim_time, Iapp_rand)
ylabel('I_{app} [mA]'); xlabel('Time [us]');

figure(4); clf; hold on;
title('Time evolution of V(X,T) in response to I_{app}');
xlabel('T'); ylabel('V [mV]'); legend('show');

figure(5); clf; hold on;
subplot(1,2,1);
xlabel('Frequency (Hz)')
ylabel('Power/Frequency difference (dB/Hz)')

figure(5); clf; hold on;
subplot(1,2,2);
xlabel('Frequency (Hz)')
ylabel('P/F diff (dB/Hz)')
xlim([0 50]);

injects = linspace(142,173,5);

for j=1:numel(injects)
    inject = injects(j);

    % construct U vector
    U_mat = @(t) makeU(t,inject,N,mu,sigma,seed);
    
    for i=1:numel(sim_time)
        t = sim_time(i);
        mat_U = U_mat(t);
        figure(3); plot(t, mat_U(inject,1),'.');
    end

    % Solve for voltage over time
    % system of differential equations
    dvdt = @(t,v) A*v + B_mat*U_mat(t);

    sim_T = sim_time/(Rm*Cm);

    [t,v] = ode23(@(t,v) dvdt(t,v), sim_time, zeros(N,1), options);

    figure(4);
    plot(sim_T,v(:,inject),'DisplayName','V_{injection}');
    plot(sim_T,v(:,1),'DisplayName','V_{soma}');

    % spectral analysis for I_app input
    fs = size(sim_time,2);
    psa = [inject,1];    
    
    figure(5+j); clf; hold on;
    set(gcf,'units','points','position',[100,100,1000,400])
    
    figure(5+j); subplot(1,2,1); hold on;
    x = Iapp_rand';
    n = size(x,1);
    xdft = fft(x);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/n:fs/2;
    plot(freq,10*log10(psdx),'DisplayName', 'Input current')
    xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)')    

    figure(5+j); subplot(1,2,2); hold on;
    powers = {};
    for i=1:numel(psa)
        ccmprt = psa(i);
        x = v(:,ccmprt);
        n = size(x,1);
        xdft = fft(x);
        xdft = xdft(1:n/2+1);
        psdx = (1/(fs*n)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:fs/n:fs/2;
        power = 10*log10(psdx);
        powers{end+1} = power;
        plot(freq, power, 'DisplayName', 'cmprt '+string(ccmprt))
    end
    xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
    grid on; legend('show')
    saveas(gcf, strcat('spectra_',num2str(inject),'_Iapp',outdate,'.png'));

    figure(5); subplot(1,2,1);
    plot(freq, powers{2}-powers{1});

    figure(5); subplot(1,2,2);
    plot(freq, powers{2}-powers{1});
    
end

figure(3); saveas(gcf, strcat('Iapp',outdate,'.png'));

figure(4);
saveas(gcf, strcat('output_Iapp',outdate,'.png'));

figure(5); suptitle('Periodogram Using FFT')
saveas(gcf, strcat('spectradiff_Iapp',outdate,'.png'));

%% Find steady state voltages in response to current step

% V = - inv(A)*B*U;
% V_ = transpose(V);  % column vector of ss voltages over compartments

%% Plot steady state
% figure(6); clf; hold on;
% 
% for i=1:numel(branches)
%     plot(position{i},V_(branches{i}),'-');
%     % text(position{i},V_(branches{i})+0.001, cellstr(string(branches{i})),'FontSize',8);
% end
% 
% xlabel('Electrontonic distance from soma'); ylabel('Steady state voltage [mV]');
% title('Steady-state voltage along cables: Injection into branch '+string(inj_cmprt));
% saveas(gcf, 'ss_voltage.png');

%% Set injection compartment
% simulations take a long time so only set 1
inj_cmprt = 173;

%% Reponse to synaptic input: construct matrix B
E_AMPA = 0;  % mV
E_GABA = -75;  % mV

B_ = zeros(N, 2*N);
for i=1:numel(cmprt)
    ccmprt = cmprt(i);
    B_(ccmprt, 2*ccmprt-1) = (E_AMPA - E_rest)/(E_AMPA - E_GABA);
    B_(ccmprt, 2*ccmprt) = (E_GABA - E_rest)/(E_AMPA - E_GABA);
end

B = B_ ./ repmat(Cm_matrix,1,2);

%% Tune synaptic input
tp_AMPA = 0.05e3;  % ms -> us
tp_GABA = 1e3;  % ms -> us
Gs_AMPA = 40e-12;  % pS -> S
Gs_GABA = 400e-12;  % pS -> S (10x higher)
D_AMPA = 1e8;  % channel/cm^2; check
D_GABA = 1e8;  % channel/cm^2; check

Gp_AMPA = @(r,dl) Gs_AMPA * D_AMPA * 2*pi*r*dl;
Gp_GABA = @(r,dl) Gs_GABA * D_GABA * 2*pi*r*dl;

g_AMPA_ = @(t,r,dl) (t/tp_AMPA).*exp(1-(t/tp_AMPA)).*Gp_AMPA(r,dl);
g_GABA_ = @(t,r,dl) (t/tp_GABA).*exp(1-(t/tp_GABA)).*Gp_GABA(r,dl);

figure(50); clf;

% Plot AMPA and GABA synaptic inputs on separate axes
subplot(2,1,2); hold on;
yyaxis left
plot(sim_time, g_AMPA_(sim_time, radius(inj_cmprt), length(inj_cmprt)),...
     'DisplayName','AMPA');
ylabel('AMPA conductance [S]')
yyaxis right
plot(sim_time, g_GABA_(sim_time, radius(inj_cmprt), length(inj_cmprt)),...
     'DisplayName','GABA');
ylabel('GABA conductance [S]'); xlabel('Time [us]');
legend('show'); title('Normalized synaptic inpuits');

% Plot AMPA and GABA synaptic inputs on the same axis
subplot(2,1,1); hold on
plot(sim_time, g_AMPA_(sim_time, radius(inj_cmprt), length(inj_cmprt)),...
     'DisplayName','AMPA');
plot(sim_time, g_GABA_(sim_time, radius(inj_cmprt), length(inj_cmprt)),...
     'DisplayName','GABA');
ylabel('Conductance [S]'); xlabel('Time [us]');
legend('show'); title('Synaptic inputs')

%% Create synaptic input trains

% Specify input times
AMPA_inputt = rand(1,1000) .* end_time;
% AMPA_inputt = [0] .* 1e2;
GABA_inputt = [] .* 1e-5;

g_AMPA = @(t,inputt,r,dl) ((t-inputt)/tp_AMPA).*exp(1-((t-inputt)/tp_AMPA))...
           .*Gp_AMPA(r,dl);
g_GABA = @(t,inputt,r,dl) ((t-inputt)/tp_GABA).*exp(1-((t-inputt)/tp_GABA))...
           .*Gp_GABA(r,dl);
       
AMPA_cond = zeros(size(sim_time));

for i=1:numel(AMPA_inputt)
    inputt = AMPA_inputt(i);
    cond = @(t) subplus(g_AMPA(t,inputt,radius(inj_cmprt),length(inj_cmprt)));
    AMPA_cond = AMPA_cond + cond(sim_time);
end

figure(51); clf; hold on;
plot(sim_time, AMPA_cond);
plot([AMPA_inputt;AMPA_inputt], [zeros(size(AMPA_inputt));ones(size(AMPA_inputt))*1e-9], 'k-');
set(gca,'TickDir','out') % draw the tick marks on the outside
ylabel('Conductance [S]'); xlabel('Time [us]');
set(gcf,'units','points','position',[100,100,1000,400])
saveas(gcf, strcat('input',outdate,'.png'));

%% Insert synapses: construct matrix G(t)
G_ = @(t) make_G_(t,inj_cmprt,N,radius,length,tp_AMPA,tp_GABA,Gs_AMPA,Gs_GABA,D_AMPA,D_GABA,AMPA_inputt,GABA_inputt);
G = @(t) G_(t) ./ Cm_matrix;

% check changes in conductance in G matrix
figure(52); clf; hold on;
for i=1:numel(sim_time)
    t = sim_time(i);
    G_mat = G(t);
    plot(t, G_mat(inj_cmprt,inj_cmprt),'.');
end

ylabel('Conductance changes [S]'); xlabel('Time [us]');
set(gcf,'units','points','position',[100,100,1000,400])
saveas(gcf, strcat('g',outdate,'.png'));
%% Response to synaptic input: construct matrix U(t)
U = @(t) make_U(t,inj_cmprt,N,radius,length,tp_AMPA,tp_GABA,Gs_AMPA,Gs_GABA,D_AMPA,D_GABA,AMPA_inputt,GABA_inputt);

% check synaptic input train is properly constructed in U
figure(53); clf; hold on;
for i=1:numel(sim_time)
    t = sim_time(i);
    U_mat = U(t);
    plot(t, U_mat(2*inj_cmprt-1,1),'.');
end

ylabel('Conductance [S]'); xlabel('Time [us]');
set(gcf,'units','points','position',[100,100,1000,400])
saveas(gcf,strcat('u',outdate,'.png'));
%% Solve for voltage over time

% system of differential equations
dVdt = @(t,V) A*V + B*U(t) + G(t)*V;

sim_T = sim_time/(Rm*Cm);

[t,V] = ode23(@(t,V) dVdt(t,V), sim_time, zeros(N,1),options);

figure(54); clf; hold on;
plot(sim_T,V(:,inj_cmprt),'DisplayName','V_{injection}');
plot(sim_T,V(:,1),'DisplayName','V_{soma}');
title('Time evolution of V(X,T)');
xlabel('T'); ylabel('V [mV]'); legend('show');
saveas(gcf, strcat('output_syn',outdate,'.png'));

%% spectral analysis
fs = size(sim_time,2);
psa = [inj_cmprt,1];
figure(55); clf; hold on;

x = AMPA_cond';
n = size(x,1);
xdft = fft(x);
xdft = xdft(1:n/2+1);
psdx = (1/(fs*n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/n:fs/2;
subplot(1,2,1)
plot(freq,10*log10(psdx),'DisplayName', 'Input current')

xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

subplot(1,2,2); hold on;
for i=1:numel(psa)
    ccmprt = psa(i);
    x = v(:,ccmprt);
    n = size(x,1);
    xdft = fft(x);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/n:fs/2;
    plot(freq,10*log10(psdx),'DisplayName', 'cmprt '+string(ccmprt))
end

grid on
suptitle('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
legend('show')
set(gcf,'units','points','position',[100,100,1000,400])
saveas(gcf, strcat('spectra_syn',outdate,'.png'));

%% save data
outfile = strcat('data_',outdate);
save(outfile, 'sim_time', 'Iapp_rand', 'v', 'AMPA_inputt', 'GABA_inputt', 'AMPA_cond', 'V');