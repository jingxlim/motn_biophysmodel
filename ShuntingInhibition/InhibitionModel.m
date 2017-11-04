%%%%

neuron = load('/Users/Simon/Google Drive/Graduate School/Hopkins/Models of the Neuron/Homeworks/Project1/motn_biophysmodel/ShuntingInhibition/202-2-23nj.CNG.txt');
axon = find(neuron(:,2)==2);
neuron(axon,:) = [];

%Define parameters
Ri = 100; %Ohm - cm
Rm = 1e4; %Ohm - cm^2
Cm = 1; %microF/cm^2
Er = -70; %Mv
Ee = 60; %Mv
Ei = Er; %Mv

num_compartments = size(neuron,1);

%Convert micrometers to centimeters
neuron(:,3:6) = neuron(:,3:6)./10000;

lambda = sqrt((Rm./Ri).*(neuron(:,6)./2));

l = zeros(num_compartments,1);
for n = 1:num_compartments
    p = neuron(n,7);
    if p < 1
        cx = neuron(n,3);
        cy = neuron(n,4);
        cz = neuron(n,5);
        dx = 0-cx;
        dy = 0-cy;
        dz = 0-cz;
    else
        dx = neuron(p,3)-neuron(n,3);
        dy = neuron(p,4)-neuron(n,4);
        dz = neuron(p,5)-neuron(n,5);
    end
    l(n) = sqrt(dx^2 + dy^2 + dz^2);
end

if any(l==0)
    i0 = l==0;
    if neuron(i0,2)==1 && length(neuron(:,2)==1)>1
        l(i0) = mean(l(neuron(:,2)==1 & ~i0));
    else
        l(i0) = 0.0005;
    end
end

L = l./lambda;

disp('Is the compartment size constraint met for all compartments?');
if any(L>0.1)
    disp('No, it is not met for all compartments.');
else
    disp('Yes!');
end

%Plot neuron morphology
%plot_tree(neuron);

Ee = ones(num_compartments,1)*Ee;
Ei = ones(num_compartments,1)*Ei;
Er = ones(num_compartments,1)*Er;

n = ones(num_compartments,1);
[A,B,R,~,C] = make_compartmental_matrices_inhibition(Ri,Rm,Cm,Er,Ee,Ei,neuron(:,6),l,n,neuron(:,7));

%Set up input vector u
u_excite = zeros(2*num_compartments,1);
u_inhibit = zeros(2*num_compartments,1);
ge = zeros(num_compartments,1);
gi = zeros(num_compartments,1);

ge(27) = 0.001; % near
%ge(80) = 0.001; % far

gi(25) = 0.001; % On path
%gi(82) = 0.001; % Off path

u_excite(1:2:end) = ge;

u_inhibit(1:2:end) = ge;
u_inhibit(2:2:end) = gi;

%Set up G Matrix
G_excite = zeros(num_compartments);
G_inhibit = zeros(num_compartments);
for j = 1:num_compartments
    G_excite(j,j) = -(ge(j))/C(j);
    G_inhibit(j,j) = -(ge(j)+gi(j))/C(j);
end

t_shift_excite = 100;
t_shift_inhibit = 100;

DvpDt = @(t,vp,u,G) A*vp + B*make_u(t,t_shift_excite,t_shift_inhibit,27,25,num_compartments) + ...
    make_G(make_u(t,t_shift_excite,t_shift_inhibit,27,25,num_compartments),C)*vp;
tp = 0:1:1e4;

v0 = zeros(num_compartments,1);
[tr,Vp] = ode23(DvpDt,tp,v0); 
V = zeros(size(Vp));

Vss_excite = -inv(A+G_excite)*B*u_excite;
Vss_inhibit = -inv(A+G_inhibit)*B*u_inhibit;
x = zeros(num_compartments,1);
children = cell(num_compartments,1);
for n = 1:num_compartments
    p = neuron(n,7);
    children{n} = neuron((neuron(:,7)==n),1);
    if p > 0 && neuron(n,2)~=1
        x(n) = x(p) + l(n);
    end
    V(:,n) = Vp(:,n)*(Ee(n)-Ei(n))+Er(n);
    Vss_excite(n) = Vss_excite(n)*(Ee(n)-Ei(n))+Er(n);
    Vss_inhibit(n) = Vss_inhibit(n)*(Ee(n)-Ei(n))+Er(n);
end
X = x./lambda;

% %Plot voltage based on distance from the soma
% figure; subplot(1,2,1); grid on; hold on;
% for n = 1:num_compartments
%     parent = neuron(n,7);
%     if parent >0
%         Xp = [X(parent) X(n)];
%         Vsep = [Vss_excite(parent) Vss_excite(n)];
%         Vsip = [Vss_inhibit(parent) Vss_inhibit(n)];
%     else
%         Xp = X(n);
%         Vsep = Vss_excite(n);
%         Vsip = Vss_inhibit(n);
%     end
%     plot(Xp,Vsep,'b.-');
%     plot(Xp,Vsip,'k.-'); %For Off-path
% %     plot(Xp,Vsip,'r.-'); %For On-path
% end
% xlabel('Electrostatic Distance from the Soma');
% ylabel('Voltage [mV]');
% 
% subplot(1,2,2); grid on; hold on;
% plot(Vss_excite,'b.-');
% plot(Vss_inhibit,'k.-'); %For Off-path
% % plot(Vss_inhibit,'r.-'); %For On-path
% xlabel('Compartment Number');
% ylabel('Voltage [mV]');
% 
% suptitle('Steady State: Excitatory Input at 27, Inhibitory at 82 (Off-path, near)');
% legend('Excitation Only','Excitation + Inhibition','Location','Best');
% 
