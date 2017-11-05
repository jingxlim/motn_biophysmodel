%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 1
% Topic 3: temporal response to transient inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% file & plot (not fancy)
dend = load('whale.txt');
% plot_tree(dend);

%% comp model
t = 2; x = 3; y = 4; z = 5; r = 6; p = 7;
Ri = 100; Rm = 10000; Cm = 1; Iapp = 1e-9;
n = size(dend,1);
l = comp_len(dend);
[A,B,L] = make_compartmental_matrices(Ri,Rm,Cm,dend(:,r),l,n,dend(:,p));

%% transient inputs
% step current
t_n = 50;
t_span = linspace(1,250,t_n);
t_shift = 25;
c = 25;

% conductances
ge = zeros(n,1);
% gi = zeros(n,1);

% ge(27) = 0.001;

% % response
% u = zeros(n,t_n);
% for i = 1:t_n
%     u(:,i) = make_u(t_span(i),t_shift,c,n)*Iapp;
% end

dvdt = @(t,v) A*v + B*make_u(t,t_shift,c,n) + make_G(make_u(t,t_shift,c,n),Cm)*v;
[tp,v] = ode23(dvdt, t_span, zeros(n,1));

v = v';

% plot for various compartments
figure('Name','Models of the Neuron Project: Topic 3')
hold on
plot(t_span,v(1,:),'m')
plot(t_span,v(c,:),'g')
plot(t_span,Iapp,'c')
plot(t_span,v(c-20,:),'b')
plot(t_span,v(c+30,:),'b')
legend('Soma','Compartment w/Iapp','Iapp','Other compartments')
xlabel('Time (us)'); ylabel('Voltage (mV)');
