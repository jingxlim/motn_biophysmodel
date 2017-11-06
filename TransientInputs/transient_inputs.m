%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 1
% Topic 3: temporal response to transient inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% file & plot (not fancy)
dend = load('whale.txt');

%Remove axon compartments
isaxon = dend(:,2)==2;
dend(isaxon,:) = [];

%Convert micrometers to centimeters
dend(:,3:6) = dend(:,3:6)./1e4;

% plot_tree(dend);

%% Jing?
jing = 1
if jing
    options = odeset('OutputFcn',@odewbar);
else
    options = odeset();
end

%% variables
% comp model
t = 2; x = 3; y = 4; z = 5; r = 6; p = 7;
Ri = 100; Rm = 10000; Cm = 1; Iapp = 1e-7;
n = size(dend,1);
l = comp_len(dend);

% transient inputs
t_n = 100;
t_span = linspace(1,5e4,t_n);
t_shift = 0;
step_dur = 100;
c = 25;

injects = round(linspace(142,173,10));
max_amp = zeros(size(injects));
for i=1:numel(injects)
    c = injects(i);
    % step input
    [A,B] = make_compartmental_matrices_simple(Ri,Rm,Cm,dend(:,r),l,n,dend(:,p));

    % numerical (ode)
    dvdt_step = @(t,v) A*v + B*make_u_step(t,t_shift,step_dur,c,n)*Iapp;
    [tp,v_step] = ode23(@(tp,v_step) dvdt_step(tp,v_step), t_span, zeros(n,1), options);
    v_step = v_step';

    % ss
    u_step = zeros(n,1);
    u_step(c) = Iapp;
    v_step_ss = -inv(A)*B*u_step;

    % plot for various compartments (step)
    % change variables to plot
    T = tp./(Rm*Cm);

    % numerical (ode)
    figure('Name','Models of the Neuron Project: Topic 3')
    hold on
    plot(T,v_step(1,:),'m', 'DisplayName', 'Soma')
    plot(T,v_step(c,:),'g', 'DisplayName', 'cmprt '+string(c))
    % plot(T,v_step(c-2,:),'k')
    % plot(T,v_step(c+2,:),'b')
    title('Voltage Given Unit Step Input')
    legend('show')
    xlabel('Time'); ylabel('Voltage (mV)');
    
    max_amp(i)= max(v_step(1,:));
    saveas(gcf, strcat('output_step_',num2str(c),'.png'))
end

figure(); plot(injects, max_amp, 'o-');
xlabel('Time'); ylabel('Peak voltage (mV)');
saveas(gcf, strcat('amp_step_',num2str(c),'.png'))

% ss
% Note from Simon: I suggest plotting steady state the way I did it for the
% shunting inhibition problem. The plot that was here didn't make sense.
% figure('Name','Models of the Neuron Project: Topic 3')
% hold on
% plot(T,v_step_ss(1),'m')
% plot(T,v_step_ss(c),'g')
% plot(T,v_step_ss(c-2),'k')
% plot(T,v_step_ss(c+2),'b')
% title('SS Voltage Given Unit Step Input')
% legend('Soma','Compartment w/Iapp',...
%     'Grandparent compartment','Grandchild compartment')
% xlabel('Time'); ylabel('Voltage (mV)');


%% not sure what this trash is but delete later if everything else works
% ge = zeros(n,1);

% ge(27) = 0.001;

% % response
% u = zeros(n,t_n);
% for i = 1:t_n
%     u(:,i) = make_u(t_span(i),t_shift,c,n)*Iapp;
% end

%% time varying synaptic conductance change
injects = round(linspace(142,173,10));
max_amp_s = zeros(size(injects));
for i=1:numel(injects)
    c = injects(i);
    [A,B,L,C] = make_compartmental_matrices(Ri,Rm,Cm,dend(:,r),l,n,dend(:,p));


    % numerical (ode)
    dvdt_exp = @(t,v) A*v + B*make_u_exp(t,t_shift,c,n)...
        + make_G(make_u_exp(t,t_shift,c,n),C)*v;
    [tp,v_exp_s] = ode23(@(tp, v_exp_s) dvdt_exp(tp, v_exp_s), t_span, zeros(n,1), options);
    v_exp_s = v_exp_s';
    % adjust v_exp (it's actually v*)
    Er = ones(n,1)*(-70);           % mV
    Ee = ones(n,1)*(60);            % mV
    Ei = Er;                        % mV
    v_exp = zeros(size(v_exp_s));
    for i = 1:n
        v_exp(i,:) = v_exp_s(i,:)*(Ee(i)-Ei(i))+Er(i);
    end

    % ss
    % jx: no steady state for transient inputs
    % % intialize
    % u_exp = zeros(2*n,1);
    % g = zeros(n,1);
    % % put in values
    % g(c) = 0.001;
    % u_exp(1:2:end) = g;
    % G = diag(-g./C);
    % % ss voltage calculation
    % v_exp_ss = -inv(A+G)*B*u_exp;

    % plot for various compartments (exp)
    % change variables
    T = tp./(Rm*Cm);

    % numerical (ode)
    figure('Name','Models of the Neuron Project: Topic 3')
    hold on
    plot(T,v_exp(1,:),'m', 'DisplayName', 'Soma')
    plot(T,v_exp(c,:),'g', 'DisplayName', 'cmprt '+string(c))
    % plot(t_span,v_exp(c-2,:),'k')
    % plot(t_span,v_exp(c+2,:),'b')
    title('Voltage Given Decaying Exponential w/Time Varying Conductance Change')
    legend('show')
    xlabel('Time (us)'); ylabel('Voltage (mV)');
    
    max_amp_s(i)= max(v_step(1,:));
    saveas(gcf, strcat('output_exp_',num2str(c),'.png'))
end

figure(); plot(injects, max_amp, 'o-');
xlabel('Time'); ylabel('Peak voltage (mV)');
saveas(gcf, strcat('amp_exp_',num2str(c),'.png'))

% ss
% figure('Name','Models of the Neuron Project: Topic 3')
% hold on
% plot(t_span,v_exp_ss(1),'m')
% plot(t_span,v_exp_ss(c),'g')
% plot(t_span,v_exp_ss(c-2),'k')
% plot(t_span,v_exp_ss(c+2),'b')
% title('SS Voltage Given Decaying Exponential w/Time Varying Conductance Change')
% legend('Soma','Compartment w/Iapp',...
%     'Grandparent compartment','Grandchild compartment')
% xlabel('Time (us)'); ylabel('Voltage (mV)');
