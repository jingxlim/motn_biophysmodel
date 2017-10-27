%%%%

neuron = load('/Users/Simon/Google Drive/Graduate School/Hopkins/Models of the Neuron/Homeworks/Project1/motn_biophysmodel/ShuntingInhibition/HP72N6B.CNG.txt');
axon = find(neuron(:,2)==2);
neuron(axon,:) = [];

Ri = 100; %Ohm - cm
Rm = 10000; %Ohm - cm^2
Cm = 1; %microF/cm^2
Iapp = 10e-9; %mA
Er = -50; %Mv
Ee = 50; %Mv
Ei = -70; %Mv

%Convert micrometers to centimeters
neuron(:,3:6) = neuron(:,3:6)./10000;

lambda = sqrt((Rm./Ri).*(neuron(:,6)./2));

l = zeros(length(neuron),1);
for n = 1:length(neuron)
    p = neuron(n,7);
    if p < 1
        dx = 0-neuron(n,3);
        dy = 0-neuron(n,4);
        dz = 0-neuron(n,5);
    else
        dx = neuron(p,3)-neuron(n,3);
        dy = neuron(p,4)-neuron(n,4);
        dz = neuron(p,5)-neuron(n,5);
    end
    l(n) = sqrt(dx^2 + dy^2 + dz^2);
end

L = l./lambda;

disp('Is the compartment size constraint met for all compartments?');
if any(L>0.1)
    disp('No, it is not met for all compartments.');
else
    disp('Yes!');
end

n = ones(length(neuron),1);
[A,B,R,L] = make_compartmental_matrices_inhibition(Ri,Rm,Cm,Er,Ee,Ei,neuron(:,6),l,n,neuron(:,7));

%Set up input vector u



%Set up G Matrix



%dvp/dt = @(vp,t) A*vp(t) + B*u










