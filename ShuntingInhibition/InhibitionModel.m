%%%%

neuron = load('/Users/Simon/Google Drive/Graduate School/Hopkins/Models of the Neuron/Homeworks/Project1/motn_biophysmodel/ShuntingInhibition/HP72N6B.CNG.txt');
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

Ee = ones(num_compartments,1)*Ee;
Ei = ones(num_compartments,1)*Ei;
Er = ones(num_compartments,1)*Er;

n = ones(num_compartments,1);
[A,B,R,~,C] = make_compartmental_matrices_inhibition(Ri,Rm,Cm,Er,Ee,Ei,neuron(:,6),l,n,neuron(:,7));

%Set up input vector u
u = zeros(2*num_compartments,1);
ge = zeros(num_compartments,1);
gi = zeros(num_compartments,1);

ge(7:10) = 1;
gi(35) = 0.1;

u(1:2:end) = ge;
u(2:2:end) = gi;

%Set up G Matrix
G = zeros(num_compartments);
for j = 1:num_compartments
   G(j,j) = -(ge(j)+gi(j))/C(j); 
end

DvpDt = @(t,vp) A*vp + B*u + G*vp;
tp = 0:0.00001:0.1;

v0 = zeros(num_compartments,1);
[tr,Vp] = ode23(DvpDt,tp,v0);

V = zeros(size(Vp));
x = zeros(num_compartments,1);
for n = 1:num_compartments
    p = neuron(n,7);
    if p > 0
        x(n) = x(p) + l(n);
    end
    V(:,n) = Vp(:,n)*(Ee(n)-Ei(n))+Er(n);
end
X = x./lambda;






