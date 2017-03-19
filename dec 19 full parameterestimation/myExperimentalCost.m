function J = myExperimentalCost(P,a,pas_data,h_data,params_rest,params_exer,averageH_rest,averageH_exer,N,ts,h,alpha,beta)

params_rest(1) = P(3);
params_rest(30) = P(4);
params_rest(31) = P(5);
params_exer(1) = P(3);
params_exer(30) = P(4);
params_exer(31) = P(5);

x1 = myEquilibriumSolver(params_rest,averageH_rest,40); % calculate rest equilibrium
xz = myEquilibriumSolver(params_exer,averageH_exer,40); % calculate exercise equilibrium

w1 = P(8); % set weight on the magnitude of control 1 (deriv of heartrate) 
w2 = P(9); % set weight on the magintude of control 2 (deriv of ventilation)
R = diag([w1,w2]); % construct matrix of weights on magnitude of controls
Q = zeros(14,14); % construct matrix of weights on controlled states
Q(10,10) = P(7); % monitor difference of PaCO2 from nominal value
B = zeros(14,2);
B(9,1) = 1; % the 1st control is the derivative of the 9th variable, heart rate
B(14,2) = 1; % the 2nd control is the derivative of the 14th variable, ventilation

x_nom = zeros(14,1); % construct vector of nominal values
x_nom(10) = 40; % set nominal value of PaCO2

%% case of constant workload %%%%%%%
Q(1,1) = P(6); % monitor difference of Pas from nominal value
x_nom(1) = xz(1);  % set nominal value of Pas 

%% set workload, apesk and Rp

W = @(t)myTrans(t,ts,0,75);
A = @(t)myTrans(t,ts,P(1),P(2));
R_p = @(t)myTrans(t,ts,1.5446,0.3);
w = 'constant';

a1rest = a(1);
a1exer = a(2);
a2exer = a(3);

%% optimization stuff
f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params_rest,params_exer,ts,W,A,R_p)); % construct vector of the dynamics of states
jac = @(t,x)(myDerivativesTimeVarying(t,x,params_rest,params_exer,ts,A,R_p)); % construct matrix of the jacobian of f
CV = @(z)myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, xz,x_nom, params_rest,params_exer,ts, alpha,beta, W,A); % construct cost vector
z0 = [repmat(x1',N,1); zeros(N,14)]; % construct vector of initial conditions
options = optimset('Display','off','MaxFunEvals',1e+008,'MaxIter',1e+008);

z = lsqnonlin(CV,z0,[],[],options); % use least squares to get states
x = [x1,z(1:N,:)']; % each column is one time step, vector of states
% p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)]; % each column is one time step, vector of adjoints
% u = -inv(R)*B'*p; % solve for control

K=zeros(1,N+1);
L=zeros(1,N+1);

for i=2:(N/3)+1 % in total, there are N+1 data points, but first data point is at t=0, so we start at n=2
    K(i) = a1rest*(x(1,i)-pas_data(i))^2;
end
for i=(N/3)+2:N+1
    K(i) = a1exer*(x(1,i)-pas_data(i))^2;
    L(i) = a2exer*(x(9,i)-h_data(i))^2;
end

J = [K L((N/3)+2:N+1)];
end

