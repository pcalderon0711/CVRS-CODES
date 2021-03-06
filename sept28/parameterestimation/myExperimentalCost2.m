function J = myExperimentalCost2(P,a,pas_data,h_data,params_rest,params_exer,averageH_rest,averageH_exer,N,ts,h,paramset_idx)

alpha = 20;
beta = 5;
w1 = 0.01; % set weight on the magnitude of control 1 (deriv of heartrate) 
w2 = 0.01; % set weight on the magintude of control 2 (deriv of ventilation)
Q = zeros(14,14); % construct matrix of weights on controlled states
Q(1,1) = 1; % monitor difference of Pas from nominal value
Q(10,10) = 1; % monitor difference of PaCO2 from nominal value
B = zeros(14,2);
B(9,1) = 1; % the 1st control is the derivative of the 9th variable, heart rate
B(14,2) = 1; % the 2nd control is the derivative of the 14th variable, ventilation

for i=1:length(paramset_idx)
    curr = paramset_idx(i);
    if curr == 1
        params_rest(1) = P(i);
        params_exer(1) = P(i);
    elseif curr == 2
        params_rest(2) = P(i);
        params_exer(2) = P(i);
    elseif curr == 3
        params_rest(3) = P(i);
        params_exer(3) = P(i);
    elseif curr == 4
        params_rest(4) = P(i);
        params_exer(4) = P(i);
    elseif curr == 5
        params_rest(5) = P(i);
        params_exer(5) = P(i);
    elseif curr == 6
        params_rest(6) = P(i);
        params_exer(6) = P(i);
    elseif curr == 7
        params_rest(7) = P(i);
        params_exer(7) = P(i);
    elseif curr == 8
        params_rest(8) = P(i);
        params_exer(8) = P(i);
    elseif curr == 9
        params_rest(9) = P(i);
        params_exer(9) = P(i);
     elseif curr == 10
        params_rest(10) = P(i);
        params_exer(10) = P(i);
    elseif curr == 11
        params_rest(11) = P(i);
        params_exer(11) = P(i);
    elseif curr == 12
        params_rest(12) = P(i);
        params_exer(12) = P(i);
    elseif curr == 13
        params_rest(13) = P(i);
        params_exer(13) = P(i);
    elseif curr == 14
        params_rest(14) = P(i);
        params_exer(14) = P(i);
    elseif curr == 15
        params_rest(15) = P(i);
        params_exer(15) = P(i);
    elseif curr == 16
        params_rest(16) = P(i);
        params_exer(16) = P(i);
    elseif curr == 17
        params_rest(17) = P(i);
        params_exer(17) = P(i);
    elseif curr == 18
        params_rest(18) = P(i);
        params_exer(18) = P(i);
    elseif curr == 19
        params_rest(19) = P(i);
        params_exer(19) = P(i);
    elseif curr == 20 %Vtot
        params_rest(21) = P(i);
        params_exer(21) = P(i);
     elseif curr == 21 %Rprest
        params_rest(22) = P(i);
    elseif curr == 22 %Rpexer
        params_exer(22) = P(i);
    elseif curr == 23 %Apeskrest
        params_rest(23) = P(i);
    elseif curr == 24 %Apeskexer
        params_exer(23) = P(i);
    elseif curr == 25
        params_rest(24) = P(i);
        params_exer(24) = P(i);
    elseif curr == 26
        params_rest(25) = P(i);
        params_exer(25) = P(i);
    elseif curr == 27
        params_rest(26) = P(i);
        params_exer(26) = P(i);
    elseif curr == 28
        params_rest(27) = P(i);
        params_exer(27) = P(i);
    elseif curr == 29
        params_rest(28) = P(i);
        params_exer(28) = P(i);
    elseif curr == 30
        params_rest(29) = P(i);
        params_exer(29) = P(i);
    elseif curr == 31
        params_rest(30) = P(i);
        params_exer(30) = P(i);
    elseif curr == 32
        params_rest(31) = P(i);
        params_exer(31) = P(i);
    elseif curr == 33
        params_rest(32) = P(i);
        params_exer(32) = P(i);
    elseif curr == 34
        params_rest(33) = P(i);
        params_exer(33) = P(i);
    elseif curr == 35 %qas
        Q(1,1) = P(i);
    elseif curr == 36 %qaco2
        Q(10,10) = P(i);
    elseif curr == 37
        w1 = P(i);
    elseif curr == 38
        w2 = P(i);
    elseif curr == 39
        alpha = P(i);
    elseif curr == 40
        beta = P(i); 
    end
end

R = diag([w1,w2]); % construct matrix of weights on magnitude of controls

x1 = myEquilibriumSolver(params_rest,averageH_rest,40); % calculate rest equilibrium
xz = myEquilibriumSolver(params_exer,averageH_exer,40); % calculate exercise equilibrium

x_nom = zeros(14,1); % construct vector of nominal values
x_nom(10) = 40; % set nominal value of PaCO2
x_nom(1) = xz(1);  % set nominal value of Pas 

%% set workload, apesk and Rp

W = @(t)myTrans(t,ts,0,75);
A = @(t)myTrans(t,ts,params_rest(23),params_exer(23));
R_p = @(t)myTrans(t,ts,params_rest(22),params_exer(22));
w = 'constant';

a1rest = a(1);
a1exer = a(2);
a2exer = a(3);

%% optimization stuff
f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params_rest,params_exer,ts,W,A,R_p)); % construct vector of the dynamics of states
jac = @(t,x)(myDerivativesTimeVarying(t,x,params_rest,params_exer,ts,A,R_p)); % construct matrix of the jacobian of f
CV = @(z)myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, xz,x_nom, params_rest,params_exer,ts, alpha,beta, W,A); % construct cost vector
z0 = [repmat(x1',N,1); zeros(N,14)]; % construct vector of initial conditions
options = optimset('Display','off','MaxFunEvals',1e+008,'MaxIter',1e+008,'Algorithm','levenberg-marquardt');

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
J = sqrt(J);
end

