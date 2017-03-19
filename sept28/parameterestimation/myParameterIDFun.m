function myParameterIDFun(paramsetID,P0)

close all

load('movingaverageHandPas.mat');
H_experiment = MOVINGAVERAGEHEART;
Pas_experiment = MAPMV;
averageH_rest = MOVINGAVERAGEHEART(1);
averageH_exer = Hexerave;

% param_names = {'c_as', 'c_vs', 'c_ap', 'c_vp', 'c_l', 'c_r', 'R_l', 'R_r', ...
%     'kappa', 'alpha_l', 'alpha_r', 'beta_l', 'beta_r', 'gamma_l', 'gamma_r',...
%     'M_O2', 'M_CO2', 'rho_O2', 'rho_CO2', 'V_tot', 'R_p_rest', 'R_p_exer', 'A_pesk_rest', ...
%     'A_pesk_exer', 'P_IO2', 'P_ICO2', 'V_AO2', 'V_ACO2', 'V_TO2', 'V_TCO2', 'K_CO2', 'k_CO2',...
%     'K_a1', 'K_a2', 'q_as', 'q_aco2', 'w1', 'w2', 'alpha', 'beta'};

paramset_idx = myIdentifyParam(paramsetID);

a = [33.11,3.3,5.87];

params_rest = myLoader('parameters_rest.txt','p'); % load rest parameters
params_exer = myLoader('parameters_exer.txt','p'); % load exercise parameters
N = 150; % set number of time points
T = 15; % set terminal time
ts = 5;
h = T/N; % calculate length of each subinterval
t = linspace(0,T,N+1); % partition time interval [0,T] into N equal subintervals

options = optimset('Display','iter','MaxFunEvals',1e+008,'MaxIter',1e+008);
J = @(P)myExperimentalCost2(P,a,Pas_experiment,H_experiment,params_rest,params_exer,averageH_rest,averageH_exer,N,ts,h,paramset_idx);
[P,fval,exitflag,output] = lsqnonlin(J,P0,[],[],options); % use least squares to get states
% [P,fval,exitflag,output] = fminsearch(J,P0,options); % use least squares to get states

save('param_estimation.mat','t','T','N','ts','P0','paramsetID','P','fval','exitflag','output','params_rest','params_exer');
%%

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
W = @(t)myTrans(t,ts,0,75);
A = @(t)myTrans(t,ts,params_rest(23),params_exer(23));
R_p = @(t)myTrans(t,ts,params_rest(22),params_exer(22));
w = 'constant';

f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params_rest,params_exer,ts,W,A,R_p)); % construct vector of the dynamics of states
jac = @(t,x)(myDerivativesTimeVarying(t,x,params_rest,params_exer,ts,A,R_p)); % construct matrix of the jacobian of f
CV = @(z)myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, xz,x_nom, params_rest,params_exer,ts, alpha,beta, W,A); % construct cost vector
z0 = [repmat(x1',N,1); zeros(N,14)]; % construct vector of initial conditions
options = optimset('Display','off','MaxFunEvals',1e+008,'MaxIter',1e+008);

z = lsqnonlin(CV,z0,[],[],options); % use least squares to get states
x = [x1,z(1:N,:)']; % each column is one time step, vector of states
p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)]; % each column is one time step, vector of adjoints
u = -inv(R)*B'*p; % solve for control

if strcmp(w,'constant')
    prefix = 'constant';
else
    prefix = 'timevar';
end

save(sprintf('%s_%s_control(N%f,T%f).mat',prefix,w,N,T),'t','T','x','p','u','fval','exitflag','output','params_rest','params_exer','ts','N','z0','w','W');

myOptimalControlPlotterTimeVarying(t,T,u,x,p,N,alpha,beta,prefix,w,Pas_experiment,H_experiment); % plot the states and controls
myOptimalControlFig(t,u,x,p,N,alpha,beta,prefix,w,Pas_experiment,H_experiment);

end
