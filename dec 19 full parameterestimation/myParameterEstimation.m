close all

load('movingaverageHandPas.mat');
H_experiment = MOVINGAVERAGEHEART;
Pas_experiment = MAPMV;
averageH_rest = MOVINGAVERAGEHEART(1);
averageH_exer = Hexerave;

param_names = {'c_as', 'c_vs', 'c_ap', 'c_vp', 'c_l', 'c_r', 'R_l', 'R_r', ...
    'kappa', 'alpha_l', 'alpha_r', 'beta_l', 'beta_r', 'gamma_l', 'gamma_r',...
    'M_O2', 'M_CO2', 'rho_O2', 'rho_CO2', 'V_tot', 'R_p_rest', 'A_pesk_rest', ...
    'A_pesk_exer', 'P_IO2', 'P_ICO2', 'V_AO2', 'V_ACO2', 'V_TO2', 'V_TCO2', 'K_CO2', 'k_CO2',...
    'K_a1', 'K_a2', 'q_as', 'q_aco2', 'w1', 'w2'};

a = [33.11,3.3,5.87];

params_rest = myLoader('parameters_rest.txt','p'); % load rest parameters
params_exer = myLoader('parameters_exer.txt','p'); % load exercise parameters
N = 150; % set number of time points
T = 15; % set terminal time
ts = 5;
h = T/N; % calculate length of each subinterval
t = linspace(0,T,N+1); % partition time interval [0,T] into N equal subintervals
alpha = 20;
beta = 5;

c_as0 = 0.01016;
K_CO20 = 0.0057;
k_CO20 = 0.224;

%Apeskrest,Apeskexer,cas,Kco2,kco2,qas,qaco2,w1,w2
P0 = [177.682,270,c_as0,K_CO20,k_CO20,1,1,0.01,0.01];
options = optimset('Display','iter','MaxFunEvals',1e+008,'MaxIter',1e+008);
J = @(P)myExperimentalCost(P,a,Pas_experiment,H_experiment,params_rest,params_exer,averageH_rest,averageH_exer,N,ts,h,alpha,beta);
[P,fval,exitflag,output] = lsqnonlin(J,P0,[],[],options); % use least squares to get states
% [P,fval,exitflag,output] = fminsearch(J,P0,options); % use least squares to get states

save('param_estimation.mat','t','T','N','alpha','beta','ts','P0','P','fval','exitflag','output');
%%

params_rest = myLoader('parameters_rest.txt','p'); % load rest parameters
params_exer = myLoader('parameters_exer.txt','p'); % load exercise parameters
params_rest(1) = P(3);
params_rest(30) = P(4);
params_rest(31) = P(5);
params_exer(1) = P(3);
params_exer(30) = P(4);
params_exer(31) = P(5);
x1 = myEquilibriumSolver(params_rest,averageH_rest,40); % calculate rest equilibrium
xz = myEquilibriumSolver(params_exer,averageH_exer,40); % calculate exercise equilibrium
w1 = P(8);
w2 = P(9);
R = diag([w1,w2]); % construct matrix of weights on magnitude of controls
Q = zeros(14,14); % construct matrix of weights on controlled states
Q(10,10) = P(7); % monitor difference of PaCO2 from nominal value
B = zeros(14,2);
B(9,1) = 1; % the 1st control is the derivative of the 9th variable, heart rate
B(14,2) = 1; % the 2nd control is the derivative of the 14th variable, ventilation
x_nom = zeros(14,1); % construct vector of nominal values
x_nom(10) = 40; % set nominal value of PaCO2
Q(1,1) = P(6); % monitor difference of Pas from nominal value
x_nom(1) = xz(1);  % set nominal value of Pas 
W = @(t)myTrans(t,ts,0,75);
A = @(t)myTrans(t,ts,P(1),P(2));
R_p = @(t)myTrans(t,ts,1.5446,0.3);
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

save(sprintf('%s_%s_control(N%f,a%f,b%f,T%f).mat',prefix,w,N,alpha,beta,T),'t','T','x','p','u','fval','exitflag','output','params_rest','params_exer','ts','N','alpha','beta','z0','w','W');

myOptimalControlPlotterTimeVarying(t,T,u,x,p,N,alpha,beta,prefix,w,Pas_experiment,H_experiment); % plot the states and controls
myOptimalControlFig(t,u,x,p,N,alpha,beta,prefix,w,Pas_experiment,H_experiment);
