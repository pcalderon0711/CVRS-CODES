close all

params_rest = myLoader('parameters_rest.txt','p'); % load rest parameters
params_exer = myLoader('parameters_exer.txt','p'); % load exercise parameters
x1 = myEquilibriumSolver(params_rest,78.5,40); % calculate rest equilibrium
xz = myEquilibriumSolver(params_exer,107,40); % calculate exercise equilibrium

N = 150; % set number of time points
alpha = 20; % set penalization parameter on O2
beta = 5; % set penalization parameter on CO2

T = 15; % set terminal time
ts = 5;
h = T/N; % calculate length of each subinterval
t = linspace(0,T,N+1); % partition time interval [0,T] into N equal subintervals

w1 = 0.01; % set weight on the magnitude of control 1 (deriv of heartrate) 
w2 = 0.01; % set weight on the magintude of control 2 (deriv of ventilation)
R = diag([w1,w2]); % construct matrix of weights on magnitude of controls
Q = zeros(14,14); % construct matrix of weights on controlled states
Q(10,10) = 1; % monitor difference of PaCO2 from nominal value
B = zeros(14,2);
B(9,1) = 1; % the 1st control is the derivative of the 9th variable, heart rate
B(14,2) = 1; % the 2nd control is the derivative of the 14th variable, ventilation

x_nom = zeros(14,1); % construct vector of nominal values
x_nom(10) = 40; % set nominal value of PaCO2

%% case of constant workload %%%%%%%
% Q(1,1) = 1; % monitor difference of Pas from nominal value
% x_nom(1) = xz(1);  % set nominal value of Pas 

%% set workload, apesk and Rp

% W = @(t)myTrans(t,ts,0,75);
% A = @(t)myTrans(t,ts,177.682,270);
% R_p = @(t)myTrans(t,ts,1.5446,0.3);
% w = 'constant';

% W = @(t)(37.5*sin(pi*t)+37.5);
% apesk_ave = (177.682+270)/2;
% A = @(t)(apesk_ave+0.5*(270-177.682)*sin(pi*t)); %A_pesk function
% rp_ave = (1.5446+0.3)/2;
% R_p = @(t)(rp_ave+0.5*(.3-1.5446)*sin(pi*t));
% w = 'sin';

% W = @(t)myTrans(t,ts,0,37.5*sin(pi*(t-ts-0.5))+37.5);
% apesk_ave = (177.682+270)/2;
% A = @(t)myTrans(t,ts,177.682,(apesk_ave+0.5*(270-177.682)*sin(pi*(t-ts-0.5))));
% rp_ave = (1.5446+0.3)/2;
% R_p = @(t)myTrans(t,ts,1.5446,(rp_ave+0.5*(.3-1.5446)*sin(pi*(t-ts-0.5))));
% w = 'sinEC';
% W = @(t)(37.5*sin(pi*(t-0.5))+37.5);
% apesk_ave = (177.682+270)/2;
% A = @(t)(apesk_ave+0.5*(270-177.682)*sin(pi*(t-0.5))); %A_pesk function
% rp_ave = (1.5446+0.3)/2;
% R_p = @(t)(rp_ave+0.5*(.3-1.5446)*sin(pi*(t-0.5)));
% w = 'sinEC';
%  


W = @(t)myTrans(t,ts,0,mySquarewave(t-ts,0,75));
A = @(t)myTrans(t,ts,177.682,mySquarewave(t,177.682,270)); %A_pesk function
R_p = @(t)myTrans(t,ts,1.5446,mySquarewave(t,1.5446,0.3));
w = 'square';

% W = @(t)myLongSquare(t,0,75);
% A = @(t)myLongSquare(t,177.682,270); %A_pesk function
% R_p = @(t)myLongSquare(t,1.5446,0.3);
% w = 'longsquare';
% 
% W = @(t)myPlateau(t,0,75);
% A = @(t)myPlateau(t,177.682,270); %A_pesk function
% R_p = @(t)myPlateau(t,1.5446,0.3);
% w = 'plateau';

% W = @(t)mySawtooth(t,0,75);
% A = @(t)mySawtooth(t,177.682,270); %A_pesk function
% R_p = @(t)mySawtooth(t,1.5446,0.3);
% w = 'sawtooth';

% W = @(t)myRamp(t,0,75);
% A = @(t)myRamp(t,177.682,270); %A_pesk function
% R_p = @(t)myRamp(t,1.5446,0.3);
% w = 'ramp';

% W = @(t)myReverseRamp(t,0,75);
% A = @(t)myReverseRamp(t,177.682,270); %A_pesk function
% R_p = @(t)myReverseRamp(t,1.5446,0.3);
% w = 'reverseramp';

% a = @(w)myLinInterpolater(w,177.682,270);
% r = @(w)myLinInterpolater(w,1.5446,0.3);
% W = @(t)myUF(t,75,20,40,30,50);
% A = @(t)myUF(t,a(75),a(20),a(40),a(30),a(50)); %A_pesk function
% R_p = @(t)myUF(t,r(75),r(20),r(40),r(30),r(50));
% w = 'ultrafiltration';

%% optimization stuff
f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params_rest,params_exer,ts,W,A,R_p)); % construct vector of the dynamics of states
jac = @(t,x)(myDerivativesTimeVarying(t,x,params_rest,params_exer,ts,A,R_p)); % construct matrix of the jacobian of f
J = @(z)myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, xz, x_nom,params_rest,params_exer,ts, alpha,beta, W,A); % construct cost vector
z0 = [repmat(x1',N,1); zeros(N,14)]; % construct vector of initial conditions

% options = optimset('Display','iter','MaxFunEvals',1e+008,'MaxIter',1e+008);
options = optimset('Display','iter','MaxFunEvals',1e+008,'MaxIter',1e+008);

% [z,fval,exitflag,output] = fminunc(J,z0,options);
[z,fval,exitflag,output] = lsqnonlin(J,z0,[],[],options); % use least squares to get states
x = [x1,z(1:N,:)']; % each column is one time step, vector of states
p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)]; % each column is one time step, vector of adjoints
u = -inv(R)*B'*p; % solve for control

if strcmp(w,'constant')
    prefix = 'constant';
else
    prefix = 'timevar';
end

cd('mats');
save(sprintf('%s_%s_control(N%f,a%f,b%f,T%f).mat',prefix,w,N,alpha,beta,T),'t','T','x','p','u','fval','exitflag','output','params_rest','params_exer','ts','N','alpha','beta','z0','w','W');
cd ..;

myOptimalControlPlotterTimeVarying(t,T,u,x,p,N,alpha,beta,prefix,w); % plot the states and controls
myOptimalControlFig(t,u,x,p,N,alpha,beta,prefix,w)
