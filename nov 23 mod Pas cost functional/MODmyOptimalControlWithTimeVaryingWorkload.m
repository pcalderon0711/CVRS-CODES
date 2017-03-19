close all

params_rest = myLoader('parameters_rest.txt','p');
params = myLoader('parameters_exer.txt','p');
x1 = myEquilibriumSolver(params_rest,78.5,40);
xz = myEquilibriumSolver(params,107,40);

N = 80;
alpha = 0;
beta = 0;

T = 10;
h = T/N;
t = linspace(0,T,N+1);

w1 = 0.01;
w2 = 0.01;
R = diag([w1,w2]);
Q = zeros(14,14);
Q(1,1) = 1;
Q(10,10) = 1;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;

%x_nom = zeros(14,1);
%x_nom(1) = xz(1);
%x_nom(10) = 40;
x_nom = zeros(14,N+1);
x_nom(10,1:N+1) = 40*ones(1,N+1);
 
x_nom(1,1:N+1) = xz(1)*ones(1,N+1);
for i=1:N+1
    j = i-1;
    if 50<=j<70 || ((90<=j)&&(j<=110)) || ((130<=j)&&(j<=150))
        x_nom(1,i) = xz(1);
    end
end
      
% W = @(t)(75);
% A = @(t)(270);
% w = 'constant';

% W = @(t)(37.5*sin(pi*t)+37.5);
% apesk_ave = (177.682+270)/2;
% A = @(t)(apesk_ave+0.5*(270-177.682)*sin(pi*t)); %A_pesk function
% rp_ave = (1.5446+0.3)/2;
% R_p = @(t)(rp_ave+0.5*(1.5446-.3)*sin(pi*t));
% w = 'sin';
 
% W = @(t)mySquarewave(t,0,75);
% A = @(t)mySquarewave(t,177.682,270); %A_pesk function
% R_p = @(t)mySquarewave(t,1.5446,0.3);
% w = 'square';

% W = @(t)mySawtooth(t,0,75);
% A = @(t)mySawtooth(t,177.682,270); %A_pesk function
% R_p = @(t)mySawtooth(t,1.5446,0.3);
% w = 'sawtooth';

% W = @(t)myRamp(t,0,75);
% A = @(t)myRamp(t,177.682,270); %A_pesk function
% R_p = @(t)myRamp(t,1.5446,0.3);
% w = 'ramp';

% a = @(w)myLinInterpolater(w,177.682,270);
% r = @(w)myLinInterpolater(w,1.5446,0.3);
% W = @(t)myUF(t,75,20,40,30,50);
% A = @(t)myUF(t,a(75),a(20),a(40),a(30),a(50)); %A_pesk function
% R_p = @(t)myUF(t,r(75),r(20),r(40),r(30),r(50));
% w = 'ultrafiltration';

W = @(t)myLongSquare(t,0,75);
A = @(t)myLongSquare(t,177.682,270); %A_pesk function
R_p = @(t)myLongSquare(t,1.5446,0.3);
% A = @(t)(270); %A_pesk function
% R_p = @(t)(0.3);
w = 'longsquare';

f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params,W,A,R_p));
jac = @(t,x)(myDerivativesTimeVarying(t,x,params,A,R_p));
J = @(z)MODmyCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha,beta, W,A);
z0 = [repmat(x1',N,1); zeros(N,14)];

options = optimset('Display','iter','MaxFunEvals',1e+008,'MaxIter',1e+008);
% options = optimset('Display','iter','MaxFunEvals',1e+006,'MaxIter',1e+006);

[z,fval,exitflag,output] = lsqnonlin(J,z0,[],[],options);
x = [x1,z(1:N,:)']; % each column is one time step
% p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom(:,N+1))]; % each column is one time step
p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)]; % each column is one time step
u = -inv(R)*B'*p;

if strcmp(w,'constant')
    prefix = 'constant';
else
    prefix = 'NOPENPASIN2MINtimevar';
end

cd('mats');
save(sprintf('%s_%s_control(N%f,a%f,b%f,T%f).mat',prefix,w,N,alpha,beta,T),'t','T','x','p','u','fval','exitflag','output','params','N','alpha','beta','z0','w','W');
cd ..;

myOptimalControlPlotterTimeVarying(t,T,u,x,p,N,alpha,beta,prefix,w);
