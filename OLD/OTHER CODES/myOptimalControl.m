%function [x,fval,exitflag,output] = myOptimalControl()

params_rest = myLoader('parameters_rest.txt','p');
params = myLoader('parameters_exer.txt','p');
% x1 = myLoader('equilibrium_rest_old.txt','i');
% xN = myLoader('equilibrium_exer_old.txt','i');
x1 = myEquilibriumSolver(params_rest,78.5,40);
xN = myEquilibriumSolver(params,107,40);

N = 20;
alpha = 0;
beta = 0;

T = 10;
h = T/N;
t = linspace(0,T,N+1);

w1 = 0.01;
w2 = 0.01;
R = diag([w1,w2]);
Q = zeros(14,14);
% Q(1,1) = 1; %penalize Pas
Q(10,10) = 1;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;

x_nom = zeros(14,1);
% x_nom(1) = 122.93; %Pas
x_nom(10) = 40;

f = @(x,u)(myModelWithControl(0,x,u,params));
jac = @(x)(myDerivatives(x,params));

%J = @(z)myCost(z, f, jac, B, R, Q, h, x1, x_nom);
% J = @(z)myCostWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta);
J = @(z)myCostVectorWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta);
z0 = [repmat(x1',N,1); zeros(N,14)];
initz = 'rest';
% xGuess = (x1+xN)/2;
% z0 = [repmat(xGuess',N,1); zeros(N,14)];
% initz = 'ave';
% z0 = [repmat(xN',N,1);zeros(N,14)];
% initz = 'mix_exer';

options = optimset('Display','iter','MaxFunEvals',10000*2*14*N,'MaxIter',1000*2*14*N, 'TolFun',10e-9);
[z,fval,exitflag,output] = lsqnonlin(J,z0,[],[],options);
% [z,fval,exitflag,output] = fminunc(J,z0,options);
x = [x1,z(1:N,:)']; % each column is one time step
p = [z(N+1:2*N,:)',Q*(x(:,N+1)-x_nom)]; % each column is one time step
u = -pinv(R)*B'*p;

cd('mats');
save(sprintf('control(N%f,a%f,b%f,%s).mat',N,alpha,beta,initz),'t','x','p','u','fval','exitflag','output','params','N','alpha','beta','z0');
cd ..;

myOptimalControlPlotter(t,u,x,p,N,alpha,beta,initz);
%J = myCostChecker(z, f, jac, B, R, Q, h, x1, x_nom);