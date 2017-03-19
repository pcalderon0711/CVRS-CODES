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
Q(10,10) = 1;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;

x_nom = zeros(14,1);
x_nom(10) = 40;

f = @(x,u)(myModelWithControl(0,x,u,params));
jac = @(x)(myDerivatives(x,params));

%J = @(z)myCost(z, f, jac, B, R, Q, h, x1, x_nom);
J = @(z)myCostWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta);
%J = @(z)myCostVectorWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta);

z0 = [repmat(x1',N,1); zeros(N,14)];
initz = 'rest';
z0 = z0(:);
% xGuess = (x1+xN)/2;
% z0 = [repmat(xGuess',N,1); zeros(N,14)];
% initz = 'ave';
z0_m = [repmat(xN',N,1);zeros(N,14)];
z0_m = z0_m(:);
z0_a = (z0 + z0_m)/2;
initz = 'mix_exer';

lb = [repmat([60,3,9,2,0,0,0,0,50,35,80,0,0,0],N,1);repmat(-10,N,14)];
ub = [repmat([140,8,18,15,100,100,10,10,200,45,110,1,1,30],N,1);repmat(10,N,14)];

init=cell(10,1);
for i=1:7
    points = (ub-lb).*rand(size(ub))+lb;
    points = points(:);
    init{i} = CustomStartPointSet(points');
    if i==1
        a = points;
    end
end

start=(ub-lb).*rand(size(ub))+lb;
start=start(:);

init{8}=CustomStartPointSet(z0');
init{9}=CustomStartPointSet(z0_a');
init{10}=CustomStartPointSet(z0_m');

%[z,fval,exitflag,output] = lsqnonlin(J,z0,[],[],options);
opts = optimoptions(@fminunc,'MaxIter',1e4,'MaxFunEvals',1e6,'TolFun',1e-9,'PlotFcns',@optimplotfval);
% opts2 = optimset('MaxIter',1e4,'MaxFunEvals',1e6,'PlotFcns',@optimplotfval,'Display','iter-detailed');

% lb = [zeros(N,14); repmat(-1,N,14)];
% ub = [repmat(200,N,14); repmat(1,N,14)];

% problem = createOptimProblem('fmincon','objective',J,'x0',z0,'lb',lb,'ub',ub,'options',opts);
problem = createOptimProblem('fminunc','objective',J,'x0',start,'options',opts);
ms = MultiStart('Display','iter','UseParallel',true);
% gs = GlobalSearch('Display','iter');

%[z,fval,exitflag,output] = lsqnonlin(J,z0,[],[],options);
% [z,fval,exitflag,output] = run(J,z0,opts);

[x,f,exitflag,output] = fminunc(problem)
% [x,f,exitflag,output] = fminsearch(J,z0,opts2)


% [x,f,exitflag,output,solutions] = run(ms,problem,init)