params_rest = myLoader('parameters_rest.txt','p');
params = myLoader('parameters_exer.txt','p');
x1 = myEquilibriumSolver(params_rest,78.5,40);
%x1 = x1*1.1;

N = 80;
T = 500;
h = T/N;

% W = @(t)(75);
% A = @(t)(270);
% w = 'constant';

W = @(t)(37.5*sin(pi*t)+37.5);
A = @(t)(177.682+(270-177.682)*sin(pi*t));

f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params,W,A));
[T,Y] = ode15s(f,[0 T], x1,[]);
