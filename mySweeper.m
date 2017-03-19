function z = mySweeper(alpha,beta,N)

params_rest = myLoader('parameters_rest.txt','p');
params = myLoader('parameters_exer.txt','p');
x1 = myEquilibriumSolver(params_rest,78.5,40);
xN = myEquilibriumSolver(params,107,40);

T = 10;
h = T/N;

w1 = 0.01;
w2 = 0.01;
R = diag([w1,w2]);
Q = zeros(14,14);
Q(10,10) = 1;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;

W = @(t)(37.5*sin(pi*t)+37.5);
apesk_ave = (177.682+270)/2;
A = @(t)(apesk_ave+0.5*(270-177.682)*sin(pi*t)); %A_pesk function
rp_ave = (1.5446+0.3)/2;
R_p = @(t)(rp_ave+0.5*(.3-1.5446)*sin(pi*t));
w = 'sin';

x_nom = zeros(14,1);
x_nom(10) = 40;

f = @(t,x)(myModelWithControlAndTimeVaryingWorkload(t,x,params,W,A,R_p));
jac = @(t,x)(myDerivativesTimeVarying(t,x,params,W,A,R_p));
J = @(z)myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha,beta, W,A);
z0 = [repmat(x1',N,1); zeros(N,14)];

options = optimset('Display','iter','MaxFunEvals',10000*2*14*N,'MaxIter',1000*2*14*N, 'TolFun',1e-6);
z = lsqnonlin(J,z0,[],[],options);
end

