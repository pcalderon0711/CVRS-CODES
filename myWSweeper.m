function z = myWSweeper(w,alpha,beta,N)

params_rest = myLoader('parameters_rest.txt','p');
params = myLoader('parameters_exer.txt','p');
x1 = myEquilibriumSolver(params_rest,78.5,40);

T = 10;
h = T/N;

w1 = 1;
w2 = 1;
R = diag([w1,w2]);
Q = zeros(14,14);
Q(10,10) = w;
B = zeros(14,2);
B(9,1) = 1;
B(14,2) = 1;

x_nom = zeros(14,1);
x_nom(10) = 40;

f = @(x,u)(myModelWithControl(0,x,u,params));
jac = @(x)(myDerivatives(x,params));

J = @(z)myCostVectorWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta);
z0 = [repmat(x1',N,1); zeros(N,14)];

options = optimset('Display','iter','MaxFunEvals',10000*2*14*N,'MaxIter',1000*2*14*N, 'TolFun',10e-9);
z = lsqnonlin(J,z0,[],[],options);
end

