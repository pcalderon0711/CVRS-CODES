function u = myOptimalControlLenhart(y0,params)

% OPTIMAL CONTROL SOLVER
% Implements the forward backward sweep in Lenhart
% Input
% 		y0		: initial condition
%		params  : parameter list

% Output
% 		u		: optimal control (u1 is first row, u2 is second row)

y0 = myLoader('equilibrium_rest.txt','i');
params = myLoader('parameters_exer.txt','p');

delta = 0.001;
terminal = 10;
h = 0.01;
N = terminal/h;
t= linspace(0,terminal,N+1);

u = zeros(2, N+1);
y = zeros(14, N+1);
y(:,1) = y0;
lambda = zeros(14,N+1);

w1 = 1;
w2 = 1;
R = diag([w1,w2]);
Q = zeros(14,14);
Q(10,10) = 1;
B = zeros(14,2);
B(10,1) = 1;
B(14,2) = 1;

y_nom = zeros(14,1);
y_nom(10) = 40;

test = -1;

while test < 0
    old_u = u;
    old_y = y;
    old_lambda = lambda;
    
    % solve y forward in time
    for i=1:N
        i
        %y(:,i+1) = myRK4Solver(i, y,lambda,u,'f',h,params,B,Q,y_nom);
        y(:,i+1) = myBWEulerSolver(i, y,lambda,u,'f',h,params,B,Q,y_nom);
    end
    plot(t,y)
    pause
    
    % initial condition for lambda is dependent on last value for y
    lambda(:,N+1) = Q*(y(:,N+1)-y_nom);
    
    % solve lambda backward in time
    for i=1:N
        j = N+2-i;
        j
        %lambda(:,j-1) = myRK4Solver(j,y,lambda,u,'b',h,params,B,Q,y_nom);
        lambda(:,j-1) = myBWEulerSolver(j,y,lambda,u,'b',h,params,B,Q,y_nom);
    end    
    
    % compute optimal control accdg to Habib
    u1 = -pinv(R)*B'*lambda;

    % u is average of prev and current optimal controls
    u = 0.5*(u1 + old_u);
    plot(t,u)
    pause
    
    %%%% use stopping condition in Lenhart %%%%
    temp1=min(delta*sum(abs(u),2) - sum(abs(old_u - u),2));
    temp2=min(delta*sum(abs(y),2) - sum(abs(old_y - y),2));
    temp3=min(delta*sum(abs(lambda),2) - sum(abs(old_lambda - lambda),2));
    test = min(temp1,min(temp2,temp3));
end
end