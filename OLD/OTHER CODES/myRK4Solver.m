function next = myRK4Solver(i,y,lambda,u,mode,h,params,B,Q,y_nom)
t=0; % myModel is used in another function requiring t, ignore here

if mode == 'f' % forward
    y_i = y(:,i);
    k1 = myModel(t,y_i,params) + B*u(:,i);
    k2 = myModel(t,y_i+k1*h/2,params) + B*0.5*(u(:,i)+u(:,i+1));
    k3 = myModel(t,y_i+k2*h/2,params) + B*0.5*(u(:,i)+u(:,i+1));
    k4 = myModel(t,y_i+k3*h,params) + B*u(:,i+1);
    next = y_i + (h/6)*(k1+2*k2+2*k3+k4);
elseif mode == 'b' % backward
    lambda_i = lambda(:,i);
    jacobian = myDerivatives(y(:,i),params);
    k1 = -Q*(y(:,i)-y_nom) - jacobian'*lambda_i;
    k2 = -Q*(0.5*(y(:,i)+y(:,i-1))-y_nom) - jacobian'*(lambda_i-k1*h/2);
    k3 = -Q*(0.5*(y(:,i)+y(:,i-1))-y_nom) - jacobian'*(lambda_i-k2*h/2);
    k4 = -Q*(y(:,i-1)-y_nom) - jacobian'*(lambda_i-k3*h);
    next = lambda_i + (h/6)*(k1+2*k2+2*k3+k4);
end
end

