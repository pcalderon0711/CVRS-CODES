function J = myCostChecker(x, p, f, jac, B, R, Q, h, x1, x_nom)

x = x(:,2:81);
p = p(:,1:80);
J = 0;
N=80;

for k=1:N
    t = k*h;
    if k==1
        v = x1;
    else
        v = x(:,k-1); %x(:,k-1) is xk
    end
    if k==N
        w = Q*(x(:,N)-x_nom);
    else
        w = p(:,k+1);
    end
    u = -pinv(R)*B'*w;
    F = v-x(:,k)+h*(f(t,x(:,k),u)+B*u);
    J = J + F'*F;
end
    
for k=1:N
    if k==1
        v = x1;
    else
        v = x(:,k-1); %x(:,k-1) is xk
    end
    if k==N
        w = Q*(x(:,N)-x_nom); % x(:,N) is x_N+1
    else
        w = p(:,k+1);
    end
    F = w - p(:,k) + h*(Q*(v-x_nom)+jac(v)'*p(:,k));
    J = J + F'*F;
end
