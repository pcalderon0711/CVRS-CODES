function J = myCost(z, f, jac, B, R, Q, h, x1, x_nom)

% z is a matrix with 2N rows and 14 columns; 

dim = size(z,1);
N = dim/2;
x = z(1:N,:)'; % col(x2,x3,...,xN+1)
p = z(N+1:2*N,:)'; % col(p1,p2,...,pN)

J = 0;

for k=1:N
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
    F = v-x(:,k)+h*(f(x(:,k),u)+B*u);
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

end