terminal = 6;
h = 0.01;
N = terminal/h;
f = @(y)(-y);
g = @(t,y)(-y);
y = zeros(1,N+1);
z = zeros(1,N+1);
y(1)=1;
z(1)=1;
t= linspace(0,terminal,N+1);

for i=1:N
    y_i = y(:,i);
    implicit = @(next)(next-(y_i/(1+h)));
    jacobian_implicit = @(next)1;
    [~, y(:,i+1), ~, ~, ~] = myNewton(implicit, jacobian_implicit,0, 1000, 1e-3);
    
    k1 = f(y_i);
    k2 = f(y_i+k1*h/2);
    k3 = f(y_i+k2*h/2);
    k4 = f(y_i+k3*h);
    z(:,i+1) = y_i + (h/6)*(k1+2*k2+2*k3+k4);
    
end
[T,Y] = ode15s(g,[0 terminal], 1);
plot(t,y,t,z,T,Y)
