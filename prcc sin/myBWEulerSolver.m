function Y = myBWEulerSolver(f, h,tspan, y0)

Y = zeros(length(tspan),length(y0));
Y(1,:) = y0;
opts = optimset('Diagnostics','off', 'Display','off');

for i=1:length(tspan)-1
    curr = Y(i,:)';
    implicit = @(next)(next - curr - h*f(tspan(i+1),next)');
    Y(i+1,:) = fsolve(implicit, curr, opts);
end    