function [T,Y] = myBWEulerSolver(f, interval, y0, N)

intlength = interval(end) - interval(1);
h = intlength/N;
T = interval(1):h:interval(end);
Y = zeros(length(T),length(y0));
Y(1,:) = y0;
opts = optimset('Diagnostics','off', 'Display','off');

for i=1:N
    i
    curr = Y(i,:)';
    implicit = @(next)(next - curr - h*f(T(i+1),next));
    Y(i+1,:) = fsolve(implicit, curr, opts);
end    