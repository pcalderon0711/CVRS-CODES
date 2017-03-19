function [k, xstar, abserr, relerr, xk] = myNewton(f, f_prime, x0, maxiter, tol)

i = 2;
dim = size(x0,1);
% Initialize vector containing the sequence of approximations
xk = zeros(i,dim);
xk(1,:)=x0;

while true
    % Assign new approx value to old approx variable.
    if i ~= 2
        x0 = x1;
    end
    % Compute new approximation.
    x1 = x0 - pinv(f_prime(x0))*f(x0);
    xk(i,1:dim) = x1';
    % Compute errors
    abserr = sum(abs(x1 - x0));
    relerr = sum(abs((x1 - x0))) / sum(abs(x1));
    % Check if current approximation satisfies the stopping conditions.
	if abserr < tol || relerr < tol || i == maxiter
        % Check if maxiter is reached and warn the user.
		if i == maxiter
			fprintf('Maximum number of iterations is reached.')
		end
		break;
	end
	i = i + 1;
end

% i-1 is the last iteration since list indexing starts with 1.
k = i - 1;
xstar = xk(i,1:dim);
xstar = xstar';