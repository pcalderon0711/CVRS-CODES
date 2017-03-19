function next = myBWEulerSolver(i,y,lambda,u,mode,h,params,B,Q,y_nom)
t=0; % myModel is used in another function requiring t, ignore here

if mode == 'f'
    y0 = myLoader('equilibrium_rest.txt','i');
    y_i = y(:,i);
    % implicit is the implicit equation we need to solve
    % next is the value at the i+1'th step, i.e. the value we need to solve
    implicit = @(next)(next-y_i-h*(myModel(t,next,params)+B*u(:,i+1)));
    % jacobian of the function f
    jacobian_f = @(x0)(myDerivatives(x0,params));
    % jacobian of implicit
    jacobian_implicit = @(x0)(eye(size(jacobian_f(x0)))-jacobian_f(x0));
    % next value cannot be obtained in closed form, use Newton's method
    % we take y0 as an initial guess for the next step
    [~, next, ~, ~, ~] = myNewton(implicit, jacobian_implicit,y0, 1000, 1e-3);
    fclose('all');

elseif mode == 'b'
    lambda_i = lambda(:,i);
    jacobian = myDerivatives(y(:,i-1),params); % at the i-1 step

    % next is the value at the i-1th step
    % next value can be obtained as a closed form expression
    next = pinv(eye(size(jacobian)) - h*jacobian')*(lambda_i + h*Q*(y(:,i-1)-y_nom));
%     if i == 2
%         next;
%         pause;
%     end

end

end


% % Sanity check
% N=30;
% h=0.2;
% y=zeros(1,31);
% y(1)=1;
% for i = 1:N
%     y(i+1)=y(i)/(1+h);
% 
% 
% end