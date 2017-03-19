function J = myCostVectorWithPenalizationAndTimeVaryingWorkload(z, f, jac, B, R, Q, h, x1, x_nom,params_rest,params_exer,t1,t2,t3,t4,t5,t6, alpha, beta, W, A)

%% myCostVectorWithPenalizationAndTimeVaryingWorkload
% Calculates the vector J of costs.
% Inputs:
%        z      : matrix with 2N (N+N) rows and 14 columns
%        f      : vector containing RHS of dot(x) = f
%        jac    : jacobian matrix of f
%        B      : matrix which reshapes controls
%        R      : matrix of weights for controls
%        Q      : matrix of weights for states
%        h      : length of time step
%        x1     : vector containing states at initial rest condition
%        x_nom  : vector containing nominal values of states
%        params : vector containing model parameters
%        alpha  : penalization weight for O2
%        beta   : penalization weight for CO2
%        W      : workload function on t
%        A      : Apesk function on t
% Output:
%        J     : vector of costs of size 2*14*N+N+N
%%
dim = size(z,1); % calculate number of rows of z
N = dim/2; % retrieve number of timesteps N
x = [x1,z(1:N,:)']; % col(x1,x2,x3,...,xN+1)
p = [z(N+1:2*N,:)',Q*(z(N,:)'-x_nom(:,N+1))]; % col(p1,p2,...,pN,pN+1)

% %% Params used to calculate penalty %%%%%
% K_CO2 = params(30);
% k_CO2 = params(31);
% K_a1 = params(32);
% K_a2 = params(33);
% M_O2 = params(16);
% M_CO2 = params(17);
% rho_O2 = params(18);
% rho_CO2 = params(19);

%% Initialize penalty and cost %%%%%
J = zeros(14*N+14*N+2*N,1);

%% Compute cost on x
for k=1:N
    t = (k-1)*h;
    J(14*(k-1)+1:14*k) = x(:,k)-x(:,k+1)+h*(f(t+h,x(:,k+1))-B*inv(R)*B'*p(:,k+1));
end
 
%% Compute cost on p
for k=1:N
    t = (k-1)*h;
    g = jac(t,x(:,k))'*p(:,k);
    J(14*N+14*(k-1)+1:14*N+14*k) = p(:,k+1)-p(:,k)+h*(Q*(x(:,k)-x_nom(:,k))+g);
end

%% Compute penalization cost %%%%%%%%%%%%%%%
for k=2:N+1 
    t = (k-1)*h;
    if (t<t1) || (t2<t)&&(t<t3) || (t4<t)&&(t<t5) || (t>t6)
        params=params_rest;
    else
        params=params_exer;
    end
    K_CO2 = params(30);
    k_CO2 = params(31);
    K_a1 = params(32);
    K_a2 = params(33);
    M_O2 = params(16);
    M_CO2 = params(17);
    rho_O2 = params(18);
    rho_CO2 = params(19);
    P_as = x(1,k);
    P_vs = x(2,k);
    P_aCO2 = x(10,k);
    P_aO2 = x(11,k);
    C_vCO2 = x(12,k);
    C_vO2 = x(13,k);
    R_s = A(t) * C_vO2;
    F_s = (P_as-P_vs)/R_s;
    MR_O2 = M_O2 + rho_O2*W(t);
    MR_CO2 = M_CO2 + rho_CO2*W(t);

    J(2*14*N+(k-1))=alpha*(MR_O2 - F_s * (K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2))^2;
    J(2*14*N+N+(k-1))=beta*(MR_CO2 + F_s * (K_CO2*P_aCO2+k_CO2-C_vCO2))^2;
end

% sum(J(2*14*N+1:2*14*N+N)),sum(J(2*14*N+N+1:end))
% sum(J(1:14*N)),sum(J(14*N+1:14*N+14*N)),sum(J(2*14*N+1:end)),sum(J)
% pause
% J = J'*J;
end