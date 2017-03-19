function J = myCostVectorWithPenalization(z, f, jac, B, R, Q, h, x1, x_nom, params, alpha, beta)

% z is a matrix with 2N rows and 14 columns; 

dim = size(z,1);
N = dim/2;
x = z(1:N,:)'; % col(x2,x3,...,xN+1)
p = z(N+1:2*N,:)'; % col(p1,p2,...,pN)

%%%%% Params used to calculate penalty %%%%%
A_pesk = params(23);
K_CO2 = params(30);
k_CO2 = params(31);
K_a1 = params(32);
K_a2 = params(33);
W = params(34);
M_O2 = params(16);
M_CO2 = params(17);
rho_O2 = params(18);
rho_CO2 = params(19);
MR_O2 = M_O2 + rho_O2*W;
MR_CO2 = M_CO2 + rho_CO2*W;

%%%%% Initialize penalty and cost %%%%%
G_O2 = 0;
G_CO2 = 0;
% J = zeros(14*N+14*N+2*N,1);
J = zeros(14*N+14*N+2*N,1);

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
    F = v-x(:,k)+h*(f(x(:,k),u)'+B*u);

    J(14*(k-1)+1:14*k)=F;
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
    J(14*N+14*(k-1)+1:14*N+14*k)=F;

end

%%%%%%% Compute penalty %%%%%%%%%%%%%%%
for k=1:N 

    P_as = x(1,k);
    P_vs = x(2,k);
    P_aCO2 = x(10,k);
    P_aO2 = x(11,k);
    C_vCO2 = x(12,k);
    C_vO2 = x(13,k);
    R_s = A_pesk * C_vO2;
        
    F_s = (1/R_s)*(P_as-P_vs);
    
    temp_O2 = (-MR_O2 + F_s * (K_a1*(1-exp(-K_a2*P_aO2))^2 - C_vO2))^2;
    temp_CO2 = (MR_CO2 + F_s * (K_CO2*P_aCO2+k_CO2-C_vCO2))^2;
    J(2*14*N+k)=alpha*temp_O2;
    J(2*14*N+N+k)=beta*temp_CO2;
end

J(1:14*N);
pause


% for i=1:N
%     J(2*14*N+N+N+i)=10*min(-p(9,i)/R(1,1),0)^2;
% end