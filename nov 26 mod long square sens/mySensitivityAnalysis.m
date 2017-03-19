function [dx_dmu, ranked_values, ranked_indices] = mySensitivityAnalysis()

load('timevar_longsquare_control(N150.000000,a20.000000,b5.000000,T15.000000).mat');
h = T/N;

params = params_rest;
params(35) = params_exer(23); %Apesk exer
params(36) = params_exer(22); %Rp exer
params(37) = params_exer(34); %W exer
len_x = size(x,1);
len_params = size(params,1);
dx_dmu = zeros(N+1,len_x,len_params);

opts = optimset('Diagnostics','off', 'Display','off');

dS = @(ix,xm)(mySensitivitySystem(ix, xm, x(:,ix), u(:,ix), params,ts,T,N));

% fill up dx_dmu for all times
for index=1:N %t = 0 to t=10, index=1 is 0, index=N+1 is T
    index
%    dS = mySensitivitySystem(index+1,dx_dmu,x,u,params_rest,params_exer,ts,T,N);
    % this is dS evaluated at the next time step
    dx = reshape(dx_dmu(index,:,:),len_x,len_params);
    implicit = @(next)(next - dx - h*dS(index+1,next));
    dx_dmu(index+1,:,:) = fsolve(implicit, dx, opts);
    
%     for i=1:len_x
%         dx = reshape(dx_dmu(index,i,:),1,len_params);
%         dS(next,index+1)
%         implicit = @(next)(next - dx - h*dS(i,:));
%         dx_dmu(index+1,i,:) = fsolve(implicit, dx, opts);
%         for j=1:len_params
%       dx_dmu(index+1,i,j) = dx_dmu(index,i,j) + h*dS(i,j);
%     end
end

for index=1:N+1 %t = 0 to t=10, index=1 is 0, index=N+1 is T
    for i=1:len_x
        for j=1:len_params
            dx_dmu(index,i,j) = dx_dmu(index,i,j) * (params(j)/x(i,index)); % nondimensional sensitivity
        end
    end
end

ranked_values = zeros(len_x,len_params);
ranked_indices = zeros(len_x,len_params);
for i=1:len_x
    max_list = zeros(len_params,1);
    for j=1:len_params
        max_list(j) = max(abs(dx_dmu(:,i,j))); %max over all time
    end
    [temp_vals,temp_indices] = sort(max_list,'descend');
    ranked_values(i,:) = temp_vals;
    ranked_indices(i,:) = temp_indices;
end

save('sensitivity_data.mat','dx_dmu', 'ranked_values', 'ranked_indices')
% ranked_values: columns are states, rows are parameters, rows are sorted
% from highest sensitivity to lowest, from top to bottom

top = 5;
leg = {'c_{as}','c_{vs}','c_{ap}','c_{vp}','c_l','c_r','R_l','R_r','\kappa','\alpha_l','\alpha_r','\beta_l','\beta_r','\gamma_l','\gamma_r','M_{O2}','M_{CO2}','\rho_{O2}','\rho_{CO2}','q_{as}','V_{tot}','R_p_rest','A_{pesk}_rest','P_{IO2}','P_{ICO2}','V_{AO2}','V_{ACO2}','V_{TO2}','V_{TCO2}','K_{CO2}', 'k_{CO2}', 'K_{a1}', 'K_{a2}', 'W_rest','A_pesk_exer','R_p_exer','W_exer'};
for i=1:14
    myPlotConfigurator(i,leg,top,ranked_indices(i,:),t,dx_dmu);
end

