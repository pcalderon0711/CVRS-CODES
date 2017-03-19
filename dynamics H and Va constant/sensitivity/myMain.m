%function myMain()
clc
clf
clear
hold off

params_rest = myLoader('parameters_rest.txt','p');
params_exer = myLoader('parameters_exer.txt','p');
y0 = myEquilibriumSolver(params_rest,78.5,40);
y0 = [y0; zeros(12*36,1)];
t1 = 5;
t2 = 15;

reltol = 6;
options = odeset('RelTol',reltol);
[T,Y] = ode15s(@(t,y) mySensitivityModel(t,y,params_rest,params_exer,t1),[0,t2], y0, options);

% Here we make the sensitivities nondimensional by multiplying by the
% parameter value and dividing by the state value
YS = zeros(size(Y));
for i=1:length(T)
    t = T(i);
    x = Y(i,1:12);
    z = reshape(Y(i,13:end),12,36);
    for j=1:12
        % if the state value is 0, we just set the sensitivity to 0
        if x(j) == 0
            z(j,:) = 0;
        else
            z(j,:) = z(j,:)/x(j);
        end
    end
    for j=1:36
        % the parameter value we multiply is dependent on whether the
        % subject is in the rest or exercise regime
        if t<t1
            params = params_rest;
        else
            params = params_exer;
        end
        % if parameter value is 0, we just set sensitivity to 0
        if params(j) == 0
            z(:,j) = 0;
        else
            z(:,j) = z(:,j)*params(j);
        end
    end
    YS(i,1:12) = x;
    YS(i,13:end) = reshape(z,12*36,1);
end

save(sprintf('%d.mat',1),'reltol','T','YS');

% indicate the top x params which we are going to plot
top = 5;
leg = {'c_{as}','c_{vs}','c_{ap}','c_{vp}','c_l','c_r','R_l','R_r','\kappa','\alpha_l','\alpha_r','\beta_l','\beta_r','\gamma_l','\gamma_r','M_{O2}','M_{CO2}','\rho_{O2}','\rho_{CO2}','q_{as}','V_{tot}','R_p','A_{pesk}','P_{IO2}','P_{ICO2}','V_{AO2}','V_{ACO2}','V_{TO2}','V_{TCO2}','K_{CO2}', 'k_{CO2}', 'K_{a1}', 'K_{a2}', 'W', 'H','dot_{VA}'};
for i=1:12
    magnitudes = zeros(36,1);
    % calculate the mean of the sensitivity of each parameter over the
    % whole time interval
    for j = 1:36
        % to account for negative sensitivity, we calculate mean of
        % absolute values
        magnitudes(j) = mean(abs(YS(:,12*j+i)));
    end
    % store the ranking of sensitivities in top_indices
    [~,top_indices] = sort(magnitudes,'descend');
    myPlotConfigurator(i,leg,top,top_indices)
end



% Read Y as:
% Y(1:12) = state variables
% Y(13:24) = sensitivities of states wrt first parameter
% Y(25:36) = sensitivityies of states wrt 2nd parameter...