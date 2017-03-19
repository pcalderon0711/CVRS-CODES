%function myMain()
clc
clf
clear
hold off

params_rest = myLoader('parameters_rest.txt','p');
params_exer = myLoader('parameters_exer.txt','p');
y0 = myEquilibriumSolver(params_rest,78.5,40);
y1 = myEquilibriumSolver(params_exer,107,40);
t1 = 5;
t2 = 15;

N_list = [150];

for i = 1:length(N_list)
    N = N_list(i);
    [T,Y] = myBWEulerSolver(@(t,y) myModel(t,y,params_rest,params_exer,t1),[0,t2], y0, N);
%     [T,Y] = ode15s(@(t,y) myModel(t,y,params_rest,params_exer,t1),[0,t2], y0, options);
    x=Y;
    x(50:end,9) = y1(9);
    x(50:end,14) = y1(14);
    x=x';
    save('150.mat','T','x','params_rest','params_exer');
end

% for i = 1:12
%     myPlotConfigurator(i,N_list);
% end