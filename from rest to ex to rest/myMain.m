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
t3 = 25;
T = 30;

N_list = [300];

for i = 1:length(N_list)
    N = N_list(i);
    [T,Y] = myBWEulerSolver(@(t,y) myModel(t,y,params_rest,params_exer,t1,t2,t3),[0,T], y0, N);
%     [T,Y] = ode15s(@(t,y) myModel(t,y,params_rest,params_exer,t1),[0,t2], y0, options);
    save(sprintf('%d.mat',i),'N','T','Y');
end

for i = 1:12
    myPlotConfigurator(i,N_list);
end