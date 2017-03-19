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

span = [0.001,0.01,0.1,0.2,0.3,0.4,0.5,1];

for i = 1:length(span)
    reltol = span(i);
    options = odeset('RelTol', reltol);
    mode = 'single';
    [T,Y] = ode15s(@(t,y) myModel(t,y,params_rest,params_exer,t1),[0,t2], y0, options);
    save(sprintf('%d.mat',i),'reltol','T','Y');
end


for i = 1:12
    myPlotConfigurator(i,span);
end