%function myMain()
clc
clf
clear
hold off

params_rest = myLoader('parameters_rest.txt','p');
params_exer = myLoader('parameters_exer.txt','p');
y0 = myEquilibriumSolver(params_rest,78.5,40);;

H = y0(9);
dotV_A = y0(14);

terminal = 10;
count=1;
mode = 'single';
% [T,Y] = ode15s(@(t,y) myModel(t,y,params_rest),[0 terminal], y0);
mod_y0 = myModifyIc(y0,H-0.5,dotV_A);
[T,Y] = ode15s(@(t,y) myModel(t,y,params_rest),[0 terminal], mod_y0);

cd base
for v=1:14
    legend_info{v,count}=sprintf('(%.2f,%.2f)',H,dotV_A);
    myPlotConfigurator;
end
save(sprintf('mod(%.2f,%.2f).mat',H,dotV_A),'T','Y');
cd ..
% terminal = 10;
% cd ../vary_H
% count = 1;
% for i=[0,5]
%     for j=[0]
%         mode = 'minus';
%         mod_y0 = myModifyIc(y0,H-i,dotV_A-j);
%         [T,Y] = ode15s(@(t,y) myModel(t,y,params),[0 terminal], mod_y0);
%         i,j
%         save(sprintf('(%.2f,%.2f).mat',H-i,dotV_A-j),'T','Y');
%         for v=1:14
%             legend_info{v,count}=sprintf('(%.2f,%.2f)',H-i,dotV_A-j);
%             myPlotConfigurator;
%         end
%         if count ~= 1
%             count = count + 1;
%             mode = 'plus';
%             mod = myModifyIc(params,H+i,dotV_A+j);
%             [T,Y] = ode45(@(t,y) myModel(t,y,mod),[0 terminal], y0);
%             save(sprintf('(%.2f,%.2f).mat',H+i,dotV_A+j),'T','Y');
%             for v=1:14
%                 legend_info{v,count}=sprintf('(%.2f,%.2f)',H+i,dotV_A+j);
%                 myPlotConfigurator;
%             end
%         end
%         count = count + 1;
%     end
% end
% 
% cd ../vary_dotV_A
% count = 1;
% for i=[0]
%     for j=[0,7.6]
%         mode = 'minus';
%         mod_y0 = myModifyIc(y0,H-i,dotV_A-j);
%         [T,Y] = ode45(@(t,y) myModel(t,y,params),[0 terminal], mod_y0);
%         i,j
%         save(sprintf('(%.2f,%.2f).mat',H-i,dotV_A-j),'T','Y');
%         for v=1:14
%             legend_info{v,count}=sprintf('(%.2f,%.2f)',H-i,dotV_A-j);
%             myPlotConfigurator;
%         end
%         if count ~= 1
%             count = count + 1;
%             mode = 'plus';
%             mod = myModifyIc(params,H+i,dotV_A+j);
%             [T,Y] = ode45(@(t,y) myModel(t,y,mod),[0 terminal], y0);
%             save(sprintf('(%.2f,%.2f).mat',H+i,dotV_A+j),'T','Y');
%             for v=1:14
%                 legend_info{v,count}=sprintf('(%.2f,%.2f)',H+i,dotV_A+j);
%                 myPlotConfigurator;
%             end
%         end
%         count = count + 1;
%     end
% end
% 
