function myOptimalControlPlotterTimeVarying(t,T,u,x,p,N,alpha,beta,prefix,w)

cd('plots'); % go to the 'plots' folder directory

u1 = u(1,:);
u2 = u(2,:);
q = plot(t,u1,t,u2);
legend('u1','u2');
title('Approximation for the controls');
xlabel('time (min)');
% ylim([0,10]);
q(1).LineWidth = 3;
q(2).LineWidth = 3;
ax = gca;
ax.XTick = 0:1:15;
ax.FontSize = 17;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_control',prefix, w, N, T, alpha, beta),'-dpng');
close;

q=plot(t,x(1,:));
title('Approximation for P_{as}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([100,130]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = [100, 105, 110, 115, 120, 125, 130];
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_P_as',prefix, w, N, T, alpha, beta),'-dpng');
close;

q=plot(t,x(2,:));
title('Approximation for P_{vs}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([4,4.6]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = [4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6];
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_P_vs',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(3,:));
title('Approximation for P_{ap}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([9,13]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 9:0.5:13;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_P_ap',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(4,:));
title('Approximation for P_{vp}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([5.5,8]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 5.5:0.5:8;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_P_vp',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(9,:));
title('Approximation for H');
xlabel('time (min)');
ylabel('1/min');
% ylim([77,85]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 77:85;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_H',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(10,:));
title('Approximation for P_{a,CO2}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([39.9,40.9]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 39.9:0.1:40.9;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_PaCO2',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(11,:));
title('Approximation for P_{a,O2}');
xlabel('time (min)');
ylabel('mmHg');
% ylim([90.5,105]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 90:5:110;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_PaO2',prefix, w,N,T, alpha, beta),'-dpng');
close;

q=plot(t,x(14,:));
title('Approximation for dotV_{A}');
xlabel('time (min)');
ylabel('liters/min');
% ylim([6,20]);
q.LineWidth = 3;
ax = gca;
ax.FontSize = 17;
ax.XTick = 0:1:15;
% ax.YTick = 6:2:20;
print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_dotVA',prefix, w,N, T,alpha, beta),'-dpng');
close;

for i=1:14
    q=plot(t,p(i,:));
    title(sprintf('Approximation for p_{%d}',i));
    xlabel('time (min)');
    % ylim([6,20]);
    q(1).LineWidth = 3;
    ax = gca;
    ax.FontSize = 17;
    ax.XTick = 0:1:15;
    % ax.YTick = 6:2:20;
    print(sprintf('%s_%s_N%d_T%d_alpha%d_beta%d_p%d',prefix,w,N,T, alpha, beta, i),'-dpng');
    close;
end

cd ..;
