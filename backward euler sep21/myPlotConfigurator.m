function myPlotConfigurator(v,leg)
hold on
for i = 1:length(leg)
    load(sprintf('%d.mat',i));
    x = plot(T,Y(:,v),'DisplayName',num2str(leg(i)));
    x.LineWidth = 4;
end
hold off

switch v
    case 1
        label='P_{as}';
        ylim([0, 200]);
    case 2
        label='P_{vs}';
        ylim([0,10]);
    case 3
        label='P_{ap}';
        ylim([0,20]);
    case 4
        label='P_{vp}';
        ylim([0,15]);
    case 5
        label='S_l';
    case 6
        label='\sigma_l';
    case 7
        label='S_r';
    case 8
        label='\sigma_r';
    case 9
        label='P_{aCO2}';
        ylim([0,60]);
    case 10
        label='P_{aO2}';
        ylim([0,150]);
    case 11
        label='C_{vCO2}';
        ylim([0,1.5]);
    case 12
        label='C_{vO2}';
        ylim([0,1.5]);
end

xlabel('t');
ylabel(label);
set(gca, 'FontSize', 15)
legend('show')

savefig(sprintf('%d',v));
print(sprintf('%d',v),'-dpng');

clf
end