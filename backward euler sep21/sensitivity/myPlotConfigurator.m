function myPlotConfigurator(v,leg,top,top_indices)
load('1.mat');
hold on
for i = 1:top
    x = plot(T,YS(:,12*top_indices(i)+v),'DisplayName',leg{top_indices(i)});
    x.LineWidth = 4;
end
hold off

switch v
    case 1
        label='P_{as}';
    case 2
        label='P_{vs}';
    case 3
        label='P_{ap}';
    case 4
        label='P_{vp}';
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
    case 10
        label='P_{aO2}';
    case 11
        label='C_{vCO2}';
    case 12
        label='C_{vO2}';
end

xlabel('t');
ylabel(label);
set(gca, 'FontSize', 15)
legend('show')

savefig(sprintf('%d',v));
print(sprintf('%d',v),'-dpng');

clf
end