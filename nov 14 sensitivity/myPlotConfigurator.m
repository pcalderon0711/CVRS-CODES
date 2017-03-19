function myPlotConfigurator(v,leg,top,ranked_indices,t,dx_dmu)

hold on
for i = 1:top
    x = plot(t,dx_dmu(:,v,ranked_indices(i)),'DisplayName',leg{ranked_indices(i)});
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
        label='H';
    case 10
        label='P_{aCO2}';
    case 11
        label='P_{aO2}';
    case 12
        label='C_{vCO2}';
    case 13
        label='C_{vO2}';
    case 14
        label='dot_VA';
end

xlabel('t');
ylabel(label);
set(gca, 'FontSize', 15)
legend('show')

savefig(sprintf('%d',v));
print(sprintf('%d',v),'-dpng');

clf
end