figure(v);
if count == 1
    plot(T,Y(:,v),'r');
elseif count == 2
    plot(T,Y(:,v),'gx');
elseif count == 3
    plot(T,Y(:,v),'b:');
elseif count == 4
    plot(T,Y(:,v),'c');
elseif count == 5
    plot(T,Y(:,v),'y');
end
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
        label='dotV_A';
end

xlabel('t');
ylabel(label);
ylim([min(Y(:,v))-0.5,max(Y(:,v))+0.5])

if strcmp(mode,'single') || (i==5 && j==0 && strcmp(mode,'plus')) || (i==0 && j==7.6 && strcmp(mode,'plus'))
    legend(legend_info{v,:});
    print(sprintf('mod%d',v),'-dpng');
end

if ~strcmp(mode,'single')
    hold on
end