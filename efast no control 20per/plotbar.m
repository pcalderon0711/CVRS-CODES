maxi = max(s_HIV.Si,[],2);
maxti = max(s_HIV.Sti,[],2);

[~,i] = sort(maxi,'ascend');
[~,ti] = sort(maxti,'ascend');

namesi = efast_var(i);
sortedi = s_HIV.Si(i,:);
bar(1:38,sortedi)
set(gca,'XTick',1:38)
set(gca,'XTickLabel',namesi)
set(gca,'XTickLabelRotation',45)
axis([0 39 0 0.5])
ylabel('eFast sensitivity')

pause
namesti = efast_var(ti);
sortedti = s_HIV.Sti(ti,:);
bar(1:38,sortedti)
set(gca,'XTick',1:38)
set(gca,'XTickLabel',namesti)
set(gca,'XTickLabelRotation',45)
axis([0 39 0 0.5])
ylabel('eFast sensitivity')