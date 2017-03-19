maxi = max(s_HIV.Si,[],2);
maxti = max(s_HIV.Sti,[],2);

[~,i] = sort(maxi,'ascend');
[~,ti] = sort(maxti,'ascend');

namesi = efast_var(i);
sortedi = s_HIV.Si(i,:);
included=setdiff(1:38,[12,14,36,37]);
bar(1:34,sortedi(included,:))
set(gca,'XTick',1:34)
set(gca,'XTickLabel',namesi(included))
set(gca,'XTickLabelRotation',45)
axis([0 35 0 0.5])
ylabel('eFast sensitivity')

pause
namesti = efast_var(ti);
sortedti = s_HIV.Sti(ti,:);
included=setdiff(1:38,[11,12,36,37]);
bar(1:34,sortedti(included,:))
set(gca,'XTick',1:34)
set(gca,'XTickLabel',namesti(included))
set(gca,'XTickLabelRotation',45)
axis([0 35 0 0.5])
ylabel('eFast sensitivity')