prcc = prcc'; %turn into 37*3
maxvals = max(prcc,[],2);
[m,i] = sort(maxvals,'ascend');
names = PRCC_var(i);
sorted = prcc(i,:);
bar(1:37,sorted)
set(gca,'XTick',1:37)
set(gca,'XTickLabel',names)
set(gca,'XTickLabelRotation',45)
axis([0 38 -1 1])
