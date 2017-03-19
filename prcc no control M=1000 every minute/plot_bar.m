prcc = prcc'; %turn into 37*3
maxvals = max(prcc,[],2);
[m,i] = sort(maxvals,'ascend');
names = PRCC_var(i);
sorted = prcc(i,:);
% bar(1:37,sorted)
included=setdiff(1:37,[24,25,35,37]);
bar(1:33,sorted(included,:))
set(gca,'XTick',1:33)
set(gca,'XTickLabel',names(included))
set(gca,'XTickLabelRotation',45)
axis([0 34 -1 1])
