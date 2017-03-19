% load('subset_5_1_at_a_time.mat')
load('subset_5_1_at_a_time_CONSTANT.mat')
valid1 = valid;
condition1 = condition';
selection1 = selection';
% load('subset_5_2_at_a_time.mat')
load('subset_5_2_at_a_time_CONSTANT.mat')
valid2 = valid(2:end,1:end);
condition2 = condition';
selection2 = selection';
condition2 = condition2(2:end);
selection2 = selection2(2:end);

selection = [selection1;selection2];
[sorted,I] = sort(selection);
selection(I);
pause

valid = [valid1;valid2]
condition = [condition1;condition2]

valid = valid(I,:);
condition = condition(I);
selection = selection(I);


labels = {}
for i=1:size(valid,1)
    for j=1:size(valid,2)
        c = valid(i,j);
        if c == 9
            r = '$\kappa$';
        elseif c == 22
            r = '$r_{p,rest}$';
        elseif c == 36
            r = '$r_{p,exer}$';
        elseif c == 23
            r = '$A_{pesk,rest}$';
        elseif c == 35
            r = '$A_{pesk,exer}$';
        elseif c == 18
            r = '$\rho_{O2}$';
        elseif c == 32
            r = '$K_{a1}$';
        else
            r = '$M_{O2}$';
        end
        labels{i,j} = r;
    end
end
f =1:length(labels)';


tab_CONSTANT = table(f',labels,condition,selection);
writetable(tab_CONSTANT);