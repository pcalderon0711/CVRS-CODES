% load('constant_constant_control(N150.000000,a20.000000,b5.000000,T15.000000).mat');
load('150.mat')
params = params_rest;
params(35) = params_exer(23); %Apesk exer
params(36) = params_exer(22); %Rp exer
params(37) = params_exer(34); %W exer

sigma0 = 0.5;
p = 5;
% 
% %%%%%
% %replace one at a time
% at_a_time = 1;
% combinations = [16,18,23,32,35];
% replacements = [9,22,36];
% 
% for i=1:length(combinations)
%     for j=1:length(replacements)
%         combinations((i-1)*3+j+1,:) = combinations(1,:);
%         combinations((i-1)*3+j+1,i) = replacements(j);
%     end
% end
% %%%%%
% 
% 
% % p = 1;
% % combinations = combnk([1:19,22:37], p); % do not include 20,21 qas Vtot
% 
% valid = [];
% condition = [];
% selection = [];
% 
% numComb = size(combinations,1);
% for i=1:numComb
%     i
%     log = zeros(37, 1);
%     for j=1:37
%         if ismember(j, combinations(i,:))
%             log(j) = 1;
%         end
%     end
%     dPasdmu = mySensitivityAnalysis(logical(log));
%     % calculate rank
%     if rank(dPasdmu'*dPasdmu) == p
%         valid(end+1,:) = combinations(i,:);
%         condition(end+1) = cond(dPasdmu'*dPasdmu);
%         sigma = (sigma0^2)*inv(dPasdmu'*dPasdmu);
%         nu = zeros(p,1);
%         for j=1:p
%             if params(combinations(i,j))==0
%                 nu(j) = 0;
%             else
%                 nu(j) = sqrt(sigma(j,j))/params(combinations(i,j));
%             end
%         end
%         selection(end+1) = norm(nu);
%     end
% end
% 
% save(sprintf('subset_%d_%d_at_a_time_CONSTANT.mat',p,at_a_time),'p','valid','condition','selection','dPasdmu');
% 
% %%%%%
% %replace two at a time
% at_a_time = 2;
% combinations = [16,18,23,32,35];
% candidates = [9,22,36];
% replacements = combnk(candidates,2);
% to_replace = combnk(combinations(1,:),2);
% for i=1:length(to_replace)
%     for j=1:length(replacements)
%         combinations((i-1)*3+j+1,:) = combinations(1,:);
%         for k=1:2
%             rep_idx = (combinations((i-1)*3+j+1,:) == to_replace(i,k));
%             combinations((i-1)*3+j+1,rep_idx) = replacements(j,k);
%         end
%     end
% end
% %%%%
% 
% valid = [];
% condition = [];
% selection = [];
% 
% numComb = size(combinations,1);
% for i=1:numComb
%     i
%     log = zeros(37, 1);
%     for j=1:37
%         if ismember(j, combinations(i,:))
%             log(j) = 1;
%         end
%     end
%     dPasdmu = mySensitivityAnalysis(logical(log));
%     % calculate rank
%     if rank(dPasdmu'*dPasdmu) == p
%         valid(end+1,:) = combinations(i,:);
%         condition(end+1) = cond(dPasdmu'*dPasdmu);
%         sigma = (sigma0^2)*inv(dPasdmu'*dPasdmu);
%         nu = zeros(p,1);
%         for j=1:p
%             if params(combinations(i,j))==0
%                 nu(j) = 0;
%             else
%                 nu(j) = sqrt(sigma(j,j))/params(combinations(i,j));
%             end
%         end
%         selection(end+1) = norm(nu);
%     end
% end
% 
% save(sprintf('subset_%d_%d_at_a_time_CONSTANT.mat',p,at_a_time),'p','valid','condition','selection','dPasdmu');
% 




combinations=[[16,18,23,32,35]];
valid = [];
condition = [];
selection = [];

numComb = size(combinations,1);
for i=1:numComb
    i
    log = zeros(37, 1);
    for j=1:37
        if ismember(j, combinations(i,:))
            log(j) = 1;
        end
    end
    dPasdmu = mySensitivityAnalysis(logical(log));
%     calculate rank
    if rank(dPasdmu'*dPasdmu) == p
        valid(end+1,:) = combinations(i,:);
        condition(end+1) = cond(dPasdmu'*dPasdmu);
        sigma = (sigma0^2)*inv(dPasdmu'*dPasdmu);
        nu = zeros(p,1);
        for j=1:p
            if params(combinations(i,j))==0
                nu(j) = 0;
            else
                nu(j) = sqrt(sigma(j,j))/params(combinations(i,j));
            end
        end
        selection(end+1) = norm(nu);
    end
end

save(sprintf('sensitive_set_CONSTANT.mat'),'p','valid','condition','selection','dPasdmu','nu');
