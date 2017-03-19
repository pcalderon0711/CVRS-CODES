load('sensitive_set_CONSTANT.mat');

% load('constant_constant_control(N150.000000,a20.000000,b5.000000,T15.000000).mat');
% load('sensitive_set.mat');

combination = [16,18,23,32,35];
%16: mo2, 18: rho2, 23: apeskrest, 32:ka1, 35:apeskexer

params = params_rest;
params(35) = params_exer(23); %Apesk exer
params(36) = params_exer(22); %Rp exer
params(37) = params_exer(34); %W exer

% sd for rest, exercise, average of rest and exer
kappel = [4.580968555327646, 4.761612614625580, 4.671290584976614];
peer = [3.426232911621728, 7.492192435571742, 5.459212673596735];
king = [4.816532069224170, 8.461770391146112, 6.639151230185141];

kappel_pm = zeros(5,3);
peer_pm = zeros(5,3);
king_pm = zeros(5,3);

for i=1:3
    sigma_kappel = kappel(i)^2*inv(dPasdmu'*dPasdmu);
    sigma_peer = peer(i)^2*inv(dPasdmu'*dPasdmu);
    sigma_king = king(i)^2*inv(dPasdmu'*dPasdmu);
    
    for j=1:p
        kappel_pm(j,i) = sqrt(sigma_kappel(j,j));
        peer_pm(j,i) = sqrt(sigma_peer(j,j));
        king_pm(j,i) = sqrt(sigma_king(j,j));
    end
end 

names = {'$M_{O2}$','$\rho_{O2}$','$A_{pesk,rest}$','$K_{a1}$','$A_{pesk,exer}$'};

% tabkappel = table(names', kappel_pm);
% tabpeer = table(names', peer_pm);
% tabking = table(names', king_pm);
% 
% writetable(tabkappel);
% writetable(tabpeer);
% writetable(tabking);
% 
tabkappel_CONSTANT = table(names', kappel_pm);
tabpeer_CONSTANT = table(names', peer_pm);
tabking_CONSTANT = table(names', king_pm);

writetable(tabkappel_CONSTANT);
writetable(tabpeer_CONSTANT);
writetable(tabking_CONSTANT);


%variance
corrected_kappel = 37.1234;
corrected_peer = 17.3781;
corrected_king = 16.5771;

ckappelsigma = corrected_kappel*inv(dPasdmu'*dPasdmu);
cpeersigma = corrected_peer*inv(dPasdmu'*dPasdmu);
ckingsigma = corrected_king*inv(dPasdmu'*dPasdmu);


for j=1:p
    corrected_kappel_pm(j) = sqrt(ckappelsigma(j,j));
    corrected_peer_pm(j) = sqrt(cpeersigma(j,j));
    corrected_king_pm(j) = sqrt(ckingsigma(j,j));
end

% tabcorrected = table(names', corrected_kappel_pm', corrected_peer_pm', corrected_king_pm');
% writetable(tabcorrected);

tabcorrected_CONSTANT = table(names', corrected_kappel_pm', corrected_peer_pm', corrected_king_pm');
writetable(tabcorrected_CONSTANT);


save('computed_uncertainty_CONSTANT.mat','kappel','peer','king','kappel_pm','peer_pm','king_pm','corrected_kappel_pm','corrected_peer_pm','corrected_king_pm');
% save('computed_uncertainty.mat','kappel','peer','king','kappel_pm','peer_pm','king_pm','corrected_kappel_pm','corrected_peer_pm','corrected_king_pm');


% standard dev * (    0.1487
%    0.0021
%    1.7745
%    0.0364
%   17.0585