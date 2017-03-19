%maxi = max(s_HIV.Si, [], 2);
%maxti = max(s_HIV.Sti, [], 2);

rank = 1:38;

si4 = s_HIV.Si(:,1);
si6 = s_HIV.Si(:,2);
si12 = s_HIV.Si(:,3);

sti4 = s_HIV.Sti(:,1);
sti6 = s_HIV.Sti(:,2);
sti12 = s_HIV.Sti(:,3);

psi4 = s_HIV.p_Si(:,1);
psi4(38) = 1000;
psi6 = s_HIV.p_Si(:,2);
psi6(38) = 1000;
psi12 = s_HIV.p_Si(:,3);
psi12(38) = 1000;

psti4 = s_HIV.p_Sti(:,1);
psti4(38) = 1000;
psti6 = s_HIV.p_Sti(:,2);
psti6(38) = 1000;
psti12 = s_HIV.p_Sti(:,3);
psti12(38) = 1000;


[~,i4] = sort(si4,'descend');
[~,i6] = sort(si6,'descend');
[~,i12] = sort(si12,'descend');

[~,ti4] = sort(sti4,'descend');
[~,ti6] = sort(sti6,'descend');
[~,ti12] = sort(sti12,'descend');

namesi4 = efast_var(i4);
namesi6 = efast_var(i6);
namesi12 = efast_var(i12);

namesi4 = strcat('$',namesi4,'$');
namesi6 = strcat('$',namesi6,'$');
namesi12 = strcat('$',namesi12,'$');

namesti4 = efast_var(ti4);
namesti6 = efast_var(ti6);
namesti12 = efast_var(ti12);

namesti4 = strcat('$',namesti4,'$');
namesti6 = strcat('$',namesti6,'$');
namesti12 = strcat('$',namesti12,'$');

sort_si4 = si4(i4);
sort_si6 = si6(i6);
sort_si12 = si12(i12);

sort_sti4 = sti4(ti4);
sort_sti6 = sti6(ti6);
sort_sti12 = sti12(ti12);

sort_psi4 = psi4(i4);
sort_psi6 = psi6(i6);
sort_psi12 = psi12(i12);

sort_psti4 = psti4(ti4);
sort_psti6 = psti6(ti6);
sort_psti12 = psti12(ti12);

sig_si4 = sort_psi4 < 0.05;
sig_si6 = sort_psi6 < 0.05;
sig_si12 = sort_psi12 < 0.05;

sig_sti4 = sort_psti4 < 0.05;
sig_sti6 = sort_psti6 < 0.05;
sig_sti12 = sort_psti12 < 0.05;

tabi4 = table(rank', namesi4', sort_si4, sort_psi4, sig_si4);
tabi6 = table(rank',namesi6', sort_si6, sort_psi6, sig_si6);
tabi12 = table(rank', namesi12',sort_si12, sort_psi12, sig_si12);

tabti4 = table(rank',namesti4', sort_sti4, sort_psti4, sig_sti4);
tabti6 = table(rank', namesti6',sort_sti6, sort_psti6, sig_sti6);
tabti12 = table(rank', namesti12',sort_sti12, sort_psti12, sig_sti12);

writetable(tabi4);
writetable(tabi6);
writetable(tabi12);
writetable(tabti4);
writetable(tabti6);
writetable(tabti12);