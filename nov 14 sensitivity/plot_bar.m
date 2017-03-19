% var={'c_{as}', 'c_{vs}', 'c_{ap}', 'c_{vp}', ...
% 'c_l','c_r', 'R_l','R_r', '\kappa','\alpha_l', ...
% '\alpha_r', '\beta_l', '\beta_r', '\gamma_l', ...
% '\gamma_r', 'M_{O2}', 'M_{CO2}', '\rho_{O2}', ...
% '\rho_{CO2}', 'R_{p,rest}', 'A_{pesk,rest}', ...
% 'P_{IO2}', 'P_{ICO2}', 'V_{AO2}', 'V_{ACO2}',...
% 'V_{TO2}', 'V_{TCO2}', 'K_{CO2}', 'k_{CO2}', ...
% 'K_{a1}', 'K_{a2}', 'A_{pesk,exer}', 'R_{p,exer}', ...
% 'H_{rest}', 'H_{exer}', 'dot_{VA,rest}', 'dot_{VA,exer}'};%
% var{37}='W_{exer}'
% var{36}='R_{p,exer}'
% var{35}='A_{pesk,exer}'
% var{20}='q_{as}'
% var{21}='V_{tot}'
% var{22}='R_{p,rest}'
% var{23}='A_{pesk,rest}'
% var{24}='P_{I,O2}'
% var{25}='P_{I,CO2}'
% var{26}='V_{A,O2}'
% var{27}='V_{A,CO2}'
% var{28}='V_{T,O2}'
% var{29}='V_{T,CO2}'
% var{30}='K_{CO2}'
% var{30}='K_{O2}'
% var{31}='K_{CO2}'
% var{32}='K_{a1}'
% var{33}='K_{a2}'
% var{34}='W_{rest}'

top = 15;

idx = ranked_indices(1,:);
val = ranked_values(1,:);
included=setdiff(1:37,[20,21,34,37]);
x = var(idx(included));
y = val(included);
bar(1:top,fliplr(y(1:top)))
set(gca,'XTick',1:top)
set(gca,'XTickLabel',fliplr(x(1:top)))
set(gca,'XTickLabelRotation',45)
axis([0 top+1 0 1.4])
