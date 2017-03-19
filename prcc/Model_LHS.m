%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
clear all;
close all;

%% Sample size N
runs=100;

%% LHS MATRIX  %%
Parameter_settings_LHS;

c_as_LHS=LHS_Call(1e-2, c_as, 2e-2, 2e-3,runs,'norm'); % baseline = 10
c_vs_LHS=LHS_Call(0.4, c_vs, 0.6, 0.05 ,runs,'norm'); % baseline = 2e-2
c_ap_LHS=LHS_Call(0.02, c_ap, 0.05, 0.005, runs,'norm'); % baseline = 3e-2
c_vp_LHS=LHS_Call(0.1, c_vp, 0.2, 0.05 ,runs,'norm'); % baseline = 2.4e-5
c_l_LHS=LHS_Call(0.01, c_l, 0.03, 0.005, runs,'norm'); % baseline = 3e-3
c_r_LHS= LHS_Call(0.03 , c_r , 0.05 , 0.005 , runs , 'norm');  % baseline = 0.24
R_l_LHS=LHS_Call(0.1,R_l,0.3, 0.05 ,runs,'norm'); % baseline = 1200
R_r_LHS=LHS_Call(0,R_r,0.1, 0.01 ,runs,'norm'); % baseline = 2.4
kappa_LHS=LHS_Call(0.03,kappa,0.07, 0.005 ,runs,'norm'); % dummy parameter
alpha_l_LHS=LHS_Call(25, alpha_l, 35, 5 ,runs,'norm'); % baseline = 10
alpha_r_LHS=LHS_Call(25, alpha_r, 35, 5 ,runs,'norm'); % baseline = 2e-2
beta_l_LHS=LHS_Call(20, beta_l, 30, 5, runs,'norm'); % baseline = 3e-2
beta_r_LHS=LHS_Call(1, beta_r, 2, 0.5 ,runs,'norm'); % baseline = 2.4e-5
gamma_l_LHS=LHS_Call(-2, gamma_l, -1, 0.5, runs,'norm'); % baseline = 3e-3
gamma_r_LHS= LHS_Call(-2, gamma_r , -1 , 0.5 , runs , 'norm');  % baseline = 0.24
M_O2_LHS=LHS_Call(0.1, M_O2,1, 0.05 ,runs,'norm'); % baseline = 1200
M_CO2_LHS=LHS_Call(0.1,M_CO2,1, 0.05 ,runs,'norm'); % baseline = 2.4
rho_O2_LHS=LHS_Call(0.005,rho_O2,0.015, 0.0005 ,runs,'norm'); % dummy parameter
rho_CO2_LHS=LHS_Call(0.005,rho_CO2, 0.015, 0.0005 ,runs,'norm'); % baseline = 10
R_p_rest_LHS=LHS_Call(1.2,R_p_rest,1.8, 0.1 ,runs,'norm'); % baseline = 2.4e-5
A_pesk_rest_LHS=LHS_Call(160, A_pesk_rest, 190, 10, runs,'norm'); % baseline = 3e-3
P_IO2_LHS= LHS_Call(140 , P_IO2 , 160 , 10 , runs , 'norm');  % baseline = 0.24
P_ICO2_LHS=LHS_Call(0,P_ICO2,1, 0.05 ,runs,'norm'); % baseline = 1200
V_AO2_LHS=LHS_Call(2,V_AO2,3, 0.5 ,runs,'norm'); % baseline = 2.4
V_ACO2_LHS=LHS_Call(2,V_ACO2,4, 0.5 ,runs,'norm'); % dummy parameter
V_TO2_LHS=LHS_Call(4, V_TO2, 8, 1 ,runs,'norm'); % baseline = 10
V_TCO2_LHS=LHS_Call(13, V_TCO2, 17, 5 ,runs,'norm'); % baseline = 2e-2
K_CO2_LHS=LHS_Call(0.001, K_CO2, 0.01, 0.05, runs,'norm'); % baseline = 3e-2
k_CO2_LHS=LHS_Call(0.2,k_CO2,0.3, 0.05 ,runs,'norm'); % baseline = 2.4e-5
K_a1_LHS=LHS_Call(0.1, K_a1, 0.3, 0.05, runs,'norm'); % baseline = 3e-3
K_a2_LHS= LHS_Call(0 , K_a2, 0.1 , 0.05 , runs , 'norm');  % baseline = 0.24
A_pesk_exer_LHS=LHS_Call(250,A_pesk_exer,300, 10 ,runs,'norm'); % baseline = 2.4
R_p_exer_LHS=LHS_Call(0.1,R_p_exer,0.5, 0.05 ,runs,'norm'); % dummy parameter


%% LHS MATRIX and PARAMETER LABELS
LHSmatrix=[c_as_LHS c_vs_LHS c_ap_LHS c_vp_LHS c_l_LHS c_r_LHS ...
              R_l_LHS R_r_LHS kappa_LHS  alpha_l_LHS alpha_r_LHS ...
              beta_l_LHS beta_r_LHS gamma_l_LHS gamma_r_LHS M_O2_LHS M_CO2_LHS ...
              rho_O2_LHS rho_CO2_LHS R_p_rest_LHS A_pesk_rest_LHS P_IO2_LHS ...
              P_ICO2_LHS V_AO2_LHS V_ACO2_LHS V_TO2_LHS V_TCO2_LHS K_CO2_LHS ...
              k_CO2_LHS K_a1_LHS K_a2_LHS A_pesk_exer_LHS R_p_exer_LHS];

for i=1:100
i
ff=myEquilibriumSolver(LHSmatrix,i,78.5,40,0);
ff(1)
end

load('flatdata.mat');
ts = 5;
opts = optimset('Diagnostics','off', 'Display','off');
h = 0.1;

for r=1:runs %Run solution x times choosing different values
    f=@ODE_LHS;
    r
    LHSmatrix(r,:);
    A = zeros(length(tspan),14);
    A(1,:) = myEquilibriumSolver(LHSmatrix,r,78.5,40,0);
    for index=1:length(tspan)-1
        f = @ODE_LHS;
        if index == 1
            prev = A(1,:);
        end
        index
        implicit = @(next)(next - prev - h*f(tspan(index+1),next,u, LHSmatrix,r,tspan,ts));
        A(index+1,:) = fsolve(implicit, prev, opts);
        prev = A(index+1,:);
    end
        
    
%     [t,y]=ode15s(@(t,y)f(t,y,LHSmatrix,x,runs),tspan,y0,[]); 
%     A=[t y]; % [time y]

    %% Save the outputs at ALL time points [tspan]
    %T_lhs(:,x)=Anew(:,1);
    %CD4_lhs(:,x)=Anew(:,2);
    %T1_lhs(:,x)=Anew(:,3);
    %T2_lhs(:,x)=Anew(:,4);
    %V_lhs(:,x)=Anew(:,5);
    
    %% Save only the outputs at the time points of interest [time_points]:
    %% MORE EFFICIENT
    Pas_lhs(:,r)=A(time_points+1,1);
    Paco2_lhs(:,r)=A(time_points+1,10);
end
%% Save the workspace
save Model_LHS.mat;
% CALCULATE PRCC
[prcc sign sign_label]=PRCC(LHSmatrix,Pas_lhs,time_points,PRCC_var,0.05);
