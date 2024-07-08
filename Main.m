clearvars, close all, clc
%% Reaction study
Tvalues = (600:50:750)+273.15;                                              %[K]
V_BIG = 10000;                                                              % [m3]
Fin_h = 2000;                                                               % [kmol/h]
Fin_m = 100;                                                                % [kmol/h]
Fin_t = 400;                                                                % [kmol/h]
Fin_b = 0;                                                                  % [kmol/h]
Fin_d = 0;                                                                  % [kmol/h]
Fin_TOT = [Fin_h;Fin_m;Fin_t;Fin_b;Fin_d];                                  % [kmol/h]

for i = 1 : length(Tvalues)
    T = Tvalues(i);
    [Vout,Fout] = ode45(@(V,F)Kinetics(V,F,T), [0:1:V_BIG],Fin_TOT);
    Fout_h(:,i) = Fout(:,1);
    Fout_m(:,i) = Fout(:,2);
    Fout_t(:,i) = Fout(:,3);
    Fout_b(:,i) = Fout(:,4);
    Fout_d(:,i) = Fout(:,5);
    Selectivity(:,i) = Fout_b(:,i)./(Fin_t*ones(size(Fout_t(:,i)))-Fout_t(:,i));
    Conversion(:,i) = (Fin_t.*ones(size(Fout_t(:,i)))-Fout_t(:,i))./Fin_t.*ones(size(Fout_t(:,i)));
end


%% Couple Kinetics and Material Balances - 11 March 2024
% DoFs
SFvalues = 0.02:0.02:1;                                                     %[-]
F1h = zeros(length(Tvalues),length(SFvalues));
F1m = zeros(length(Tvalues),length(SFvalues));
F2t = zeros(length(Tvalues),length(SFvalues));
INh = zeros(length(Tvalues),length(SFvalues));
INm = zeros(length(Tvalues),length(SFvalues));
INt = zeros(length(Tvalues),length(SFvalues));
Tt = zeros(length(Tvalues),length(SFvalues));
Dd = zeros(length(Tvalues),length(SFvalues));
Bb = zeros(length(Tvalues),length(SFvalues));
RVh = zeros(length(Tvalues),length(SFvalues));
RVm = zeros(length(Tvalues),length(SFvalues));
Vh = zeros(length(Tvalues),length(SFvalues));
Vm = zeros(length(Tvalues),length(SFvalues));
Rh = zeros(length(Tvalues),length(SFvalues));
Rm = zeros(length(Tvalues),length(SFvalues));
V = zeros(length(Tvalues),length(SFvalues));

for i = 1:length(Tvalues)
    T = Tvalues(i);
    % Guess made by the SF = 100% and the basic solution of the MBs in the
    % same order the solutions are made
    guess = [2000 100 200 2000 100 400 200 6 265 1000 200 1000 200 0 0 50];

    for j = length(SFvalues):(-1):1
        SF = SFvalues(j);
        solution = fsolve(@(unk)MaterialBalances(unk,SF,T), guess);
        F1h(i,j) = solution(1);    F1m(i,j) = solution(2);    F2t(i,j) = solution(3);    INh(i,j) = solution(4);    INm(i,j) = solution(5);
        INt(i,j) = solution(6);    Tt(i,j) = solution(7);     Dd(i,j) = solution(8);     Bb(i,j) = solution(9);     RVh(i,j) = solution(10);
        RVm(i,j) = solution(11);   Vh(i,j) = solution(12);    Vm(i,j) = solution(13);    Rh(i,j) = solution(14);    Rm(i,j) = solution(15);
        V(i,j) = solution(16);
        guess = solution;    
    end
end
%% Evaluation of Adiabatic Temperature - 18 March 2024
Nin= [INh(:,end) INm(:,end) INt(:,end) zeros(length(Tvalues),1) zeros(length(Tvalues),1)].*3600; %[kmol/s] 
Nout = [RVh(:,end) RVm(:,end) Tt(:,end) Bb(:,end) 0*Dd(:,end)];
for i=1:length(Tvalues)
    no=Nin(i,:); %line vector    
    Tin=Tvalues(i);
    x0 = Tin;
    Toutsolution=fsolve(@(Tout)DeltaT_ad(Tout,no,Tin),x0);
    DeltaTadiabitc(i,:)=Toutsolution-Tin;
end

%% EP2 Calculation - 18 March 2024
Keys = {'Benzene','Toluene','Hydrogen','Byphenil', 'Methane'};
Cost_sell = [12.5, 8.8, 2.1, 7.4, 2.1];                                     % [€/kmol]
Cost_sell_Materials = containers.Map(Keys,Cost_sell);
delH_burn = [1.41, 1.68, 0.123, 2.688, 0.383];                              % [MBtu/kmol]
Cost_burn_Materials = containers.Map(Keys, delH_burn*4);                    % [€/MBtu]

% 1st Scenario [M€/y] - SELL BIPHENYL
EP2_sell = 1e-6*8000*(Cost_sell_Materials('Benzene')*Bb + Cost_sell_Materials('Byphenil')*Dd...
    - Cost_sell_Materials('Hydrogen')*F1h - Cost_sell_Materials('Methane')*F1m...
    - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogen')*Vh...
    + Cost_burn_Materials('Methane')*Vm);

% 2nd Scenario [M€/y] - BURN BIPHENYL
EP2_burn = 1e-6*8000*(Cost_sell_Materials('Benzene')*Bb + Cost_burn_Materials('Byphenil')*Dd...
    - Cost_sell_Materials('Hydrogen')*F1h - Cost_sell_Materials('Methane')*F1m...
    - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogen')*Vh...
    + Cost_burn_Materials('Methane')*Vm);
%% Excercise 4 - Calculating EP3 - 25 March 2024

% Reactor
M_and_S = 1100;
Dep_time=5;                                                                 %[y] Depreciation time.
H_D_ratio = 10;     % Can be varied from 6 to 10
D = ones(4,50);
for i=1:4
    for j=1:50
        D(i,j) = (4*V(i,j)/(pi*H_D_ratio))^(1/3);
    end
end
H = H_D_ratio*D;
Fp_reactor_out = 1.45;
Material = {'Carbon Steel', 'Stainless Steel 316','Monel', 'Titanium'};
Fm_val = [1, 3.67,6.34, 7.89];
Fm = containers.Map(Material,Fm_val);
Fm_desired = double(Fm('Stainless Steel 316'));
Fc = Fm_desired*Fp_reactor_out;
IC_reac = 1.15*M_and_S/280*101.9*(D.^1.066).*(H.^0.802)*(2.18+Fc)/1e6;      %[M€] I.C. of the reactor at various T and split factor.
IC_reac_dep=IC_reac/Dep_time;                                               %[M€/y]


% Compressor

Pin_IMP=449;                                                                %[psi]
Pout_IMP=537;                                                               %[psi]
R=8.314;                                                                    %[J/mol/K]
T1=273.15+35;                                                               %[K]
Eff=0.9;                                                                    % Mechanical efficiency.
Eff_Elect=0.9;                                                              % Electric efficiency.
Fc_Comp=1;
Elec_cost=0.061095;                                                         %[€/kWh]

% CAPEX evaluation (compressor)
Beta_IMP=Pout_IMP/Pin_IMP;
yRh=Rh./(Rm+Rh);
yRm=Rm./(Rm+Rh);

gamma_h=yRh.*0.29;                                                          % [-]
gamma_m=yRm.*0.23;                                                          % [-]
gamma_mix=gamma_h+gamma_m;                                                  % [-]

Duty_Id=R*T1*(Beta_IMP.^gamma_mix-1)./gamma_mix;                            %[J/mol] % Ideal specific duty.
Duty_Comp=Duty_Id/Eff;                                                      %[J/mol]
Duty_Elect=Duty_Comp/Eff_Elect;                                             %[J/mol]

Rmix=(Rh+Rm)*1000/3600;                                                     %[mol/s] Total flowrate to be compressed.

Power_Comp=Rmix.*Duty_Elect;                                                %[W] Actual electric power required for the compressor.
Power_Comp_IMP=Power_Comp/1000*1.341;                                       %[bph] Converting from W to bhp.

IC_Comp=M_and_S/280*517.5*Power_Comp_IMP.^0.82*(2.11+Fc_Comp)/1e6;          %[M€] CAPEX of the compressor.
CAPEX_Comp_Dep=IC_Comp/Dep_time;                                            %[M€/y] Compressor CAPEX with depreciation time accounted for.


% OPEX evaluation (compressor)

Duty_Comp_Total=Power_Comp*8000/1000;                                       %[kWh/y] Duty of the compressor in a year.
OPEX_Comp=Elec_cost*Duty_Comp_Total/1e6;                                    %[M€/y] OPEX of the compressor.


% EP3 evaluation
EP3_burn = EP2_burn - IC_reac_dep - CAPEX_Comp_Dep - OPEX_Comp;

%% Data collection of the best case for each temperature

index_600=find(EP3_burn==max(EP3_burn(1,:)));
index_650=find(EP3_burn==max(EP3_burn(2,:)));
index_700=find(EP3_burn==max(EP3_burn(3,:)));
index_750=find(EP3_burn==max(EP3_burn(4,:)));

SF_best=[SFvalues(index_600) SFvalues(index_650) SFvalues(index_700) SFvalues(index_750)]';
F1h_best=[F1h(1,index_600) F1h(2,index_650) F1h(3,index_700) F1h(4,index_750)]';
F1m_best=[F1m(1,index_600) F1m(2,index_650) F1m(3,index_700) F1m(4,index_750)]';
F2t_best=[F2t(1,index_600) F2t(2,index_650) F2t(3,index_700) F2t(4,index_750)]';
Rh_best=[Rh(1,index_600) Rh(2,index_650) Rh(3,index_700) Rh(4,index_750)]';
Rm_best=[Rm(1,index_600) Rm(2,index_650) Rm(3,index_700) Rm(4,index_750)]';
Tt_best=[Tt(1,index_600) Tt(2,index_650) Tt(3,index_700) Tt(4,index_750)]';
V_best=[V(1,index_600) V(2,index_650) V(3,index_700) V(4,index_750)]';
D_best=[D(1,index_600) D(2,index_650) D(3,index_700) D(4,index_750)]';
H_best=[H(1,index_600) H(2,index_650) H(3,index_700) H(4,index_750)]';
Temperatures = Tvalues'-273.15;

F1_best=F1h_best+F1m_best;
Rtot_best=Rm_best+Rh_best;
xmR=Rm_best./Rtot_best;
xhR=Rh_best./Rtot_best;
clc;
disp(table(Temperatures, SF_best,F1h_best,F1m_best,F2t_best,Rh_best,Rm_best,Tt_best,V_best,D_best,H_best))
disp(table(F1_best,Rtot_best,xhR,xmR))


EP2_best=[EP2_burn(index_600) EP2_burn(index_650) EP2_burn(index_700) EP2_burn(index_750)]; %Maximum values of EP2 at various T.
EP3_best=[EP3_burn(index_600) EP3_burn(index_650) EP3_burn(index_700) EP3_burn(index_750)]; %Maximum values of EP3 at various T.

%% Calculating EP4

% COLUMNS
T_S = [24/12 18/12 12/12];                                                  % [ft] tray spacing.
Fs = [1 1.4 2.2];                                                           % Fs at various tray spacing (24, 18 and 12 in).
H_Top_Bottom = 5*100/2.54/12;                                               % [ft]
Ft = 0;                                                                     % Ft for sieve trays.
Fm_CS = 0;                                                                  % Fm for Carbon Steel.
Fm_SS = 1.7;                                                                % Fm for Stainless Steel. ONLY FOR THE STABILIZER.


% Stabilizer
N_stabilizer = [3 3 3 3];                                                   % Number of trays for the stabilizer at various T.
D_stabilizer = [3 3 3 3.5];                                                 % [ft] diameter of the stabilizer at various T.
HETP = [1.335 1.334 1.333 1.332];                                           % [ft] HETP of the stabilizer at various T. For the price I can use 18 in since it's the closest.
Fp_vessel_stab = 1;
Fm_vessel_stab = 3.67;
Fc_vessel_stab = Fp_vessel_stab*Fm_vessel_stab;

H_tot_stabilizer = (N_stabilizer-1)*T_S(2)+H_Top_Bottom;                    % [ft] Height at various T.
Fc_stabilizer = Fs(2)+Ft+Fm_SS;
IC_stabilizer_trays = M_and_S/280*4.7.*D_stabilizer.^1.55.*H_tot_stabilizer.*Fc_stabilizer/1e6; % [M€] IC of the stabilizer trays at various T.
IC_stabilizer_vessel = M_and_S/280*101.9*(D_stabilizer.^1.066).*(H_tot_stabilizer.^0.802)*(2.18+Fc_vessel_stab)/1e6; % [M€] IC of the stabilizer "vessel" at various T.
IC_stabilizer = IC_stabilizer_trays + IC_stabilizer_vessel;                 % [M€] Total IC of the stabilizer at various T.

% Product column
N_product = [51 52 52 53];                                                  % Number of trays for the product column at various T.
D_product_24 = [7.5 7.5 8 8];                                               % Diameter of the product column at various T IF the tray spacing is 24 in.
D_product_18 = [8 8 8.5 8.5];                                               % Diameter of the product column at various T IF the tray spacing is 18 in.
D_product_12 = [8.5 9.5 9.5 9.5];                                           % Diameter of the product column at various T IF the tray spacing is 12 in.
D_product = [D_product_24; D_product_18; D_product_12];
Fp_vessel_prod = 1;
Fm_vessel_prod = 1;
Fc_vessel_prod = Fp_vessel_prod*Fm_vessel_prod;

for i=1:length(T_S)
    H_tot_product(i,:) = (N_product-1).*T_S(i)+H_Top_Bottom;                % [ft] Height at various T.
    Fc_product(i) = Fs(i)+Ft+Fm_CS;
    IC_product_trays(i,:) = M_and_S/280*4.7.*D_product(i,:).^1.55.*H_tot_product(i,:).*Fc_product(i)/1e6; % [M€] IC of the product column trays at various T.
    IC_product_vessel(i,:) = M_and_S/280*101.9*(D_product(i,:).^1.066).*(H_tot_product(i,:).^0.802)*(2.18+Fc_vessel_prod)/1e6; % [M€] IC of the product column "vessel" at various T.
end

IC_product = IC_product_trays + IC_product_vessel;                          % [M€] Total IC of the product column at various T.

for i=1:length(Tvalues)
    IC_best_tray_product(i) = min(IC_product(:,i));    
    [row,col] = find(IC_product(:,i)==IC_best_tray_product(i));              % Finding the position of the lowest I.C for each T.
    Index_tray_sizing_prod(i) = row;                                         % Finding the tray sizing corresponding to the lowest I.C.
    Best_tray_prod(i) = T_S(Index_tray_sizing_prod(i));                      % [ft] Vector with the best tray spacing for the product column at various T.
end

% Recycle column
N_recycle = [75 74 76 80];                                                  % Number of trays for the recycle column at various T.
D_recycle_24 = [3.5 4 3.5 4.5];                                             % Diameter of the recycle column at various T IF the tray spacing is 24 in.
D_recycle_18 = [3.5 4 4 4.5];                                               % Diameter of the recycle column at various T IF the tray spacing is 18 in.
D_recycle_12 = [3.5 4 4.5 4.5];                                             % Diameter of the recycle column at various T IF the tray spacing is 12 in.
D_recycle = [D_recycle_24; D_recycle_18; D_recycle_12];
Fp_vessel_rec = 1;
Fm_vessel_rec = 1;
Fc_vessel_rec = Fp_vessel_rec*Fm_vessel_rec;

for i=1:length(T_S)
    H_tot_recycle(i,:) = (N_recycle-1).*T_S(i)+H_Top_Bottom;                % [ft] Height at various T.
    Fc_recycle(i) = Fs(i)+Ft+Fm_CS;
    IC_recycle_trays(i,:) = M_and_S/280*4.7.*D_recycle(i,:).^1.55.*H_tot_recycle(i,:).*Fc_recycle(i)/1e6; % [M€] IC of the recycle column trays at various T.
    IC_recycle_vessel(i,:) = M_and_S/280*101.9*(D_recycle(i,:).^1.066).*(H_tot_recycle(i,:).^0.802)*(2.18+Fc_vessel_rec)/1e6; % [M€] IC of the recycle column "vessel" at various T.
end

IC_recycle = IC_recycle_trays + IC_recycle_vessel;                          % [M€] Total IC of the recycle column at various T.

for i=1:length(Tvalues)
    IC_best_tray_recycle(i) = min(IC_recycle(:,i));    
    [row,col] = find(IC_recycle(:,i)==IC_best_tray_recycle(i));             % Finding the position of the lowest I.C for each T.
    Index_tray_sizing_rec(i) = row;                                         % Finding the tray sizing corresponding to the lowest I.C.
    Best_tray_rec(i) = T_S(Index_tray_sizing_rec(i));                       % [ft] Vector with the best tray spacing for the recycle column at various T.
end

IC_all_columns = IC_stabilizer+IC_best_tray_product+IC_best_tray_recycle;   % [M€] IC of all the columns in the process at various T.

% HEAT EXCHANGERS
Fm_CS_CS = 1;
Fm_CS_SS = 2.81;

Fd_fixed = 0.8;
Fd_kettle = 1.35;

Fp = 0;

% Condensers:

% Stabilizer:
A_cond_stab = [84.5527 83.7097 98.2651 96.5164];                            % [ft2] Heat exchange area for the condenser of the stabilizer at various T.
Fc_cond_stab = (Fd_fixed+Fp)*Fm_CS_SS;
IC_cond_stab = M_and_S/280*101.3.*A_cond_stab.^0.65.*(2.29+Fc_cond_stab)/1e6; % [M€] IC of the condenser of the stabilizer.

H2O_flow_stab = [0.3904 0.3895 0.4592 0.4545];                              % [USgal/s]
OPEX_cond_stab = 0.06/1000*H2O_flow_stab*3600*8000/1e6;                     % [M€/y] OC of the condenser of the stabilizer.

% Product column:
A_cond_prod = [2791.1538 2786.6763 2848.8197 2907.4507];                    % [ft2] Heat exchange area for the condenser of the product column at various T.
Fc_cond_prod = (Fd_fixed+Fp)*Fm_CS_CS;
IC_cond_prod = M_and_S/280*101.3.*A_cond_prod.^0.65.*(2.29+Fc_cond_prod)/1e6; % [M€] IC of the condenser of the product column.

H2O_flow_prod = [18.3961 18.3763 18.7960 19.1879];                          % [USgal/s]
OPEX_cond_prod = 0.06/1000*H2O_flow_prod*3600*8000/1e6;                     % [M€/y] OC of the condenser of the product column.

% Recycle column:
A_cond_rec = [225.5487 290.5370 333.4838 403.1514];                         % [ft2] Heat exchange area for the condenser of the recycle column at various T.
Fc_cond_rec = (Fd_fixed+Fp)*Fm_CS_CS;
IC_cond_rec = M_and_S/280*101.3.*A_cond_rec.^0.65.*(2.29+Fc_cond_rec)/1e6;  % [M€] IC of the condenser of the recycle column.

H2O_flow_rec = [2.6429 3.4094 3.9191 4.7447];                               % [USgal/s]
OPEX_cond_rec = 0.06/1000*H2O_flow_rec*3600*8000/1e6;                       % [M€/y] OC of the condenser of the recycle column.


A_cond = [A_cond_stab; A_cond_prod; A_cond_rec];                            % [ft2] heat exchange area of all the condensers at various T.
IC_cond_tot = IC_cond_stab+IC_cond_prod+IC_cond_rec;                        % [M€] IC of all the condensers at various T.
OPEX_cond_tot = OPEX_cond_stab+OPEX_cond_prod+OPEX_cond_rec;                % [M€/y] OC of all the condensers at various T.


% Reboilers

% Stabilizer:
A_reb_stab = [223.4027 248.6523 267.9519 296.3425];                         % [ft2] Heat exchange area for the reboiler of the stabilizer at various T.
Fc_reb_stab = (Fd_kettle+Fp)*Fm_CS_CS;
IC_reb_stab = M_and_S/280*101.3.*A_reb_stab.^0.65.*(2.29+Fc_reb_stab)/1e6;  % [M€] IC of the reboiler of the stabilizer.

Steam_reb_stab = [0.9043 1.0066 1.0847 1.1996];                             % [lb/s] Steam flowrate for the reboiler of the stabilizer at various T.
OPEX_reb_stab = 1.65/1000*Steam_reb_stab*3600*8000/1e6;                     % [M€/y] OC of the reboiler of the stabilizer.

% Product column:
A_reb_prod = [1514.2432 1771.5513 1833.8795 1881.3790];                     % [ft2] Heat exchange area for the reboiler of the product column at various T.
Fc_reb_prod = (Fd_kettle+Fp)*Fm_CS_CS;
IC_reb_prod = M_and_S/280*101.3.*A_reb_prod.^0.65.*(2.29+Fc_reb_prod)/1e6;  % [M€] IC of the reboiler of the product column.

Steam_reb_prod = [6.1297 7.1713 7.4236 7.6159];                             % [lb/s] Steam flowrate for the reboiler of the product column at various T.
OPEX_reb_prod = 1.65/1000*Steam_reb_prod*3600*8000/1e6;                     % [M€/y] OC of the reboiler of the product column.

% Recycle column
A_reb_rec = [273.4940 345.4657 393.4299 471.3628];                          % [ft2] Heat exchange area for the reboiler of the recycle column at various T.
Fc_reb_rec = (Fd_kettle+Fp)*Fm_CS_CS;
IC_reb_rec = M_and_S/280*101.3.*A_reb_rec.^0.65.*(2.29+Fc_reb_rec)/1e6;     % [M€] IC of the reboiler of the product column.

Steam_reb_rec = [1.1071 1.3985 1.5926 1.9081];                              % [lb/s] Steam flowrate for the reboiler of the recycle column at various T.
OPEX_reb_rec = 1.65/1000*Steam_reb_rec*3600*8000/1e6;                       % [M€/y] OC of the reboiler of the recycle column.


A_reb = [A_reb_stab; A_reb_prod; A_reb_rec];                                % [ft2] heat exchange area of all the reboilers at various T.
IC_reb_tot = IC_reb_stab+IC_reb_prod+IC_reb_rec;                            % [M€] IC of all the reboiler at various T.
OPEX_reb_tot = OPEX_reb_stab+OPEX_reb_prod+OPEX_reb_rec;                    % [M€/y] OC of all the reboilers at various T.

Separation_cost = (IC_all_columns+IC_cond_tot + IC_reb_tot)/Dep_time + ...
OPEX_cond_tot+OPEX_reb_tot;                                                 % [M€/y] CAPEX+OPEX of the separation section.

EP4_best = EP3_best-Separation_cost;                                        % [M€/y] Evaluation of the EP4 at various T.


%% GRAPHICS

figure(1)
% Amount of hydrogen flowing along the reactor
subplot(2,2,1)
plot(Vout,Fout_h)
title('(a) Hydrogen molar flow rate vs. Reactor volume')
xlabel('Volume [m^3]')
grid on
ylabel('Flow rate [kmol/h]')
legend('600°C','650°C','700°C','750°C')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(2000,1950,str,'FontSize',8)
subplot(2,2,2)
% Amount of methane flowing along the reactor
plot(Vout,Fout_m)
grid on
title('(b) Methane molar flow rate vs. Reactor volume')
xlabel('Volume [m^3]')
ylabel('Flow rate [kmol/h]')
legend('600°C','650°C','700°C','750°C','location','southeast')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(3000,250,str,'FontSize',8)
subplot(2,2,3)
% Amount of toluene flowing along the reactor
plot(Vout,Fout_t)
grid on
title('(c) Toluene molar flow rate vs. Reactor volume')
xlabel('Volume [m^3]')
ylabel('Flow rate [kmol/h]')
legend('600°C','650°C','700°C','750°C')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(4000,150,str,'FontSize',8)
subplot(2,2,4)
% Amount of benzene flowing along the reactor
plot(Vout,Fout_b)
grid on
title('(d) Benzene molar flow rate vs. Reactor volume')
xlabel('Volume [m^3]')
ylabel('Flow rate [kmol/h]')
legend('600°C','650°C','700°C','750°C')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(5000,150,str,'FontSize',8)

% Amount of biphenyl flowing along the reactor
figure(2)
plot(Vout,Fout_d)
grid on
title('Biphenyl molar flow rate vs. Reactor volume')
xlabel('Volume [m^3]')
ylabel('Flow rate [kmol/h]')
legend('600°C','650°C','700°C','750°C','location','southeast')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(3500,20,str,'FontSize',8)

% Conversion versus volume
figure(3)
plot(Vout,Conversion*100)
title('Conversion of toluene vs. Reactor volume')
xlabel('Volume [m^3]'); ylabel('Conversion (%)'); grid on
legend('600°C','650°C','700°C','750°C','location','best')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(2000,80,str,'FontSize',8)

figure(4)
% Selectivity versus volume
plot(Vout,Selectivity*100)
title('Selectivity of the dealkylation vs. Reactor volume')
xlabel('Volume [m^3]'); ylabel('Selectivity (%)'); grid on
legend('600°C','650°C','700°C','750°C','location','best')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(2000,96,str,'FontSize',8)


% Selectivity of the dealkylation vs Conversion of toluene
figure(5)
plot(Conversion*100,Selectivity*100)
title('Selectivity of the dealkylation vs. Conversion of toluene')
grid on
xlabel('Conversion (%)')
ylabel('Selectivity (%)')
yline(96)
legend('600°C','650°C','700°C','750°C','location','southwest')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(20,40,str,'FontSize',8)
str096 = {'96'};
text(-5.5,96,str096)


% Conversion vs. Temperature, with a selectivity of 96%
Conv_t = 100*(1 - Tt./INt);
figure(6)
plot(Temperatures,Conv_t(:,end),'-.o')
title('Conversion (%) vs. Temperature, with a selectivity of 96%')
grid on
xlabel('Temperature [°C]')
ylabel('Conversion (%)')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.2,0.4,str,'FontSize',8)


% Reactor volume vs. Temperature, with a selectivity of 96%
figure(7)
plot(Temperatures,V(:,end),'-.o')
title('Reactor volume vs. Temperature, with a selectivity of 96%')
grid on
xlabel('Temperature [°C]')
ylabel('Reactor volume [m3]')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.2,0.4,str,'FontSize',8)

% Adiabatic Temperature Difference vs. Temperature
figure(8)
plot(Temperatures,DeltaTadiabitc,'-.o')
title('Adiabatic Temperature Difference vs. Temperature')
grid on; xlabel('Temperature [°C]'); ylabel('Adiabatic Temperature Difference [°C]')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.2,0.4,str,'FontSize',8)

% Amount of hydrogen in the vent stream vs. Split Factor
Vent_frac = Vh./(Vh+Vm);                                                    % [///] Hydrogen molar fraction in the vent stream.
figure(9)
plot(SFvalues,Vent_frac(1,:))
hold on
for i=2:4
    plot(SFvalues,Vent_frac(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Hydrogen (molar fraction)')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,0.3,str,'FontSize',8)
legend('600°C','650°C','700°C','750°C','location','Northwest')
title('Amount of hydrogen in the vent stream vs. Split Factor')

% EP2 (both "sell" and "burn" scenarios) vs. Amount of hydrogen in the vent
figure(10)
plot(Vent_frac(1,:),EP2_sell(1,:),':',color=[0 0.4470 0.7410])
hold on
plot(Vent_frac(2,:),EP2_sell(2,:),':',color=[0.8500 0.3250 0.0980])
plot(Vent_frac(3,:),EP2_sell(3,:),':',color=[0.9290 0.6940 0.1250])
plot(Vent_frac(4,:),EP2_sell(4,:),':',color=[0.4940 0.1840 0.5560])

plot(Vent_frac(1,:),EP2_burn(1,:),color=[0 0.4470 0.7410])
plot(Vent_frac(2,:),EP2_burn(2,:),color=[0.8500 0.3250 0.0980])
plot(Vent_frac(3,:),EP2_burn(3,:),color=[0.9290 0.6940 0.1250])
plot(Vent_frac(4,:),EP2_burn(4,:),color=[0.4940 0.1840 0.5560])
grid on
ylim([0 +Inf])
xlabel('H2 (molar fraction)'); ylabel('EP2 (M€/y)')
legend('600°C (Sell)','650°C (Sell)','700°C (Sell)','750°C (Sell)',...
    '600°C (Burn)','650°C (Burn)','700°C (Burn)','750°C (Burn)',...
    'location','northwest','NumColumns',2,'location','northeast')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.05,1,str,'FontSize',8)
title('EP2 vs. Amount of hydrogen in the vent stream')
hold off

% EP2 (both "sell" and "burn" scenarios) vs. Split factor
figure(11)
plot(SFvalues,EP2_sell(1,:),':',color=[0 0.4470 0.7410])
hold on
plot(SFvalues,EP2_sell(2,:),':',color=[0.8500 0.3250 0.0980])
plot(SFvalues,EP2_sell(3,:),':',color=[0.9290 0.6940 0.1250])
plot(SFvalues,EP2_sell(4,:),':',color=[0.4940 0.1840 0.5560])

plot(SFvalues,EP2_burn(1,:),color=[0 0.4470 0.7410])
plot(SFvalues,EP2_burn(2,:),color=[0.8500 0.3250 0.0980])
plot(SFvalues,EP2_burn(3,:),color=[0.9290 0.6940 0.1250])
plot(SFvalues,EP2_burn(4,:),color=[0.4940 0.1840 0.5560])
grid on
ylim([0 +Inf])
xlabel('Split Factor'); ylabel('EP2 (M€/y)')
legend('600°C (Sell)','650°C (Sell)','700°C (Sell)','750°C (Sell)',...
    '600°C (Burn)','650°C (Burn)','700°C (Burn)','750°C (Burn)',...
    'location','northwest','NumColumns',2,'location','northeast')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.05,1,str,'FontSize',8)
title('EP2 vs. Split factor')
hold off

% EP2 (both "sell" and "burn" scenarios) vs. Temperature
figure(12)
Graph1 = plot(Tvalues-273.15, EP2_sell(:,1),'-.',color=[1 0.1+1/80 0]);
hold on
for i=5:5:50
    plot(Tvalues-273.15,EP2_sell(:,i),'-.',color=[1 0.1+i/80 0])
end
grid on

Graph2 = plot(Tvalues-273.15, EP2_burn(:,1),color=[1 0.1+1/80 0]);
hold on
for i=5:5:50
    plot(Tvalues-273.15,EP2_burn(:,i),color=[1 0.1+i/80 0])
end
grid on
yline(0)
xlabel('Temperature [°C]'); ylabel('EP2 (M€/y)')
title('EP2 vs. Temperature');
str1 = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(610,8,str1,'FontSize',8)
str2 = {'0.02', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1',};
x_str2 = [751 751 751 751 751 751 751 751 751 751 751];
y_str2 = [6 4 2.2 -0.4 -2.7 -5 -7.4 -9.8 -12.1 -14.5 -17];
text(x_str2,y_str2,str2,'FontSize',5,'color','blue')
str3 = {'SPLIT', 'FACTOR', 'VALUES:'};
x_str2 = [751 751 751];
y_str3 = [10 9 8];
text(x_str2,y_str3,str3,'FontSize',6,'color','blue')
Legendset = [Graph1 Graph2];
lgd=legend(Legendset,'Sell','Burn','location','Southwest');
hold off

% EP2 (both "sell" and "burn" scenarios) vs. Conversion of toluene
figure(13)
Graph1 = plot(Conv_t(:,1),EP2_sell(:,1),'-.',color=[0 0 0]);
hold on
for i=5:5:50
    plot(Conv_t(:,i),EP2_sell(:,i),'-.',color=[0.3*i/50 0.4+i/100 0.4+i/300])
end

Graph2 = plot(Conv_t(:,1),EP2_burn(:,1),color=[0 0 0]);
hold on
for i=5:5:50
    plot(Conv_t(:,i),EP2_burn(:,i),color=[0.3*i/50 0.4+i/100 0.4+i/300])
end
grid on
yline(0)
xlabel('Conversion (%)'); ylabel('EP2 (M€/y)')
title('EP2 vs. Conversion of toluene')
str1 = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(67,8,str1,'FontSize',8)
str2 = {'0.02', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1',};
x_str2 = [91 81 78 77 76 75.5 75.1 74.9 74.8 74.7 74.7];
y_str2 = [6.5 5 3 0.8 -1.2 -3.3 -5.5 -7.5 -9.5 -11.6 -14];
text(x_str2,y_str2,str2,'FontSize',5,'color','blue')
Legendset = [Graph1 Graph2];
lgd=legend(Legendset,'Sell','Burn','location','Southeast');
hold off

% Reactor volume vs. Split factor
figure(14)
plot(SFvalues,V(1,:))
hold on
for i=2:4
    plot(SFvalues,V(i,:))
end
hold off
grid on
xlabel('Split Factor'); ylabel('Reactor Volume [m^3]')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,1e4,str,'FontSize',8)
legend('600°C','650°C','700°C','750°C')
title('Reactor volume vs. Split factor')

% Reactor diameter vs. Split factor
figure(15)
plot(SFvalues,D(1,:))
hold on 
for i=2:4
    plot(SFvalues,D(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Reactor Diameter [m]')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,10,str,'FontSize',8)
legend('600°C','650°C','700°C','750°C')
title('Reactor diameter vs. Split factor')
hold off

% I.C. of the reactor vs. Split factor
figure(16)
plot(SFvalues,IC_reac(1,:))
hold on
for i=2:4
    plot(SFvalues,IC_reac(i,:))
end
grid on
xlabel('Split Factor'); ylabel('IC reactor [M€]')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,1.5,str,'FontSize',8)
legend('600°C','650°C','700°C','750°C')
title('Installation cost of the reactor vs. Split factor')
hold off

% Recycle stream vs. Split factor
Rtot = Rh + Rm;                                                             % [kmol/h]
figure(17)
plot(SFvalues,Rtot(1,:))
hold on
for i=2:4
    plot(SFvalues,Rtot(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Recycle stream [kmol/h]')
legend('600°C','650°C','700°C','750°C')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,8000,str,'FontSize',8)
title('Recycle stream vs. Split factor')

% I.C. of the compressor vs. Split factor
figure(18)
plot(SFvalues,IC_Comp(1,:))
hold on
for i=2:4
    plot(SFvalues,IC_Comp(i,:))
end
hold off
title('Installation cost of the compressor vs. Split factor')
xlabel('Split factor'); ylabel('I.C. of the compressor [M€]')
grid on
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,2,str,'FontSize',8)
legend('600°C','650°C','700°C','750°C')

% O.C. of the compressor vs. Split factor
figure(19)
plot(SFvalues,OPEX_Comp(1,:))
hold on
for i=2:4
    plot(SFvalues,OPEX_Comp(i,:))
end
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.5,0.6,str,'FontSize',8)
title('Operating cost of the compressor vs. Split factor')
xlabel('Split factor'); ylabel('O.C. of the compressor [M€/y]')
grid on
legend('600°C','650°C','700°C','750°C')
hold off

% EP3 ("burn" scenario) vs. Split factor
figure(20)
plot(SFvalues,EP3_burn(1,:),color=[0 0.4470 0.7410])
hold on
plot(SFvalues,EP3_burn(2,:),color=[0.8500 0.3250 0.0980])
plot(SFvalues,EP3_burn(3,:),color=[0.9290 0.6940 0.1250])
plot(SFvalues,EP3_burn(4,:),color=[0.4940 0.1840 0.5560])
xlabel('Split Factor'); ylabel('EP3 [M€/y]')
grid on
ylim([0 5])
title('EP3 vs. Split factor')
legend('600°C','650°C','700°C','750°C','location','northeast')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(0.05,0.5,str,'FontSize',8)
hold off

% Separation section total cost vs. Temperature
figure(21)
plot(Tvalues-273.15,Separation_cost,'-o')
title('Separation section total cost vs. Temperature')
xlabel('Temperature [°C]'); ylabel('Separation section cost [M€/y]')
grid on
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(650,0.56,str,'FontSize',8)

% EP3 vs. Temperature
figure(22)
plot(Tvalues-273.15,EP3_best,'-o')
title('EP3 vs. Temperature')
xlabel('Temperature [°C]'); ylabel('EP3 [M€/y]')
grid on
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(650,4.22,str,'FontSize',8)

% EP4 vs. Temperature
figure(23)
plot(Tvalues-273.15,EP4_best,'-o')
title('EP4 vs. Temperature')
xlabel('Temperature [°C]'); ylabel('EP4 [M€/y]')
grid on
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(650,3.61,str,'FontSize',8)

% Economic Potential vs. Temperature
figure(24)
plot(Tvalues-273.15,EP4_best,'-o')
title('Economic Potential vs. Temperature')
hold on
plot(Tvalues-273.15,EP2_best,'--x')
hold on
plot(Tvalues-273.15,EP3_best,'-.^')
xlabel('Temperature [°C]'); ylabel('Economic potential [M€/y]')
grid on
legend('EP4','EP2','EP3','location','northwest')
str = {'© Matteo Robbiano &','Mohammad Sina Ghanbari Pakdehi'};
text(625,5,str,'FontSize',8)
hold off
