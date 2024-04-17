clearvars, close all, clc
%% Couple Kinetics and Material Balances - 11 March 2024
% DoFs
Tvalues = (600:50:750)+273.15;              %[K]
SFvalues = 0.02:0.02:1;                     %[-]

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
R = Rh + Rm;            %[kmol/h]
figure(1)
plot(SFvalues,R(1,:))
hold on
for i=2:4
    plot(SFvalues,R(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Recycle Stream')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('Split Factor vs Recycle Stream')
%% Plotting H2 fraction in vent stream versus Split Factor - 18 March 2024
Vent_frac = Vh./(Vh+Vm);        % Hydrogen fraction in vent stream
figure(2)
plot(SFvalues,Vent_frac(1,:))
hold on
for i=2:4
    plot(SFvalues,Vent_frac(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Vh [kmol/hr]')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('Amount of Hydrogen in Vent Stream Versus Split Factor')
%% EP2 Calculation - 18 March 2024
Keys = {'Benzene','Toluene','Hydrogene','Byphenil', 'Methane'};
Cost_sell = [12.5, 8.8, 2.1, 7.4, 2.1];
Cost_sell_Materials = containers.Map(Keys,Cost_sell);
delH_burn = [1.41, 1.68, 0.123, 2.688, 0.383];
Cost_burn_Materials = containers.Map(Keys, delH_burn*4);

% 1st Scenario [M€/8000hr]
EP2_sell = 1e-6*8000*(Cost_sell_Materials('Benzene')*Bb + Cost_sell_Materials('Byphenil')*Dd - Cost_sell_Materials('Hydrogene')*(F1h+F1m) - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogene')*Vh + Cost_burn_Materials('Methane')*Vm);
% 2nd Scenario [M€/8000hr]
EP2_burn = 1e-6*8000*(Cost_sell_Materials('Benzene')*Bb + Cost_burn_Materials('Byphenil')*Dd - Cost_sell_Materials('Hydrogene')*(F1h+F1m) - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogene')*Vh + Cost_burn_Materials('Methane')*Vm);


% Plot EP2 (both "sell" and "burn") vs. Split factor, by imposing a selectivity of at least 96%
figure(3)
subplot(1,2,1);
plot(SFvalues,EP2_sell(1,:))
hold on
plot(SFvalues,EP2_sell(2,:))
plot(SFvalues,EP2_sell(3,:))
plot(SFvalues,EP2_sell(4,:))
grid on
xlabel('Split Factor'); ylabel('EP (M€/year)')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('EP2 for sell scenario')

subplot(1,2,2);
plot(SFvalues,EP2_burn(1,:))
hold on
plot(SFvalues,EP2_burn(2,:))
plot(SFvalues,EP2_burn(3,:))
plot(SFvalues,EP2_burn(4,:))
grid on
xlabel('Split Factor'); ylabel('EP (M€/8000hr)')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('EP2 for burn scenario')

% Plot  EP2 (both "sell" and "burn") vs. Temperature, as a function of the split factor
figure(4)
subplot(1,2,1);
plot(Tvalues, EP2_sell(:,1))
hold on
for i=2:40
    plot(Tvalues,EP2_sell(:,i))
end
grid on
xlabel('Temperature [K]'); ylabel('EP2 (M€/year)')
title('EP2 for sell scenario')

subplot(1,2,2);
plot(Tvalues, EP2_burn(:,1))
hold on
for i=2:40
    plot(Tvalues,EP2_burn(:,i))
end
grid on
xlabel('Temperature [K]'); ylabel('EP2 (M€/year)')
title('EP2 for burn scenario')

% Conversion versus selectivity graph
Conv_t = 100*(1 - Tt./INt);
figure(5)
subplot(1,2,1);
plot(Conv_t(1,:),EP2_sell(1,:))
hold on
for i=2:4
    plot(Conv_t(i,:),EP2_sell(i,:))
end
xlabel('Conversion (%)'); ylabel('EP2 (M€/year)')
grid on
title('EP2 for sell scenario')
subplot(1,2,2);
plot(Conv_t(1,:),EP2_burn(1,:))
hold on
for i=2:4
    plot(Conv_t(i,:),EP2_burn(i,:))
end
grid on
xlabel('Conversion (%)'); ylabel('EP2 (M€/year)')
title('EP2 for burn scenario')

%% Excercise 4 - Calculating EP3 - 25 March 2024
% This is for the Reactor
M_and_S = 1100;
H_D_ratio = 10; % Can be varied from 6 to 10
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
IC_reac = 1.15*M_and_S/280*101.9*(D.^1.066).*(H.^0.802)*(2.18+Fc)/1e6; %[M€]
IC_reac_dep=IC_reac/5; %[M€/y]

figure(6)
plot(SFvalues,V(1,:))
hold on
for i=2:4
    plot(SFvalues,V(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Reactor Volume [m^3]')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')

figure(7)
plot(SFvalues,D(1,:))
hold on 
for i=2:4
    plot(SFvalues,D(i,:))
end
grid on
xlabel('Split Factor'); ylabel('Reactor Diameter [m^3]')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')

figure(8)
plot(SFvalues,IC_reac(1,:))
hold on
for i=2:4
    plot(SFvalues,IC_reac(i,:))
end
grid on
xlabel('Split Factor'); ylabel('IC reactor [M€]')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')


% This is for the Compressor

M_S=1110;
Pin_IMP=449;            %[psi]
Pout_IMP=537;           %[psi]
R=8.314;                %[J/mol/K]
T1=273.15+35;           %[K]
Eff=0.9;                % Mechanical efficiency.
Eff_Elect=0.9;          % Electric efficiency.
Fc_Comp=1;
Elec_cost=0.061095;     %[€/kWh]

% CAPEX evaluation (compressor)
Beta_IMP=Pout_IMP/Pin_IMP;
yRh=Rh./(Rm+Rh);
yRm=Rm./(Rm+Rh);

gamma_h=yRh.*0.29;
gamma_m=yRm.*0.23;

gamma_mix=gamma_h+gamma_m;

Duty_Id=R*T1*(Beta_IMP.^gamma_mix-1)./gamma_mix; %[J/mol] % Ideal specific duty.

Duty_Comp=Duty_Id/Eff; %[J/mol]

Duty_Elect=Duty_Comp/Eff_Elect; %[J/mol]

Rmix=(Rh+Rm)*1000/3600; %[mol/s] Total flowrate to be compressed.

Power_Comp=Rmix.*Duty_Elect; %[W] Actual electric power required for the compressor.

Power_Comp_IMP=Power_Comp/1000*1.341; %[bph] Converting from W to bhp.

IC_Comp=M_S/280*517.5*Power_Comp_IMP.^0.82*(2.11+Fc_Comp)/1e6; %[M€] CAPEX of the compressor.

CAPEX_Comp_Dep=IC_Comp/5; %[M€/y] Compressor CAPEX with depreciation time accounted for.


% OPEX evaluation (compressor)

Duty_Comp_Total=Power_Comp*8000/1000; %[kWh/y] Duty of the compressor in a year.

OPEX_Comp=Elec_cost*Duty_Comp_Total/1e6; %[M€/y] OPEX of the compressor.

% Graphics
figure(9)
plot(SFvalues,IC_Comp(1,:))
hold on
for i=2:4
    plot(SFvalues,IC_Comp(i,:))
end
% plot(SFvalues,IC_Comp(2,:),'g')
% hold on
% plot(SFvalues,IC_Comp(3,:),'b')
% hold on
% plot(SFvalues,IC_Comp(4,:),'k')
hold off
title('I.C. of the compressor vs. Split factor')
xlabel('Split factor'); ylabel('I.C. of the compressor [M€/y]')
grid on
legend('600°C','650°C','700°C','750°C')


figure(10)
plot(SFvalues,OPEX_Comp(1,:))
hold on
for i=2:4
    plot(SFvalues,OPEX_Comp(i,:))
end
% plot(SFvalues,OPEX_Comp(2,:),'g')
% hold on
% plot(SFvalues,OPEX_Comp(3,:),'b')
% hold on
% plot(SFvalues,OPEX_Comp(4,:),'k')
% hold off
title('O.C. of the compressor vs. Split factor')
xlabel('Split factor'); ylabel('O.C. of the compressor [M€/y]')
grid on
legend('600°C','650°C','700°C','750°C')

% Calculating EP3 in burn scenario
EP3_burn = EP2_burn - IC_reac_dep - CAPEX_Comp_Dep - OPEX_Comp;
% Calculating EP3 in sell scenario
EP3_sell = EP2_sell - IC_reac_dep - CAPEX_Comp_Dep - OPEX_Comp;

figure(11)
subplot(1,2,1)
plot(SFvalues,EP3_sell(1,:))
hold on
for i=2:4
    plot(SFvalues,EP3_sell(i,:))
end
xlabel('Split Factor'); ylabel('EP3 of Process[M€/y]')
grid on
title('Sell Scenario')

subplot(1,2,2)
plot(SFvalues,EP3_burn(1,:))
hold on
for i=2:4
    plot(SFvalues,EP3_burn(i,:))
end
xlabel('Split Factor'); ylabel('EP3 of Process[M€/y]')
grid on
title('Burn Scenario')

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
Temperatures = Tvalues';
clc;
disp(table(Temperatures, SF_best,F1h_best,F1m_best,F2t_best,Rh_best,Rm_best,Tt_best,V_best,D_best,H_best))
