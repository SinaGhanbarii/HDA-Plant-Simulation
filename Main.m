clearvars
close all
clc

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
plot(SFvalues,R)
xlabel('Split Factor [-]')
ylabel('Recycle, R [kmol/h]')

Vent_frac = Vh./(Vh+Vm);        % Hydrogen fraction in vent stream

figure(2)
plot(SFvalues,Vent_frac(1,:))
hold on
plot(SFvalues,Vent_frac(2,:))
plot(SFvalues,Vent_frac(3,:))
plot(SFvalues,Vent_frac(4,:))
grid on
xlabel('Split Factor'); ylabel('Vh [kmol/hr]')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('Amount of Hydrogen in Vent Stream Versus Split Factor')
%% EP Calculation - 18 March 2024
Keys = {'Benzene','Toluene','Hydrogene','Byphenil', 'Methane'};
Cost_sell = [12.5, 8.8, 2.1, 7.4, 2.1];
Cost_sell_Materials = containers.Map(Keys,Cost_sell);
delH_burn = [1.41, 1.68, 0.123, 2.688, 0.383];
Cost_burn_Materials = containers.Map(Keys, delH_burn*4);

% 1st Scenario
EP2_sell = 8000*(Cost_sell_Materials('Benzene')*Bb + Cost_sell_Materials('Byphenil')*Dd - Cost_sell_Materials('Hydrogene')*(F1h+F1m) - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogene')*Vh + Cost_burn_Materials('Methane')*Vm);
EP2_burn = 8000*(Cost_sell_Materials('Benzene')*Bb + Cost_burn_Materials('Byphenil')*Dd - Cost_sell_Materials('Hydrogene')*(F1h+F1m) - Cost_sell_Materials('Toluene')*F2t + Cost_burn_Materials('Hydrogene')*Vh + Cost_burn_Materials('Methane')*Vm);

% 
figure(3)
subplot(1,2,1);
plot(SFvalues,EP2_sell(1,:))
hold on
plot(SFvalues,EP2_sell(2,:))
plot(SFvalues,EP2_sell(3,:))
plot(SFvalues,EP2_sell(4,:))
grid on
xlabel('Split Factor'); ylabel('EP (M$/8000hr)')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('EP2 for sell scenario')
% figure(4)
subplot(1,2,2);
plot(SFvalues,EP2_burn(1,:))
hold on
plot(SFvalues,EP2_burn(2,:))
plot(SFvalues,EP2_burn(3,:))
plot(SFvalues,EP2_burn(4,:))
grid on
xlabel('Split Factor'); ylabel('EP (M$/8000hr)')
legend('873.15 K', '923.15 K', '973.15 K', '1023.15 K')
title('EP2 for burn scenario')

figure(4)
subplot(1,2,1);
plot(Tvalues, EP2_sell(:,1))
hold on
for i=2:40
    plot(Tvalues,EP2_sell(:,i))
end
grid on
xlabel('Temperature [K]'); ylabel('EP (M$/8000hr)')
title('EP2 for sell scenario')

% figure(6)
subplot(1,2,2);
plot(Tvalues, EP2_burn(:,1))
hold on
for i=2:40
    plot(Tvalues,EP2_burn(:,i))
end
grid on
xlabel('Temperature [K]'); ylabel('EP (M$/8000hr)')
title('EP2 for burn scenario')

