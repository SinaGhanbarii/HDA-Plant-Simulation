function eqs = MaterialBalances(unk)
% Input Data:
x_M = 0.05; % Mole fraction of Methane
s_b = 0.96; % Benzene Selectivity
p_b = 265; % Benzene amount in P [kmol/h]
HTR = 5; % Hydrogen to Toluene ratio

% Unknowns:
F1h = unk(1);
F1m = unk(2);
F2t = unk(3);
INh = unk(4);
INm = unk(5);
INt = unk(6);
Tt = unk(7);
Dd = unk(8);
Bb = unk(9);
RVh = unk(10);
RVm = unk(11);
Vh = unk(12);
Vm = unk(13);
Rh = unk(14);
Rm = unk(15);
V = unk(16);
% Equations:

% Material Balance for mixer:
eqs(1) = F1h + Rh - INh; %(i=h)
eqs(2) = F1m + Rm - INm; %(i=m)
eqs(3) = F2t + Tt - INt; %(i=t)

% Material Balance for Reactor:
Fin = [INh INm INt 0 0];
[~, Fout] = ode45(@(V,F)Kinetics(V,F,T),[0 V], Fin);

% ??????
OUT = ones(1,5);
for i = 1: length(OUT)
    OUT(i) = Fout(end,i);
end

OUTh = Fout(end,1);
OUTm = Fout(end,2);
OUTt = Fout(end,3);
OUTb = Fout(end,4);
OUTd = Fout(end,5);

% Material Balance for seperation system:
eqs(4) = OUTh - RVh;
eqs(5) = OUTm - RVm;
eqs(6) = OUTt - RVt;
eqs(7) = OUTb - Bb;
eqs(8) = OUTd - Dd;


% Material Balance for splitter:
eqs(9) = RVh;
eqs(10) = RVm*SF - Vm;
eqs(11) = RVh - Vh - Rh;
eqs(12) = RVm - Vm - Rm;

% Spec#1
eqs(13) = F1m/(F1m+F1h)- x_M;
% Spec#2
eqs(14) = OUTb/(INt-OUTt) - s_b;
% Spec#3
eqs(15) = Bb - p_b;
% Spec#4
eqs(16) = INh/INt - HTR;

end