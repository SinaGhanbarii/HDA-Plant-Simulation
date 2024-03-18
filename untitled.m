clear all, close all, clc
% DOFs are Split Factor (SP), which is in range (0,1] and T, Which is in range [600 650 700 750] K
Tvalues = (600:50:750) + 273.15; % Temperature range [K]
SFvalues = 0.5:0.1:1;

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

for i= 1:length(Tvalues)
    T = Tvalues(i);
    % Guess values for SF = 100%
    guess = [2000 100 200 2000 100 400 200 6 265 1000 200 1000 200 0 0 1000];
    for j=length(SFvalues):-1:1
        SF = SFvalues(j);

        solution = fsolve(@(unk)MaterialBalances(unk,SF,T),guess);
        F1h(i,j) = solution(1); % [kmol/hr]
        F1m(i,j) = solution(2); % [kmol/hr]
        F2t(i,j) = solution(3); % [kmol/hr]
        INh(i,j) = solution(4); % [kmol/hr]
        INm(i,j) = solution(5); % [kmol/hr]
        INt(i,j) = solution(6); % [kmol/hr]
        Tt(i,j) = solution(7); % [kmol/hr]
        Dd(i,j) = solution(8); % [kmol/hr]
        Bb(i,j) = solution(9); % [kmol/hr]
        RVh(i,j) = solution(10); % [kmol/hr]
        RVm(i,j) = solution(11); % [kmol/hr]
        Vh(i,j) = solution(12); % [kmol/hr]
        Vm(i,j) = solution(13); % [kmol/hr]
        Rh(i,j) = solution(14); % [kmol/hr]
        Rm(i,j) = solution(15); % [kmol/hr]
        V(i,j) = solution(16); % [m^3]

        guess = solution;
    end
end

R = Rh + Rm; %[kmol/h]
plot(SF,values)
xlabel('Split Factor')
ylabel  ('Recycle, R [kmol/hr]')