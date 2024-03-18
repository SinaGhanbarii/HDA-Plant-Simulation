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