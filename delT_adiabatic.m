clearvars
close all
clc
Tvalues = (600:50:750)+273.15; %K
SFvalues = 0.8:0.1:1; %[-]

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

    T = Tvalues(i); %[K]
    Guess = [2000 100 200 2000 100 400 200 0 265 1000 200 1000 200 0 0 100];


    for j = length(SFvalues):(-1):1

        SF = SFvalues(j); %[-]

        Solution = fsolve(@(unk)MaterialBalances(unk,SF,T),Guess);

        Guess = Solution;

        F1h(i,j) = Solution(1); %[kmol/h]
        F1m(i,j) = Solution(2); 
        F2t(i,j) = Solution(3);
        INh(i,j) = Solution(4); 
        INm(i,j) = Solution(5); 
        INt(i,j) = Solution(6);
        Tt(i,j) = Solution(7); 
        Dd(i,j) = Solution(8); 
        Bb(i,j) = Solution(9); 
        RVh(i,j) = Solution(10); 
        RVm(i,j) = Solution(11); 
        Vh(i,j) = Solution(12); 
        Vm(i,j) = Solution(13);
        Rh(i,j) = Solution(14); 
        Rm(i,j) = Solution(15);
        V(i,j) = Solution(16); %m3
        
        
    end

end

%%determine, via numerical integration of the plug‐flow model of the reactor, the conversion,
% selectivity and residence time as a function of the operating temperature, by assuming the
% reactor isothermal (HINT: at this level of detail, the presence of recycles can be neglected
% in the evaluation of the initial molar flows/concentrations).

%We call the Molar flows now 
Nin=[INh(:,end) INm(:,end) INt(:,end) zeros(length(Tvalues),1) zeros(length(Tvalues),1)].*3600; %[kmol/s] why the zeros?
Nout=[RVh(:,end) RVm(:,end) Tt(:,end) Bb(:,end) 0*Dd(:,end)];
for i=1:length(Tvalues)

    no=Nin(i,:); %line vector
        
    Tin=Tvalues(i);
    xo=Tin;
    Toutsolution=fsolve(@(Tout)DeltaT_ad(Tout,no,Tin),xo);
    DeltaTadiabitc(i,:)=Toutsolution-Tin;

end

for i=1:length(Tvalues)
        T=Tvalues(i); %[K]

    figure(1)

plot(Tvalues-273.15,DeltaTadiabitc,'or')
hold on
xlabel("Temperature $^\circ C$",'Interpreter','latex')
ylabel("Adiabatic $\Delta T$ [$^\circ C$]",'Interpreter','latex')
grid on
end
function F=DeltaT_ad(Tout,no,Tin)

%Cp [KJ/(kmol*K)
cp_h=[25.399 2.0178e-02 -3.8549e-05 3.188e-08 -8.7585e-12]; % T range[250 1500]K
cp_m=[34.942 -3.9957e-02 1.9184e-04 -1.5303e-07 3.9321e-11]; % T range [50 1500]K
cp_t=[-24.097 5.2187e-01 -2.9827e-04 6.122e-08 1.2576e-12]; % T range [200 1500]K
cp_b=[-31.368 4.746e-01 -3.1137e-04 8.5237e-08 -5.0524e-12]; % T range[200 1500]K
cp_d=[-29.153 7.6716e-01 -3.4341e-04 -3.7724e-08 4.6179e-11]; % T range [200 1500]K
cp=[cp_h;cp_m;cp_t;cp_b;cp_d];
a=cp(:,1); %column vector?
b=cp(:,2);
c=cp(:,3);
d=cp(:,4);
e=cp(:,5);

Deltahf=[0 -74.85 50.00 82.93 182.09]./1e-3; %[kJ/kmol]
hVT=@(x) Deltahf'+a.*(x-298)+b./2.*(x^2-298^2)+c./3.*(x^3-298^3)+d./4.*(x^4-298^4)+e./5.*(x^5-298^5); %vettore colonna

hVTin=hVT(Tin);
hVTout=hVT(Tout);

Hin=no*hVTin;

nu=[-1 1 -1 1 0];
l=no(3); %total conversion of toluene
n=no+nu.*l;
Hout=n*hVTout;

F=Hin-Hout;


end
