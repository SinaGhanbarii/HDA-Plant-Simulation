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