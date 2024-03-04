function dFdV = Kinetics(V,F)

Fh = F(1); % [kmol/hr]
Fm = F(2); % [kmol/hr]
Ft = F(3); % [kmol/hr]
Fb = F(4); % [kmol/hr]
Fd = F(5); % [kmol/hr]

Ftot = sum([Fh Fm Ft Fb Fd]);

% Data
P = 34E+5; % Pa
A1 = 3.5E+10; % [m^(1.5)/kmol^(0.5)*s]
A2 = 2.1E+12; % [m^3/kmol/s]
E1 = 50900; %[kcal/kmol]
E2 = 60500; %[kcal/kmol]

Rgas1 = 8314; %[m^3.Pa/(mol.K)]
Rgas2 = 1.987; %[kcal/kmol.K]

Ch = (P/Rgas1/T)* (Fh/Ftot);
Ct = (P/Rgas1/T)* (Ft/Ftot);
Cb = (P/Rgas1/T)* (Fb/Ftot);



k1 = A1*exp(- E1/(Rgas2*T));
k2 = A2*exp(- E2/(Rgas2*T));

R1 = 3600*k1*Ct*sqrt(Ch);
R2 = 3600*k2* Cb^2;

dFhdV = -R1+R2;
dFmdV = +R1;
dFtdV = -R1;
dFbdV = +R1 - 2*R2;
dFddV = +R2;

dFdV = [dFhdV dFmdV dFtdV dFbdV dFddV]';

end