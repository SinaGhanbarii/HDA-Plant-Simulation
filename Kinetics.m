function dFdV = Kinetics(V,F,T)
    Fh = F(1);                      %[kmol/h]
    Fm = F(2);                      %[kmol/h]
    Ft = F(3);                      %[kmol/h]
    Fb = F(4);                      %[kmol/h]
    Fd = F(5);                      %[kmol/h]

    Ftot = Fh + Fm + Ft + Fb + Fd;  %[kmol/h]

    % --------------------------------------------------------------------
    % % Data:
    P = 34E+5;          %[Pa]
    A1 = 3.5E+10;       %[m^(1.5)/kmol^(0.5)/s]
    A2 = 2.1E+12;       %[m^3/kmol/s]
    E1 = 50900;         %[kcal/kmol]
    E2 = 60500;         %[kcal/kmol]

    Rgas1 = 8314;       %[m^3*Pa/kmol/K]
    Rgas2 = 1.987;      %[kcal/kmol/K]
    
    % --------------------------------------------------------------------
    % % Calculations:
    % Concentrations: PV*MW = nRT*MW     ==>>    Ci = PMW/RT
    Ch = (P/Rgas1/T)*(Fh/Ftot);      %[kmol/m^3]
    Ct = (P/Rgas1/T)*(Ft/Ftot);      %[kmol/m^3]
    Cb = (P/Rgas1/T)*(Fb/Ftot) ;     %[kmol/m^3]

    k1 = A1*exp(-E1/Rgas2/T);           %[m^(1.5)/kmol^(0.5)/s]
    k2 = A2*exp(-E2/Rgas2/T);           %[m^3/kmol/s]

    R1 = 3600*k1*Ct*sqrt(Ch);           %[m^(1.5)/kmol^(0.5)/h]
    R2 = 3600*k2*(Cb^2);                %[m^3/kmol/h]

    dFhdV = -R1 + R2;
    dFmdV = +R1;
    dFtdV = -R1;
    dFbdV = +R1 - 2*R2;
    dFddV = +R2;

    dFdV = [dFhdV; dFmdV; dFtdV; dFbdV; dFddV];
end
