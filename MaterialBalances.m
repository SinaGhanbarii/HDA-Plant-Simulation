function eqs = MaterialBalances(unk,SF,T)
    % ------------------------------------------------------------------------------------
    % % Input data:
    x_m = 0.05;     %[-]
    s_b = 0.96;     %[-] Benzene selectivity
    p_b = 265;      %[kmol/h]
    HTR = 5;        %[-]

    % ------------------------------------------------------------------------------------
    % % Unknowns:
    F1h = unk(1);    F1m = unk(2);    F2t = unk(3);    INh = unk(4);    INm = unk(5);
    INt = unk(6);    Tt = unk(7);     Dd = unk(8);     Bb = unk(9);     RVh = unk(10);
    RVm = unk(11);   Vh = unk(12);    Vm = unk(13);    Rh = unk(14);    Rm = unk(15);

    V = unk(16);

    % ------------------------------------------------------------------------------------
    % % Equations:
    % MBi @ Mixer:
    eqs(1) = F1h + Rh - INh;            %(i=h)
    eqs(2) = F1m + Rm - INm;            %(i=m)
    eqs(3) = F2t + Tt - INt;            %(i=t)

    % MBi @ Reactor:
    Fin = [INh INm INt 0 0];
    [~,Fout] = ode45(@(V,F)Kinetics(V,F,T),[0 V],Fin);
    OUTh = Fout(end,1);
    OUTm = Fout(end,2);
    OUTt = Fout(end,3);
    OUTb = Fout(end,4);
    OUTd = Fout(end,5);

    % MBi @ Separator:
    eqs(4) = OUTh - RVh;                %(i=h)
    eqs(5) = OUTm - RVm;                %(i=m)
    eqs(6) = OUTt - Tt;                 %(i=t)
    eqs(7) = OUTb - Bb;                 %(i=b)
    eqs(8) = OUTd - Dd;                 %(i=d)

    % MBi @ Splitter:
    eqs(9) = (RVh*SF) - Vh;             %(i=h)
    eqs(10) = (RVm*SF) - Vm;            %(i=m)
    eqs(11) = RVh - Vh - Rh;            %(i=h)
    eqs(12) = RVm - Vm - Rm;            %(i=m)

    % Specifications:
    eqs(13) = F1m/(F1m + F1h) - x_m;
    eqs(14) = OUTb/(INt - OUTt) - s_b;
    eqs(15) = Bb - p_b;
    eqs(16) = INh/INt - HTR;            %Needed because we don't know the V of the reactor

end