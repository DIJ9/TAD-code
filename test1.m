clc
Vd = -6:0.01:-1;
k = 1.2;
k1 = -1.125;
k2 = -1.372;
k3 = -1.718;
k4 = -0.655;

Pstore = zeros(1,length(Vd));
Estore = zeros(1,length(Vd));
vstore = zeros(1,length(Vd));
ab = 10000;% balloon altitude [m]
dtime = abs(ab./Vd);% time to climb [s]

for ii = 1:length(Vd)
    
    if Vd(1,ii)/vh < 0 && -2 < Vd(1,ii)/vh
        vdescent = vh*(k+k1*(Vd(1,ii)/vh)+k2*(Vd(1,ii)/vh)^2+k3*(Vd(1,ii)/vh)^3+k4*(Vd(1,ii)/vh)^4);
    else
        vdescent = vh*(-Vd(1,ii)/(2*vh)-sqrt((Vd(1,ii)/(2*vh))^2-1));
    end
    vstore(1,ii) = vdescent;
    
%     Pdescent = Ph*(Vd(1,ii)/(2*vh)-sqrt((Vd(1,ii)/(2*vh))^2-1));
    Pdescent = Ph*(vdescent/vh+Vd(1,ii)/vh);
    Tdescent = -2*rho*A*(Vd(1,ii)+vdescent)*vdescent;
    
    Pda = Pdescent/FM/etac/etam;% actual power [W]
    Pdatotal = Pda*4;
    Pstore(1,ii) = Pdatotal;
    Ed = Pdatotal*dtime(1,ii)/3600;% Total climb energy 
    Estore(1,ii) = Ed;
end
plot(Vd,Pstore)