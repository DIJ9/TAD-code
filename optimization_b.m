%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Balloon mission design %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 06/20/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Daiju uehara, Larry Matthies, Phil Tokumaru %%%%% 

clc 
clear all
%% Define planet paremeters on Titan
a = 190;% speed of sound [m/s]
rho = 5.34;% air density [kg/m^3]
g = 1.352;% gravity [m/s^2]
nu = 1.2e-6;% kinematic viscosity [m^2/s]

%% Define variables
m = 4:2:10;% vehicle mass [kg]
V = 2:0.1:14;% forward speed flight range
Vc = 1:0.1:10;% climb speed range
Vd = -6:0.01:-1;% descent speed range

%% Matrix allocation
% forward flight
timef = zeros(length(m),length(V));
rangef = zeros(length(m),length(V));
Pf = zeros(length(m),length(V));
Pif = zeros(length(m),length(V));
Ppf = zeros(length(m),length(V));
alphaf = zeros(length(m),length(V));
Thrustf = zeros(length(m),length(V));
Ef = zeros(length(m),length(V));

% climb 
rangec = zeros(length(m),length(Vc));
Pc = zeros(length(m),length(Vc));
Pic = zeros(length(m),length(Vc));
Ppc = zeros(length(m),length(Vc));
alphac = zeros(length(m),length(Vc));
Thrustc = zeros(length(m),length(Vc));
Ec = zeros(length(m),length(Vc));
vc = zeros(length(m),length(Vc));

% descent
vd = zeros(length(m),length(Vc));
Pd = zeros(length(m),length(Vc));
Ed = zeros(length(m),length(Vc));

% mass margin
marg = zeros(1,length(m));
mba = zeros(1,length(m));
Eb = zeros(1,length(m));
Eba = zeros(1,length(m));
time = zeros(1,length(m));


%% Rotor constant
Nr = 4;% number of rotor, assuming a quadcopter config
FM = 0.75;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.25;% Tip mach number
DL = 50/g;% disk loading of each rotor [kg/m^2] 

%% Efficiency coefficients
etam = 0.85;% motor efficiency 
etac = 0.95;% control efficiency 
kinduced = 1.15;% induced power factor 
etafw = FM*etam*etac;% total forward flight efficiency
    
%% Sampling state analysis
tsample = 5*60;% sampling time [s]
Psample = 20;% sampling power [W]

%% Balloon parameter
ab = 10000;% balloon altitude [m]
Vbal = 2;% wind velocity [m/s]

%% Vehicle parameters
mp = 2;% payload mass [kg]
ma = 0.5;% avionics mass [kg]

Cdmulti = 1.5;% drag multiplier

m1 = 0.003;% motor constant [kg]
m2 = 0.322;% motor constant [kg/N-m]

Ebconst = 100;% Wh/kg

%%

for jj = 1:length(m)
    N = m(1,jj)*g;% vehicle weight [N]

    % Parameters of each rotor 
    Vtip = Mtip*a;% Tip speed [m/s]
    A = m(1,jj)/Nr/DL;% disk area
    R = sqrt(A/pi);% rotor radius [m]
    omega = Vtip/R;% angular velocity [rad/s]
    RPM = 2*pi/omega*60;% angular velocity [RPM]
    dblade = 0.15;% blade mass density based on disk area [kg/m^2]
    mblade = dblade*A*4;

    % Ideal hover power 
    Th = N/Nr;%thrust of each rotor for hover
    vh = sqrt(DL*g/2/rho);% ideal induced velocity in hover
    Ph = Th*vh;% ideal induced power in hover per rotor 
    Phtotal = Ph*Nr;% Total induced power in hover

    % Actual hover power 
    Pha = Ph/etafw;% actual hover power consumption [W]
    Phatotal = Nr*Pha;% Total actual power in hover
    Esample = (Psample*tsample+Phatotal*tsample)/3600;% sampling energy [Wh]

    % Payload/fuselage design
    df = 300;% fuselage density [kg/m^3]
    rf = (m(1,jj)/df*3/4*1/pi)^(1/3);
    Af = pi*rf^2;% fuselage frontal area [m^2]

    % Energy for docking
    dtime = 5*60;% time to dock [s]
    Edock = dtime*Phatotal/3600;
    
    %% Descent state (windmill state Vd<-2vh)

    k = 1.0;
    k1 = -1.125;
    k2 = -1.372;
    k3 = -1.718;
    k4 = -0.655;
    dtime = abs(ab./Vd);% time to descent [s]

    for ii = 1:length(Vd)  
        Re = Vd(1,ii)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end 
        Ddescent = 1/2*Cdmulti*Cdbody*Vd(1,ii)^2*rho*Af;% parasite drag [N]
        Tdescent = (N)/4;
        vtemp = sqrt(Tdescent/2/rho/A);
        if Vd(1,ii)/vtemp < 0 && -2 < Vd(1,ii)/vtemp
            vdescent = vtemp*(k+k1*(Vd(1,ii)/vtemp)+k2*(Vd(1,ii)/vtemp)^2+...
                       k3*(Vd(1,ii)/vtemp)^3+k4*(Vd(1,ii)/vtemp)^4);
        else
            vdescent = vtemp*(-Vd(1,ii)/(2*vtemp)-sqrt((Vd(1,ii)/(2*vtemp))^2-1));
        end
        vd(jj,ii) = vdescent;

        Ptemp = Tdescent*vtemp;
        Pdescent = Ptemp*(vdescent/vtemp+Vd(1,ii)/vtemp);
%         Tdescent = -2*rho*A*(Vd(1,ii)+vdescent)*vdescent;

        Pda = Pdescent/etafw;% actual power [W]
        Pdatotal = Pda*4;
        Pd(jj,ii) = Pdatotal;
        
        Edtemp = Pdatotal*dtime(1,ii)/3600;% Total climb energy 
        Ed(jj,ii) = Edtemp;
    end
    Vdref = -3;
    [diff,ind] = min(abs(Vdref-Vd));
    Edpicked = Ed(jj,ind);
    Pdpicked = Pd(jj,ind);
    timedpicked = dtime(1,ind);
    
    %% Climb state
    ctime = ab./Vc;% time to climb [s]

    for kk = 1:length(Vc)
        Re = Vc(1,kk)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end 
        Dclimb = 1/2*Cdmulti*Cdbody*Vc(1,kk)^2*rho*Af;% parasite drag [N]
        Tclimb = (Dclimb+N)/4;
        vtemp = sqrt(Tclimb/2/rho/A);
        vclimb = vtemp*(-Vc(1,kk)/(2*vtemp)+sqrt(1+(Vc(1,kk)/(2*vtemp))^2));% induced velocity in climb [m/s]
        vc(jj,kk) = vclimb;
        
        Ptemp = Tclimb*vtemp;
        Pclimb = Ptemp*(Vc(1,kk)/(2*vtemp)+sqrt(1+(Vc(1,kk)/(2*vtemp))^2));% induced power in climb [W]
        Pca = Pclimb/etafw;% actual power [W]
        Pcatotal = Pca*4;% Total climb power 
        Pc(jj,kk) = Pcatotal;
        
        Ectemp = Pcatotal*ctime(1,kk)/3600;% Total climb energy 
        Ec(jj,kk) = Ectemp;
    end
    [Ecpicked,ind] = min(Ec(jj,:));
    Pcpicked = Pc(jj,ind);
    timecpicked = ctime(1,ind);
    
    %% Forward flight state
    
    % Balloon travel distance during descend, climb, sampling
    dbal = 2*(timedpicked+timecpicked+tsample);% [m]
    
    for ll = 1:length(V)
        Re = V(1,ll)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end    
        Dbody = 1/2*Cdmulti*Cdbody*V(1,ll)^2*rho*Af;% parasite drag [N]
        Ppara = Dbody*V(1,ll);% Parasite power [W]
        Drotor = Dbody/Nr;% drag per rotor [N]
        AoA = atan(Drotor/Th);% disk angle of attack [rad]
        degree = rad2deg(AoA);% AoA [deg]
        Ttotal = sqrt(Drotor^2+Th^2);% Total thrust of each rotor [N]
        mu = V(1,ll)*cos(AoA)/(Vtip);
        Ct = Ttotal/(rho*A*Vtip^2);
        lambda = sqrt(Ct/2);
        lambda_old = 10*lambda;
        iter = 0;
        while abs(lambda_old-lambda) > 10^-3 && lambda ~= 0
            lambda_old = lambda;
            lambda = lambda - (lambda-mu*tan(AoA)-Ct/2/sqrt(mu^2+lambda^2))...
                     /(1+Ct/2*lambda/((mu^2+lambda^2)^(3/2)));
            iter = iter + 1;
        end
        vi = lambda*Vtip-V(1,ll)*sin(AoA);% induced velocity in forward flight
        Pi = Ttotal*vi;% induced power without induced loss
        Profile = 0;% Assumed profile power
        Paero = Pi*Nr/etafw+Ppara+Profile;
        Paerorotor = Paero/Nr;% shaft power of each rotor 
        tlevel = dbal/(V(1,ll)-Vbal);% level flight duration to catch balloon
        timef(jj,ll) = tlevel;
        Eftemp = Paero*tlevel/3600;% Total forward flight energy
        Ef(jj,ll) = Eftemp;

        Thrustf(jj,ll) = Ttotal;
        alphaf(jj,ll) = degree;
        Pf(jj,ll) = Paero;
        Pif(jj,ll) = Pi*Nr/etafw;
        Ppf(jj,ll) = Ppara;
        
    end
    
    [Efpicked,ind] = min(Ef(jj,:));
    Pfpicked = Pf(jj,ind);
    timefpicked = timef(jj,ind);
    %%     
    % Find maximum power among hover, climb, and forward flight
    [Pmotor,I] = max([Pha Pfpicked/Nr Pcpicked/Nr Pdpicked/Nr]);

    % Compute motor mass
    Q = Pmotor/omega;% Torque of each rotor [Nm]
    mmotor = m1+m2*Q;% motor mass [kg]

    % Compute battery mass 
    mbattery = m(1,jj)-mp-ma-0.2*m(1,jj)-mmotor*Nr-mblade;
    mba(1,jj) = mbattery;

    % Compute energy
    Ebtemp = Ebconst*mbattery;% available battery energy [Wh]
    Eb(1,jj) = Ebtemp;
    Ebatemp = Ecpicked+Efpicked+Edock+Edpicked+Esample;% Total estimated energy
    Eba(1,jj) = Ebatemp;
    margin = (Ebtemp-Ebatemp)/Ebconst/m(1,jj)*100;

    marg(1,jj) = margin;
%     time(1,jj) = Eb/mean([Pmotor])*3600/3600;% mission endurance [s]
%     range(jj,ii) = time(jj,ii)*V(1,ii)/1000;% mission range [km]
%     rm(jj,ii) = range(jj,ii)/2;% mission radius [km]

    
end 
%%
[X,Y] = meshgrid(Vd,m);
figure(1)
contourf(X,Y,Ed,50,'edgecolor','none');
colorbar
fs = 12;
h = gca;
xlabel('Descent rate Vd[m/s]','FontSize',fs);
ylabel('Vehicle mass [kg] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

%%
[X,Y] = meshgrid(Vc,m);
figure(2)
contourf(X,Y,Ec,50,'edgecolor','none');
colorbar
fs = 12;
h = gca;
xlabel('Climb rate Vd[m/s]','FontSize',fs);
ylabel('Vehicle mass [kg] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
%%
figure(4)
hold on
si = 2;
plot(m,marg,'linewidth',si)
fs = 12;
h = gca;
xlabel('Vehicle mass [kg]','FontSize',fs);
ylabel('Mass margin [%]','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
%%
% figure(2)
% hold on
% si = 2;
% h1 = plot(V,Ps,'linewidth',si);
% plot(V,Pis,'linewidth',si)
% plot(V,Pps,'linewidth',si)
% hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);

%%
figure(2)
hold on
si = 2;
plot(Vd,vd(3,:),'linewidth',si)
fs = 12;
h = gca;
xlabel('Descent rate Vd[m/s]','FontSize',fs);
ylabel('Induced velocity [m/s] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
ylim([0.5 6]);

%%
figure(3)
hold on
si = 2;
plot(Vd,Pd(3,:),'linewidth',si)
fs = 12;
h = gca;
xlabel('Descent rate Vd[m/s]','FontSize',fs);
ylabel('Induced power [W] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
% legend('Total','Induced','Body drag');

% a1 = 0;
% a2 = V(1,end);
% b1 = 0;
% b2 = 80;
% xlim([a1 a2]);
ylim([-60 50]);

% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:5:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([3 1 1])

%%
figure(4)
hold on
si = 2;
plot(Vc,Ec(5,:),'linewidth',si)
fs = 12;
h = gca;
xlabel('Climb rate Vd[m/s]','FontSize',fs);
ylabel('Energy consumption [Wh] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
% legend('Total','Induced','Body drag');

% a1 = 0;
% a2 = V(1,end);
% b1 = 0;
% b2 = 80;
% xlim([a1 a2]);
% ylim([-60 50]);

% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:5:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([3 1 1])
%%
figure(5)
hold on
si = 2;
plot(V,Pf(5,:),'linewidth',si)
plot(V,Pif(5,:),'linewidth',si)
plot(V,Ppf(5,:),'linewidth',si)
h = gca;
grid off;
fs = 12;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('Power [W] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
legend('Total','Induced','Body drag');

a1 = V(1,1);
a2 = 12;
b1 = 0;
b2 = 100;
xlim([a1 a2]);
ylim([b1 b2]);

FS = 12;
xt = a1:2:a2 ;
set(h,'XTick',xt) ;
set(h,'XTickLabel',xt,'fontsize',FS) ;
yt = b1:20:b2 ; 
set(h,'YTick',yt) ;
set(h,'YTickLabel',yt,'fontsize',FS) ;
pbaspect([3 1 1])

%%
figure(6)
hold on
si = 2;
plot(V,Ef(5,:),'linewidth',si)
% plot(V,Pif(5,:),'linewidth',si)
% plot(V,Ppf(5,:),'linewidth',si)
h = gca;
grid off;
fs = 12;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('Energy [Wh] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
% legend('Total','Induced','Body drag');

a1 = V(1,1);
a2 = 12;
b1 = 0;
b2 = 80;
xlim([a1 a2]);
ylim([b1 b2]);

% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:20:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([3 1 1])
%%
% figure(2)
% hold on
% si = 2;
% subplot(2,1,1)
% plot(V,hours,'linewidth',si)
% h = gca;
% grid off;
% fs = 12;
% 
% ylabel('Endurance [hours] ','FontSize',fs);
% set(h,'FontName','Times New Roman','linewidth',si);
% set(h,'FontName','Times New Roman','linewidth',si);
% 
% subplot(2,1,2)
% plot(V,range,'linewidth',si)
% h = gca;
% grid off;
% fs = 12;
% 
% xlabel('Airspeed V [m/s]','FontSize',fs);
% ylabel('Range [km] ','FontSize',fs);
% set(h,'FontName','Times New Roman','linewidth',si);
% set(h,'FontName','Times New Roman','linewidth',si);
% 
% a1 = 0;
% a2 = V(1,end);
% b1 = 0;
% b2 = 5.5;
% xlim([a1 a2]);
% % ylim([b1 b2]);
% 
% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:0.5:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
%%
% figure(3)
% hold on
% si = 2;
% plot(V,marg,'linewidth',si)
% 
% h = gca;
% grid off;
% fs = 12;
% % set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% % title('ground config');
% xlabel('airspeed V [m/s]','FontSize',fs);
% ylabel('mass margin [%]','FontSize',fs);
% set(h,'FontName','Times New Roman','linewidth',si);
% set(h,'FontName','Times New Roman','linewidth',si);
% 
% a1 = 0;
% a2 = V(1,end);
% b1 = 0;
% b2 = 50;
% xlim([a1 a2]);
% ylim([b1 b2]);
% 
% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:20:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([2 1 1])


