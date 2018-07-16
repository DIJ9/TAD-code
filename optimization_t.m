%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Tailsitter mission design %%%%%%%%%%%%%%%%%%%%%%%
%%%% 07/09/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Daiju uehara, Larry Matthies, Phil Tokumaru %%%%% 

clc 
clear all

% Define planet paremeters on Titan
a = 190;% speed of sound [m/s]
rho = 5.34;% air density [kg/m^3]
g = 1.352;% gravity [m/s^2]
nu = 1.2e-6;% kinematic viscosity [m^2/s]

%%
% Define vehicle parameters
m = 10;% vehicle mass [kg]
N = m*g;% vehicle weight [N]
Nr = 2;% number of rotor, assuming a quadcopter config

% Hover performance assumption
FM = 0.72;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.18;% Tip mach number
DL = 50/g;% disk loading [kg/m^2] typical for toy-scale helicopter?

% Parameters of each rotor 
Vtip = Mtip*a;% Tip speed [m/s]
A = m/Nr/DL;% disk area
R = sqrt(A/pi);% rotor radius [m]
omega = Vtip/R;% angular velocity [rad/s]
RPM = 2*pi/omega*60;% angular velocity [RPM]

% Ideal hover power 
Th = N/Nr;%thrust of each rotor for hover
vh = sqrt(DL*g/2/rho);% ideal induced velocity in hover
Ph = Th*vh;% ideal induced power in hover

% Payload fuselage design
% Payload/fuselage design
df = 300;% fuselage density [kg/m^3]
rf = (m/df*3/4*1/pi)^(1/3);
Af = pi*rf^2;% fuselage frontal area [m^2]

mp = 2;% payload mass [kg]
ma = 0.5;% avionics mass [kg]
da = 2500;% avionics density [kg/m^3]
va = ma/da;% avionics volume [m^3]
de = 2500;% energy storage density [kg/m^3]

Cdmulti = 1.5;% drag multiplier

% Efficiency coefficients
etam = 0.80;% motor efficiency 
etac = 0.95;% control efficiency 

etafw = FM*etam*etac;% total forward flight efficiency
kinduced = 1.15;% induced power factor 

% Blade parameters
dblade = 0.15;% blade mass density based on disk area [kg/m^2]
mblade = dblade*A*4;
Cdblade = 0.007;% blade drag coefficient

% Support tube mass
msupport = 0.0720*R*Nr;

% Motor constant 
m1 = 0.003;% motor constant [kg]
m2 = 0.322;% motor constant [kg/N-m]

% Battery specific energy
Ebconst = 100;% Wh/kg

%%
% Level flight analysis
V = 3:0.05:15;% forward speed flight range

mbattery = 0;
% Matrix allocation
Ss = zeros(1,length(V));
time = zeros(1,length(V));
hours = zeros(1,length(V));
range = zeros(1,length(V));
rm = zeros(1,length(V));
Ps = zeros(1,length(V));
mwings = zeros(1,length(V));
mba = zeros(1,length(V));
Thrust = zeros(1,length(V));
kwss = zeros(1,length(V));
kbs = zeros(1,length(V));
vs = zeros(1,length(V));

kw = 20;% wing weight coefficient [N/m^2]

AR = 11;
CD0 = 0.03;
e = 0.75;
K = 1/pi/AR/e;
CL = sqrt(CD0/K);
counter = 0;

for ii = 1:length(V)

        Re = V(1,ii)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end
        S = 2*N/rho/CL/(V(1,ii)^2);
        mwing = kw*S;
        kws = N/S;
        Dwing = 1/2*rho*(CD0+K*CL^2)*S*V(1,ii)^2;
        Dbody = 1/2*Cdbody*V(1,ii)^2*rho*Af;% parasite drag [N]
        Ttotal = Dwing+Dbody;
        Trotor = Ttotal/Nr;
        
        vtemp = sqrt(Trotor/2/rho/A);
        vclimb = vtemp*(-V(1,ii)/(2*vtemp)+sqrt(1+(V(1,ii)/(2*vtemp))^2));% induced velocity in climb [m/s]

        Ptemp = Trotor*vtemp;
        Pclimb = Ptemp*(V(1,ii)/(2*vtemp)+sqrt(1+(V(1,ii)/(2*vtemp))^2));% induced power in climb [W]
        Pca = Pclimb/etafw;% actual power [W]
        Pcatotal = Pca*Nr;% Total climb power 
        
        Q = Pclimb/omega;% Torque of each rotor [Nm]
        mmotor = m1+m2*Q;% motor mass [kg]
        
        % Find total mass required
        mbattery = m-mmotor*Nr-mwing-ma-mp-msupport-mblade-0.1*m;
        counter = counter+1;

        Eb = Ebconst*mbattery;% battery energy [Wh]
        time(1,ii) = Eb/Pcatotal*3600;% mission endurance [s]
        hours(1,ii) = time(1,ii)/3600;% mission endurance [h]
        range(1,ii) = time(1,ii)*V(1,ii)/1000;% mission range [km]
        rm(1,ii) = range(1,ii)/2;% mission radius [km]
        mba(1,ii) = mbattery;
        kbs(1,ii) = mbattery/m;
        Thrust(1,ii) = Ttotal;
        Ps(1,ii) = Pcatotal;
        vs(1,ii) = vclimb;
        Ss(1,ii) = S;
        mwings(1,ii) = mwing;
        kwss(1,ii) = kws;

end

%%
% figure(1)
% hold on
% si = 2;
% h5 = plot(V,mba,'linewidth',si);
% plot(V,mwings,'linewidth',si)
% % plot(V,Pps,'linewidth',si)
% % hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);
% 
% h = gca;
% grid off;
% fs = 12;
% % set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% % title('ground config');
% xlabel('airspeed V [m/s]','FontSize',fs);
% ylabel('Mass [kg] ','FontSize',fs);
% set(h,'FontName','Times New Roman','linewidth',si);
% set(h,'FontName','Times New Roman','linewidth',si);
% legend('Battery mass','Wing mass');
% 
% % a1 = 0;
% % a2 = 16;
% b1 = 0;
% b2 = 10;
% % xlim([a1 a2]);
% ylim([b1 b2]);
% 
% FS = 12;
% xt = a1:1:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:10:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([2 1 1])
%%
figure(2)
hold on
si = 2;
h5 = plot(V,Ss,'linewidth',si);
% plot(V,mwings,'linewidth',si)
% plot(V,Pps,'linewidth',si)
% hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);

h = gca;
grid off;
fs = 12;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('Wing area [m^2] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
% legend('Battery mass','Wing mass');

% a1 = 0;
% a2 = 16;
% b1 = 0;
% b2 = 10;
% % xlim([a1 a2]);
% ylim([b1 b2]);
% 
% FS = 12;
% xt = a1:1:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:10:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;
% pbaspect([2 1 1])
%%
%%
figure(2)
si = 2;
subplot(2,1,1)
hold on
hhh5 = plot(V,hours,'linewidth',si);
% hhhh5 = plot(V(1,Imin),hours(1,Imin),'ok','linewidth',si);

h = gca;
grid off;
fs = 12;

ylabel('Endurance [hours] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
ylim([0 10]);

subplot(2,1,2)
hold on
hhhhh5 = plot(V,range,'linewidth',si);
% hhhhhh5 = plot(V(1,Imin),range(1,Imin),'ok','linewidth',si);

h = gca;
grid off;
fs = 12;

xlabel('Airspeed V [m/s]','FontSize',fs);
ylabel('Range [km] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

a1 = 0;
a2 = 16;
b1 = 0;
b2 = 250;
% xlim([a1 a2]);
ylim([b1 b2]);

% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:0.5:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;


