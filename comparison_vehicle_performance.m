%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Tailsitter multi-copter comparison %%%%%%%%%%%%%%
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
Nrt = 2;% number of rotor
Nrm = 4;

% Hover performance assumption
FM = 0.72;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.18;% Tip mach number
DL = 50/g;% disk loading [kg/m^2] typical for toy-scale helicopter?

% Parameters of each rotor 
Vtip = Mtip*a;% Tip speed [m/s]

At = m/Nrt/DL;% disk area
Am = m/Nrm/DL;% disk area
Rt = sqrt(At/pi);% rotor radius [m]
Rm = sqrt(Am/pi);% rotor radius [m]
omegat = Vtip/Rt;% angular velocity [rad/s]
omegam = Vtip/Rm;% angular velocity [rad/s]

RPMt = 2*pi/omegat*60;% angular velocity [RPM]
RPMm = 2*pi/omegam*60;% angular velocity [RPM]

% Ideal hover power 
Tht = N/Nrt;%thrust of each rotor for hover
vht = sqrt(DL*g/2/rho);% ideal induced velocity in hover
Pht = Tht*vht;% ideal induced power in hover
Thm = N/Nrm;%thrust of each rotor for hover
vhm = sqrt(DL*g/2/rho);% ideal induced velocity in hover
Phm = Thm*vhm;% ideal induced power in hover

% Payload/fuselage design
df = 300;% fuselage density [kg/m^3]
% rf = (m/df*3/4*1/pi)^(1/3);
rf = 0.050;
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
mbladet = dblade*At*4;
mbladem = dblade*Am*4;
Cdblade = 0.007;% blade drag coefficient

% Support tube mass
msupport = 0;

% Motor constant 
m1 = 0.003;% motor constant [kg]
m2 = 0.322;% motor constant [kg/N-m]

% Battery specific energy
Ebconst = 100;% Wh/kg

%%
% Level flight analysis
V = 1:0.05:20;% forward speed flight range

mbattery = 0;
% Matrix allocation
timet = zeros(1,length(V));
timem = zeros(1,length(V));

hourst = zeros(1,length(V));
hoursm = zeros(1,length(V));

ranget = zeros(1,length(V));
rangem = zeros(1,length(V));

Pst = zeros(1,length(V));
Psm = zeros(1,length(V));

mbat = zeros(1,length(V));
mbam = zeros(1,length(V));

Thrustt = zeros(1,length(V));
Thrustm = zeros(1,length(V));

Ss = zeros(1,length(V));
kwss = zeros(1,length(V));
% kbs = zeros(1,length(V));
% vs = zeros(1,length(V));

kw = 20;% wing weight coefficient [N/m^2]
AR = 6;
CD0 = 0.03;
e = 0.75;
K = 1/pi/AR/e;
CL = sqrt(CD0/K);
counter = 0;

for ii = 1:length(V)
% Tail sitter
        Re = V(1,ii)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.3;
        end
        S = 2*N/rho/CL/(V(1,ii)^2);
        mwing = kw*S;
        kws = N/S;
        Dwing = 1/2*rho*(CD0+K*CL^2)*S*V(1,ii)^2;
        Dbody = 1/2*Cdbody*V(1,ii)^2*rho*Af;% parasite drag [N]
        Ttotal = Dwing+Dbody;
%         Ttotal = Dwing;
        Trotor = Ttotal/Nrt;
        
        vtemp = sqrt(Trotor/2/rho/At);
        vclimb = vtemp*(-V(1,ii)/(2*vtemp)+sqrt(1+(V(1,ii)/(2*vtemp))^2));% induced velocity in climb [m/s]

        Ptemp = Trotor*vtemp;
        Pclimb = Ptemp*(V(1,ii)/(2*vtemp)+sqrt(1+(V(1,ii)/(2*vtemp))^2));% induced power in climb [W]
        Profile = 20;

        Pca = Pclimb/etafw;% actual power [W]
        Pcatotal = Pca*Nrt+Profile;% Total climb power 
        
        Q = Pca/omegat;% Torque of each rotor [Nm]
        mmotor = m1+m2*Q;% motor mass [kg]
        
        % Find total mass required
        mbattery = m-mmotor*Nrt-mwing-ma-mp-msupport-mbladet-0.1*m;
        counter = counter+1;

        Eb = Ebconst*mbattery;% battery energy [Wh]
        timet(1,ii) = Eb/Pcatotal*3600;% mission endurance [s]
        hourst(1,ii) = timet(1,ii)/3600;% mission endurance [h]
        ranget(1,ii) = timet(1,ii)*V(1,ii)/1000;% mission range [km]
        mbat(1,ii) = mbattery;
%         kbs(1,ii) = mbattery/m;
        Thrustt(1,ii) = Ttotal;
        Pst(1,ii) = Pcatotal;
%         vs(1,ii) = vclimb;
        Ss(1,ii) = S;
%         mwings(1,ii) = mwing;
        kwss(1,ii) = kws;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             evolume = mbattery/de;% energy storage volume [m^3]
%             wbvolume = evolume+va;% warm box volume [m^3]
%             wblength = (wbvolume)^(1/3);% warm box side length [m]
%             Awb = (wblength+it)^2;% warm box frontal area with insulation [m^2]
%             Afuselage = Awb+Ap;% payload + warm box [m^2]
%             lre = Afuselage/pw;% characteristic length for Re [m]
%             Re = V(1,ii)*2*rf/nu;% Reynolds number
%             if Re < 2e5
%                 Cdbody = 0.47;% drag coefficient of vehicle body
%             else
%                 Cdbody = 0.18;
%             end
            Dbody = 1/2*Cdmulti*Cdbody*V(1,ii)^2*rho*Af;% parasite drag [N]
            Ppara = Dbody*V(1,ii);% Parasite power [W]
            Drotor = Dbody/Nrm;% drag per rotor [N]
            AoA = atan(Drotor/Thm);% disk angle of attack [rad]
            degree = rad2deg(AoA);% AoA [deg]
            Ttotal = sqrt(Drotor^2+Thm^2);% Total thrust of each rotor [N]
            mu = V(1,ii)*cos(AoA)/(Vtip);
            Ct = Ttotal/(rho*Am*Vtip^2);
            lambda = sqrt(Ct/2);
            lambda_old = 10*lambda;
            iter = 0;
            while abs(lambda_old-lambda) > 10^-3 && lambda ~= 0
                lambda_old = lambda;
                lambda = lambda - (lambda-mu*tan(AoA)-Ct/2/sqrt(mu^2+lambda^2))...
                         /(1+Ct/2*lambda/((mu^2+lambda^2)^(3/2)));
                iter = iter + 1;
            end
            vi = lambda*Vtip-V(1,ii)*sin(AoA);% induced velocity in forward flight
            Pi = Ttotal*vi;% induced power without induced loss
            Paero = (Pi*Nrm+Ppara)/etafw+Profile;% no profile power because of the momentum theory based analysis
            Paerorotor = Paero/Nrm;% shaft power of each rotor 
%             mskin = 1.600*8*Afuselage;% Fuselage mass
%             rch = Pi/Ph;% ratio of forward flight induced power to hover 
            Q = Paerorotor/omegam;% Torque of each rotor [Nm]
            mmotor = m1+m2*Q;% motor mass [kg]

            % Find total mass required
            mbattery = m-mmotor*Nrm-ma-mp-msupport-mbladem-0.1*m;
            counter = counter+1;

          
        Eb = Ebconst*mbattery;% battery energy [Wh]
        timem(1,ii) = Eb/Paero*3600;% mission endurance [s]
        hoursm(1,ii) = timem(1,ii)/3600;% mission endurance [h]
        rangem(1,ii) = timem(1,ii)*V(1,ii)/1000;% mission range [km]
%         Qs(1,ii) = Q;
%         alpha(1,ii) = degree;
        mbam(1,ii) = mbattery;
        Thrustm(1,ii) = Ttotal;
        Psm(1,ii) = Paero;
%         Pis(1,ii) = Pi*Nr/etafw;
%         Pps(1,ii) = Ppara;
%         vis(1,ii) = vi;
end
%%
figure(1)
hold on
si = 2;
h5 = plot(V,mbam,'linewidth',si);
plot(V,mbat,'linewidth',si)
% plot(V,Pps,'linewidth',si)
% hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);

h = gca;
grid off;
fs = 16;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('battery mass [kg] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
legend('Multi-copter','Tail-sitter');

% a1 = 0;
% a2 = 16;
b1 = 0;
b2 = 7;
% xlim([a1 a2]);
ylim([b1 b2]);

FS = 12;

%%
figure(2)
si = 2;
subplot(2,1,1)
hold on
hhh5 = plot(V,hoursm,'linewidth',si);
hhhh5 = plot(V,hourst,'linewidth',si);

h = gca;
grid off;
fs = 12;

ylabel('Endurance [hours] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
ylim([0 Inf]);

subplot(2,1,2)
hold on
hhhhh5 = plot(V,rangem,'linewidth',si);
hhhhh5 = plot(V,ranget,'linewidth',si);

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
b2 = 500;
% xlim([a1 a2]);
ylim([b1 b2]);

%%
figure(3)
hold on
si = 2;
h5 = plot(V,Psm,'linewidth',si);
plot(V,Pst,'linewidth',si)
% plot(V,Pps,'linewidth',si)
% hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);

h = gca;
grid off;
fs = 16;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('Total power [W] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
legend('Multi-copter','Tail-sitter');

% a1 = 0;
% a2 = 16;
b1 = 0;
b2 = 7;
% xlim([a1 a2]);
% ylim([b1 b2]);
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;