%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Lander mission design %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 06/18/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Daiju uehara, Larry Matthies, Phil Tokumaru %%%%% 


clc 
clear all

% Define planet paremeters on Titan
a = 195;% speed of sound [m/s]
rho = 5.428;% air density [kg/m^3]
g = 1.35;% gravity [m/s^2]
nu = 1.23e-6;% kinematic viscosity [m^2/s]

%%
% Define vehicle parameters
m = 320;% vehicle mass [kg]
N = m*g;% vehicle weight [N]
Nr = 4;% number of rotor, assuming a quadcopter config

% Hover performance assumption
FM = 0.72;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.18;% Tip mach number
% DL = 50/g;% disk loading [kg/m^2] typical for toy-scale helicopter?

% Parameters of each rotor 
Vtip = Mtip*a;% Tip speed [m/s]
R = 0.65;% rotor radius [m]
A = pi*R^2;% disk area
omega = Vtip/R;% angular velocity [rad/s]
RPM = 2*pi/omega*60;% angular velocity [RPM]

% Ideal hover power 
Th = N/Nr;%thrust of each rotor for hover
vh = sqrt(Th/2/rho/A);% ideal induced velocity in hover
Ph = Th*vh;% ideal induced power in hover
Cth = Th/(rho*A*Vtip^2);

% Payload fuselage design
it = 0.06;% insulation thickness [m]
pl = 0.17;% payload length [m]
pw = 0.12;% payload width [m]
ph = 0.07;% payload height [m]
Ap = 0.012;% payload frontal area [m^2] or pw x ph 
mp = 200;% payload mass [kg]
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
% Cdblade = 0.007;% blade drag coefficient

% Support tube mass
% msupport = 0.0720*R*Nr;

% Motor constant 
m1 = 0.003;% motor constant [kg]
m2 = 0.322;% motor constant [kg/N-m]

% Battery specific energy
Ebconst = 100;% Wh/kg
mbattery = 30;
%%
% Forward flight analysis
V = 0.1:0.1:18;% forward speed flight range
% Optimize disk angle of attack till maximum mission radius is found
% AoAop = zeros(1,length(V));
alpha = zeros(1,length(V));
time = zeros(1,length(V));
hours = zeros(1,length(V));
range = zeros(1,length(V));
rm = zeros(1,length(V));
Ps = zeros(1,length(V));
Qs = zeros(1,length(V));
mba = zeros(1,length(V));
Thrust = zeros(1,length(V));

Eb = Ebconst*mbattery;% battery energy [Wh]
Afuselage = 0.5;% frontal area
Cdbody = 0.5;% drag coefficient of vehicle body
counter = 0;

for ii = 1:length(V)
    Dbody = 1/2*Cdbody*V(1,ii)^2*rho*Afuselage;% parasite drag [N]
    Ppara = Dbody*V(1,ii);% Parasite power [W]
    Drotor = Dbody/Nr;% drag per rotor [N]
    AoA = atan(Drotor/Th);% disk angle of attack [rad]
    degree = rad2deg(AoA);% AoA [deg]
    Ttotal = sqrt(Drotor^2+Th^2);% Total thrust of each rotor [N]
    mu = V(1,ii)*cos(AoA)/(Vtip);
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
    vi = lambda*Vtip-V(1,ii)*sin(AoA);% induced velocity in forward flight
    Pi = Ttotal*vi;% induced power without induced loss
    Profile = 100;
    Paero = kinduced*Pi*Nr+Ppara+Profile;% no profile power because of the momentum theory based analysis
    Paerorotor = Paero/Nr;% shaft power of each rotor 
%     mskin = 1.600*8*Afuselage;% Fuselage mass
    rch = Pi/Ph;% ratio of forward flight induced power to hover 
%     Q = Paerorotor/omega;% Torque of each rotor [Nm]
%     mmotor = m1+m2*Q;% motor mass [kg]

    % Find total mass required
%     mtotal = mp+ma+mskin+mbattery+mmotor*Nr+msupport+mblade+0.1*m;% total mass [kg]
    counter = counter+1;

    time(1,ii) = Eb/Paero*3600;% mission endurance [s]
    hours(1,ii) = time(1,ii)/3600;% mission endurance [h]
    range(1,ii) = time(1,ii)*V(1,ii)/1000;% mission range [km]
    rm(1,ii) = range(1,ii)/2;% mission radius [km]
    Ps(1,ii) = Paero;
    Pis(1,ii) = kinduced*Pi*Nr;
    Pps(1,ii) = Ppara;
    vis(1,ii) = vi;

%     Qs(1,ii) = Q;
    alpha(1,ii) = degree;
%     mba(1,ii) = mbattery;
    Thrust(1,ii) = Ttotal*Nr;
end
%%
figure(1)
hold on
si = 2;
plot(V,Ps/1000,'linewidth',si)
plot(V,Pis/1000,'linewidth',si)
plot(V,Pps/1000,'linewidth',si)
h = gca;
grid off;
fs = 12;
% set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% title('ground config');
xlabel('airspeed V [m/s]','FontSize',fs);
ylabel('Power [kW] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);
legend('Total','Induced','Body drag');

a1 = 0;
a2 = 18;
b1 = 0;
b2 = 5.5;
xlim([a1 a2]);
ylim([b1 b2]);

FS = 12;
xt = a1:2:a2 ;
set(h,'XTick',xt) ;
set(h,'XTickLabel',xt,'fontsize',FS) ;
yt = b1:0.5:b2 ; 
set(h,'YTick',yt) ;
set(h,'YTickLabel',yt,'fontsize',FS) ;


%%
figure(2)
hold on
si = 2;
subplot(2,1,1)
plot(V,hours,'linewidth',si)
h = gca;
grid off;
fs = 12;

ylabel('Endurance [hours] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

subplot(2,1,2)
plot(V,range,'linewidth',si)
h = gca;
grid off;
fs = 12;

xlabel('Airspeed V [m/s]','FontSize',fs);
ylabel('Range [km] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

a1 = 0;
a2 = 18;
b1 = 0;
b2 = 5.5;
xlim([a1 a2]);
% ylim([b1 b2]);

FS = 12;
xt = a1:2:a2 ;
set(h,'XTick',xt) ;
set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:0.5:b2 ; 
% set(h,'YTick',yt) ;
% set(h,'YTickLabel',yt,'fontsize',FS) ;

%%
figure(3)
hold on
si = 2;
subplot(2,1,1)
plot(V,-alpha,'linewidth',si)
h = gca;
grid off;
fs = 12;

ylabel('Disk AoA [deg] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

subplot(2,1,2)
plot(V,Thrust,'linewidth',si)
h = gca;
grid off;
fs = 12;

xlabel('Airspeed V [m/s]','FontSize',fs);
ylabel('Net thrust [N] ','FontSize',fs);
set(h,'FontName','Times New Roman','linewidth',si);
set(h,'FontName','Times New Roman','linewidth',si);

a1 = 0;
a2 = 18;
b1 = 450;
b2 = 550;
xlim([a1 a2]);
ylim([b1 b2]);

FS = 12;
xt = a1:2:a2 ;
set(h,'XTick',xt) ;
set(h,'XTickLabel',xt,'fontsize',FS) ;
yt = b1:50:b2 ; 
set(h,'YTick',yt) ;
set(h,'YTickLabel',yt,'fontsize',FS) ;