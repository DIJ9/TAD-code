%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Balloon mission design %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 06/20/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
m = 4:0.5:20;% vehicle mass [kg]
% Axial flight state (return to mothership)
Vc = 8;% climb rate [m/s] 
V = 4:1:14;% forward speed flight range
time = zeros(length(m),length(V));
hours = zeros(length(m),length(V));
range = zeros(length(m),length(V));
rm = zeros(length(m),length(V));
Ps = zeros(length(m),length(V));
marg = zeros(length(m),length(V));
mba = zeros(length(m),length(V));
Thrust = zeros(length(m),length(V));
Pis = zeros(length(m),length(V));
Pps = zeros(length(m),length(V));
vis = zeros(length(m),length(V));
Pcas = zeros(length(m),length(V));
alpha = zeros(length(m),length(V));
Phover = zeros(length(m),length(V));

for jj = 1:length(m)
N = m(1,jj)*g;% vehicle weight [N]
Nr = 4;% number of rotor, assuming a quadcopter config

% Hover performance assumption
FM = 0.75;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.25;% Tip mach number
DL = 50/g;% disk loading of each rotor [kg/m^2] 

% Efficiency coefficients
etam = 0.85;% motor efficiency 
etac = 0.95;% control efficiency 
kinduced = 1.15;% induced power factor 
etafw = FM*etam*etac;% total forward flight efficiency

% Parameters of each rotor 
Vtip = Mtip*a;% Tip speed [m/s]
A = m(1,jj)/Nr/DL;% disk area
R = sqrt(A/pi);% rotor radius [m]
omega = Vtip/R;% angular velocity [rad/s]
RPM = 2*pi/omega*60;% angular velocity [RPM]

% Ideal hover power 
Th = N/4;%thrust of each rotor for hover
vh = sqrt(DL*g/2/rho);% ideal induced velocity in hover
Ph = Th*vh;% ideal induced power in hover per rotor 
Phtotal = Ph*4;% Total induced power in hover

% Actual hover power 
Pha = Ph/FM/etam/etac;% actual hover power consumption [W]
Phatotal = 4*Pha;% Total actual power in hover

% Sampling state analysis
tsample = 5*60;% sampling time [s]
Psample = 20;% sampling power [W]
Esample = (Psample*tsample+Phatotal*tsample)/3600;% sampling energy [Wh]

% Balloon parameter
ab = 10000;% balloon altitude [m]
Vbal = 2;% wind velocity [m/s]

% Payload/fuselage design
df = 1200;% fuselage density [kg/m^3]
rf = (m(1,jj)/df*3/4*1/pi)^(1/3);
Af = pi*rf^2;% fuselage frontal area [m^2]

mp = 2;% payload mass [kg]
ma = 0.5;% avionics mass [kg]

Cdmulti = 1.5;% drag multiplier

% Blade parameters
dblade = 0.15;% blade mass density based on disk area [kg/m^2]
mblade = dblade*A*4;

% Motor constant 
m1 = 0.003;% motor constant [kg]
m2 = 0.322;% motor constant [kg/N-m]

% Battery specific energy
Ebconst = 100;% Wh/kg

%%
ctime = ab/Vc;% time to climb [s]
dtime = 5*60;% time to dock [s]
Edock = dtime*Phatotal/3600;
Vd = -8;% descent flight
detime = ab/abs(Vd);% time to descend [s]

% balloon travel distance during descend, climb, sampling
dbal = 2*(ctime+detime+tsample);% [m]
    
    for ii = 1:length(V)
        % Calculate power in forward flight 
        Re = V(1,ii)*2*rf/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end    
        Dbody = 1/2*Cdmulti*Cdbody*V(1,ii)^2*rho*Af;% parasite drag [N]
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
        Profile = 20;% Assumed profile power
        Paero = Pi*Nr/etafw+Ppara+Profile;
        Paerorotor = Paero/Nr;% shaft power of each rotor 
        tlevel = dbal/(V(1,ii)-Vbal);% level flight duration to catch balloon
        Ef = Paero*tlevel/3600;% Total forward flight energy

        % Calculate power in climb state
        vclimb = vh*(-Vc/(2*vh)+sqrt(1+(Vc/(2*vh))^2));% induced velocity in climb [m/s]
        Pclimb = Ph*(Vc/(2*vh)+sqrt(1+(Vc/(2*vh))^2));% induced power in climb [W]
        Tclimb = 2*rho*A*(Vc+vclimb)*vclimb;% thrust
        Pca = Pclimb/FM/etac/etam;% actual power [W]
        Pcatotal = Pca*4;% Total climb power 
        Ec = Pcatotal*ctime/3600;% Total climb energy 

        % Calculate power in descent state (windmill state Vd<-2vh)
        k = 1.2;
        k1 = -1.125;
        k2 = -1.372;
        k3 = -1.718;
        k4 = -0.655;
        
        if Vd/vh <= 0 && -2 <= Vd/vh
            vdescent = vh*(k+k1*(Vd/vh)+k2*(Vd/vh)^2+k3*(Vd/vh)^3+k4*(Vd/vh)^4);
        else
            vdescent = vh*(-Vd/(2*vh)-sqrt((Vd/(2*vh))^2-1));
        end
        Pdescent = Ph*(Vd/(2*vh)-sqrt((Vd/(2*vh))^2-1));
        Tdescent = -2*rho*A*(Vd+vdescent)*vdescent;

        % Find maximum power among hover, climb, and forward flight
        [Pmotor,I] = max([Pha Paerorotor Pca]);

        % Compute motor mass
        Q = Pmotor/omega;% Torque of each rotor [Nm]
        mmotor = m1+m2*Q;% motor mass [kg]

        % Compute battery mass 
        mbattery = m(1,jj)-mp-ma-0.2*m(1,jj)-mmotor*Nr-mblade;

        % Compute energy
        Eb = Ebconst*mbattery;% available battery energy [Wh]
        Eba = Ec+Ef+Edock+Esample;% Total estimated energy
        margin = (Eb-Eba)/Ebconst/m(1,jj)*100;

        marg(jj,ii) = margin;
        time(jj,ii) = Eb/Paero*3600;% mission endurance [s]
        hours(jj,ii) = time(jj,ii)/3600;% mission endurance [h]
        range(jj,ii) = time(jj,ii)*V(1,ii)/1000;% mission range [km]
        rm(jj,ii) = range(jj,ii)/2;% mission radius [km]
        Ps(jj,ii) = Paero;
        Pcas(jj,ii) = Pcatotal; 
        alpha(jj,ii) = degree;
        mba(jj,ii) = mbattery;
        Thrust(jj,ii) = Ttotal;
        Pis(jj,ii) = Pi*Nr/etafw;
        Pps(jj,ii) = Ppara;
        vis(jj,ii) = vi;
        Phover(jj,ii) = Phatotal;


    end
end
%%
[X,Y] = meshgrid(V,m);
figure(1)
contourf(X,Y,marg,50,'edgecolor','none');
colorbar
%%
figure(2)
hold on
si = 2;
h1 = plot(V,Ps,'linewidth',si);
plot(V,Pis,'linewidth',si)
plot(V,Pps,'linewidth',si)
hh5 = plot(V(1,Imin),Ps(1,Imin),'ok','linewidth',si);

%%
figure(3)
plot(m,marg)
%%
% figure(1)
% hold on
% si = 2;
% plot(V,Ps,'linewidth',si)
% plot(V,Pis,'linewidth',si)
% plot(V,Pps,'linewidth',si)
% h = gca;
% grid off;
% fs = 12;
% % set(h,'GridAlpha',0.1,'GridLineStyle','--','GridColor',[0.5 0.5 0.5]);
% % title('ground config');
% xlabel('airspeed V [m/s]','FontSize',fs);
% ylabel('Power [W] ','FontSize',fs);
% set(h,'FontName','Times New Roman','linewidth',si);
% set(h,'FontName','Times New Roman','linewidth',si);
% legend('Total','Induced','Body drag');
% 
% a1 = 0;
% a2 = V(1,end);
% b1 = 0;
% b2 = 80;
% xlim([a1 a2]);
% ylim([b1 b2]);
% 
% FS = 12;
% xt = a1:2:a2 ;
% set(h,'XTick',xt) ;
% set(h,'XTickLabel',xt,'fontsize',FS) ;
% yt = b1:5:b2 ; 
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


