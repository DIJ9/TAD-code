%%%% Titan Aerial Daughtercraft (TAD) %%%%%%%%%%%%%%%%
%%%% Momentum theory based parametric analysis %%%%%%%
%%%% Lander mission design %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 06/18/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Daiju uehara, Larry Matthies, Phil Tokumaru %%%%% 

clc 
% Define planet paremeters on Titan
a = 190;% speed of sound [m/s]
rho = 5.34;% air density [kg/m^3]
g = 1.352;% gravity [m/s^2]
nu = 1.2e-6;% kinematic viscosity [m^2/s]

%%
% Define vehicle parameters
m = 6;% vehicle mass [kg]
N = m*g;% vehicle weight [N]
Nr = 4;% number of rotor, assuming a quadcopter config

% Hover performance assumption
FM = 0.75;% rotor figure of merit is given parameter in this analysis 
Mtip = 0.2;% Tip mach number
DL = 50/g;% disk loading [kg/m^2] typical for toy-scale helicopter?

% Parameters of each rotor 
Vtip = Mtip*a;% Tip speed [m/s]
A = m/Nr/DL;% disk area
R = sqrt(A/pi);% rotor radius [m]
omega = Vtip/R;% angular velocity [rad/s]
RPM = 2*pi/omega*60;% angular velocity [RPM]

% Ideal hover power 
Th = N/4;%thrust of each rotor for hover
vh = sqrt(DL*g/2/rho);% ideal induced velocity in hover
Ph = Th*vh;% ideal induced power in hover

% Payload fuselage design
it = 0.06;% insulation thickness [m]
pl = 0.17;% payload length [m]
pw = 0.12;% payload width [m]
ph = 0.07;% payload height [m]
Ap = 0.012;% payload frontal area [m^2] or pw x ph 
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
% Forward flight analysis
V = 0.1:0.1:10;% forward speed flight range
alpha = [0:1:20]*pi/180;% disk AoA range
% Optimize disk angle of attack till maximum mission radius is found
% AoAop = zeros(1,length(V));
mbattery = 0.1:0.05:1;
Pstorage = zeros(length(V),length(mbattery));
mstorage = zeros(length(V),length(mbattery));
rmstorage = zeros(length(V),length(mbattery));

for ii = 1:length(V)
    for jj = 1:length(mbattery)
        evolume = mbattery(1,jj)/de;% energy storage volume [m^3]
        wbvolume = evolume+va;% warm box volume [m^3]
        wblength = (wbvolume)^(1/3);% warm box side length [m]
        Awb = (wblength+it)^2;% warm box frontal area with insulation [m^2]
        Afuselage = Awb+Ap;% payload + warm box [m^2]
        lre = Afuselage/pw;% characteristic length for Re [m]
        Re = V(1,ii)*lre/nu;% Reynolds number
        if Re < 2e5
            Cdbody = 0.47;% drag coefficient of vehicle body
        else
            Cdbody = 0.18;
        end
        Dbody = 1/2*Cdmulti*Cdbody*V(1,ii)^2*rho*Afuselage;% parasite drag [N]
        Ppara = Dbody*V(1,ii);% Parasite power [W]
        Drotor = Dbody/4;% drag per rotor [N]
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
        vi = lambda*Vtip;% induced velocity in forward flight
        Pi = Ttotal*(V(1,ii)*sin(AoA)+vi);% induced power without induced loss
        Paero = kinduced*Pi*4+Ppara;% no profile power because of the momentum theory based analysis
        Paerorotor = Paero/4;% shaft power of each rotor 
        mskin = 1.600*8*Afuselage;% Fuselage mass
        rch = Pi/Ph;% ratio of forward flight induced power to hover 
        Q = Paero/omega;% Torque of each rotor [Nm]
        mmotor = m2+m1*Q;% motor mass [kg]
        
        % Find total mass required
        mtotal = mp+ma+mskin+mbattery(1,jj)+mmotor*Nr+msupport+mblade+0.05*m;% total mass [kg]
        
        % Check if it is make sense 
        diff = abs(m-mtotal);
        
        Eb = Ebconst*mbattery(1,jj);% battery energy [Wh]
        time = Eb/Paero*3600;% mission endurance [s]
        hours = time/3600;% mission endurance [h]
        range = time*V(1,ii)/1000;% mission range [km]
        rm = range/2;% mission radius [km]
        mstorage(ii,jj) = mtotal;
        Pstorage(ii,jj) = Paero;
        rmstorage(ii,jj) = rm;
        

    end    
end
[X,Y] = meshgrid(mbattery,V);
surf(X,Y,rmstorage)
