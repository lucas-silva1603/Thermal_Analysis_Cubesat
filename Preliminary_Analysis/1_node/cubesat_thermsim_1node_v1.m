% Thermal simulation code for 1U CubeSat based on "Preliminary Thermal
% Analysis of Small Satellites", by Casper Versteeg and David L. Cotten
% Single node analysis

clear
clc


% Variable setup

% Inputs -----------------------------------------------------------
    
    % Orbital data
    h = 300;  % Altitude (km)
    orbit_count = 3;
    T0 = 273; % Initial temperature

    % Satellite data
    mass = 1; % Satellite's mass (kg)
    cp = 896; % Specific heat
    s_area = 0.06; % Surface area of the satellite

    % Simulation parameters
    timestep = 60; % Seconds


    %alpha =; % Absorptivity   (these 2 can be left as variables)
    %epsilon =; % Emissivity

    % Constants
    sbc = 0.000000567; % Stephan-Boltzmann constant
    G = 0.000000000066743;  % Gravitational constant
    M = 5.9722e+24; % Earth's mass
    R = 6378137; % Earth's radius (m)

% ------------------------------------------------------------------

% Beta angle
    % Can be left as variable

beta = 30;

% Critical beta angle (for the specififed altitude)

h = h*1000;
beta_cr = asind(R/(R+h));

% ------------------------------------------------------------------

% Orbital period

tau = 2*pi*sqrt(((h+R)^3)/(G*M));

% Eclipse fraction

if abs(beta) < beta_cr
    fe = (1/180)*acosd((sqrt(h^2+2*R*h)/(R+h)*cosd(beta)));
else 
    fe = 0;
end


% Albedo 

if beta<30
    a = 0.14;
else
    a = 0.19;
end

% Infrared heat (q_ir)

if beta<30
    q_ir = 228;
else
    q_ir = 218;
end

% ------------------------------------------------------------------

% Lumped mass 1 node aproximation (satellite as a sphere)

s_radius = sqrt(s_area/(4*pi)); % Sphere radius

disc_area = pi*(s_radius)^2; % Irradiated area of the sphere is a disc
A_ir = disc_area; % Infrared area
A_sun = disc_area; % Sun and albedo area
A_s = s_area; % Radaiting area

albedo = a;
q_sun = 1361; % Solar radiation (W/m^2)
q_int = 10; % Internal power (W)

% ------------------------------------------------------------------
% Simulation loop (for variable emittance)

alpha_list=[];
time_list=[];
Temp_list=[];

for alpha = 0:0.01:1
    epsilon = alpha;
    T=T0;

    time_list(end+1, 1) = 0;    
    alpha_list(end+1, 1) = alpha;
    Temp_list(end+1, 1) = T;


    for t = timestep:timestep:tau*orbit_count

        % Solar incidence multiplier (turns the sun exposure on and off)
%{
        if (tau/2)*(1+fe)<t && t<(tau/2)*(1-fe)
            s=1;
        else
            s=0;
        end
%}
        Q = q_ir*A_ir + (1+albedo)*q_sun*A_sun*s(t,tau,fe)*alpha + q_int - A_s*sbc*epsilon*T^4;
        T = T + (timestep/cp*mass)*Q;

        time_list(end+1, 1) = t;    
        alpha_list(end+1, 1) = alpha;
        Temp_list(end+1, 1) = T;

    end
end
dados = cat(2, time_list, alpha_list, Temp_list);
dados2=single(dados);
%clear time_list alpha_list Temp_list dados



%{
% Simulation loop (for variable initial temperature)

initial_temp_list=[];
time_list=[];
Temp_list=[];

for T0 = 0:1:373

    epsilon = alpha;
    T = T0;

    time_list = [time_list, 0];    
    initial_temp_list =[initial_temp_list, T0];
    Temp_list = [Temp_list, T];

    for t = timestep:timestep:tau*orbit_count

        % Solar incidence multiplier (turns the sun exposure on and off)

        if (tau/2)*(1+fe)<t && t<(tau/2)*(1-fe)
            s=1;
        else
            s=0;
        end

        Q = q_ir*A_ir + (1+albedo)*q_sun*A_sun*s*alpha + q_int - A_s*sbc*epsilon*T^4;
        T = T + (timestep/cp*mass)*Q;

        time_list = [time_list, t];
        initial_temp_list =[initial_temp_list, alpha];
        Temp_list = [Temp_list, T];


    end
end


% Simulation loop (for variable beta angle)

beta_list=[];
time_list=[];
Temp_list=[];

for beta = -90:1:90

    epsilon = alpha;
    T = T0;

    time_list = [time_list, 0];    
    beta_list =[beta_list, beta];
    Temp_list = [Temp_list, T];

    for t = timestep:timestep:tau*orbit_count

        % Solar incidence multiplier (turns the sun exposure on and off)

        if (tau/2)*(1+fe)<t && t<(tau/2)*(1-fe)
            s=1;
        else
            s=0;
        end

        Q = q_ir*A_ir + (1+albedo)*q_sun*A_sun*s*alpha + q_int - A_s*sbc*epsilon*T^4;
        T = T + (timestep/cp*mass)*Q;

        time_list = [time_list, t];  
        beta_list =[beta_list, beta];        
        Temp_list = [Temp_list, T];


    end
end




%}

