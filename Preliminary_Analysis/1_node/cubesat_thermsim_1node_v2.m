% Thermal simulation code for 1U CubeSat based on "Preliminary Thermal
% Analysis of Small Satellites", by Casper Versteeg and David L. Cotten
% Single node analysis

clear
clearvars
%clc


% Variable setup

% Inputs -----------------------------------------------------------
    
    % Orbital data
    h = 450;  % Altitude (km)
    orbit_count = 20;
    T0 = 273; % Initial temperature

    beta = 45; % Beta angle (orbit inclination from equator)


    % Satellite data
    mass = 1.33; % Satellite's mass (kg)
    cp = 896; % Specific heat
    s_area = 0.06; % Surface area of the satellite
    q_int = 5; % Internal power (W)

    % Simulation parameters
    timestep = 1; % Seconds (apenas funciona para 1 por agora) 

% ------------------------------------------------------------------

% Constants
sbc = 0.000000567; % Stephan-Boltzmann constant
G = 0.000000000066743;  % Gravitational constant
M = 5.9722e+24; % Earth's mass
R = 6378137; % Earth's radius (m)

% ------------------------------------------------------------------

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


% ------------------------------------------------------------------


% Simulation loop (for variable emittance)

%alpha_list=[];
%time_list=[];
%Temp_list=[];

alphas = linspace(0,1,100);
num_steps = (tau*orbit_count)/timestep;
time = (0:num_steps-1) * timestep;

T_over_time = zeros(floor(num_steps),floor(length(alphas)));
oc=1;




for i = 1:int64(length(alphas))
    
    alpha = alphas(i);
    epsilon = alpha;
    T=T0;

    %time_list(end+1, 1) = 0;    
    %alpha_list(end+1, 1) = alpha;
    %Temp_list(end+1, 1) = T;


    for j = 0:num_steps-1

        t = time(j+1);

        current_orbit_time = mod(t, tau);
       
        if t==0
          T_over_time(j+1,i) = T0;  
        else
            Q = q_ir*A_ir + (1+albedo)*q_sun*A_sun*s(current_orbit_time,tau,fe)*alpha + q_int - A_s*sbc*epsilon*T^4;
            T = T + (timestep/cp*mass)*Q;
            T_over_time(j+1,i) = T;
        end
               

        

    end
end

max(T_over_time, [], 'all')

time = (0:num_steps-1) * timestep;
figure;
surf(alphas, time, T_over_time);
shading interp
colormap jet
xlabel('Alpha');
ylabel('Time (s)');
zlabel('Temperature (K)');
title('Satellite Temperature Over Time and Absorptivity');
colorbar;
hold

%exportgraphics(gcf,'Graph1.png','BackgroundColor','black','Resolution',1000);

figure;
plot(time,T_over_time(:,11))
hold

%exportgraphics(gcf,'Graph2.png', 'BackgroundColor','black','Resolution',1000);


function s_v = s(t, tau, fe)

    if (((tau/2)*(1+fe)))>t

        if t>(((tau/2)*(1-fe)))
            s_v=0;
        else
            s_v=1;
        end
    else
        s_v = 1;
    end
   
end