    clear all
close all
clc
%% Initialisation
% Initialising parameters from params array
        
rho_p = 1.225;
rho_f = 1000;
mu    = 1e-3;
r = linspace(1e-6,1e-3,1000).*0.5;
g = 9.81;
sigma = 71e-3;
%% Terminal velocity expression by Stokes
for i=1:length(r)
    Ut_stokes(i) = (g * (2*r(i))^2 * abs(rho_p - rho_f)) / (18 * mu); 
end

%% Terminal velocity expression by Karamanev

% Part 1) Calculate Archimedes number first
for i=1:length(r)
    Ar = (2*r(i))^3 * g *abs(rho_p - rho_f) * rho_f / mu^2;
    if Ar > 13000
        Cd = 0.95;
    else
        Cd =  (432/Ar) * (1 + (0.047 * (Ar^(2/3)))) + (0.517/(1 + (154 * (Ar^(-1/3))))) ;
    end

    % Part 2) Calculate d_e/d_h
    Eo = g * (rho_f - rho_p) * (2*r(i))^2 / sigma;
    if Eo > 40
        debydh = 0.62;
    else
        debydh = (1 + 0.163 * (Eo^0.757))^(-1/3);
    end

    % Part 3) Calculate the terminal velocity from the expression of D. G.
    % Karamanev
    V = (4/3) * pi * r(i)^3;                                   % Volume of the bubble in consideration
    Ut_karamanev(i) = 40.3 * (debydh) * (V^(1/3) / Cd)^0.5;
end

%% Terminal velocity expression by Detsch for seawater and pure water

for i=1:length(r)
    
    d = 2*r(i)*1e6;
            
    % Terminal velocity expression by Stokes
    Ut_detsch_pure(i) = abs(((6.82e-2) + (3.82e-3 * d) + (1.83e-5 * d^2))*0.01);
%     Ut_detsch_sea(i)  = abs(((-4.17e-1) + (1.12e-2 * d) + (1.42e-6 * d^2))*0.01);

end





%% Plotting

figure(1)
loglog(r.*2,Ut_stokes)
hold all
loglog (r.*2,Ut_karamanev)
loglog (r.*2,Ut_detsch_pure)
xlabel('Bubble diameter (m)')
ylabel('Terminal rise velocity (m/sec)')
legend('Stokes','Karamanev (\sigma:71 mN/m)','Detsch')
title('Terminal velocities for: \mu_{bubble}=1.81e-5 kg m/sec, \rho_{fluid} = 1000 kg/m^3, \rho_{bubble} = 1.225 kg/m^3')
% loglog (r.*2,Ut_detsch_sea)





