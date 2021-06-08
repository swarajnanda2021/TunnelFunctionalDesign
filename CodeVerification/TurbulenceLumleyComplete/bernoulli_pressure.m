function [Pressure,  p_inlet ] = bernoulli_pressure( Q, Cross_Section_area, P_free_surface, Static_pressure )
% ORIGINALLY: [ P_liquid, P_Stat_dyn, Velocity, p_inlet ] = bernoulli( Q, Cross_Section_area, P_free_surface, Static_pressure )
% Calculation of pressure from Bernoulli

% Part 1) Calculation of the velocity profile in the tunnel circuit
Velocity = (Cross_Section_area).^(-1) * Q;                      % From continuity, Area*velocity = Q = constant

% Part 2) Calculation of the Total pressure within the tunnel
Pres_total = (0.5*1000.*(Velocity.^2)) + Static_pressure;       % Adding Static and Dynamic pressure
P_Stat_dyn = (Pres_total.*(10^(-5)))+P_free_surface;

% Part 3) Calculation of Pressure from Bernoulli
p_inlet = max(P_Stat_dyn) - min(P_Stat_dyn);                    % IMPORTANT: Pump outlet pressure is by taking pressure difference between the static and dynamic pressure at the end of the circuit and that at the beginning
Constant = p_inlet + P_Stat_dyn(1);                             % Find constant for Bernoulli
P_liquid = Constant - (P_Stat_dyn);                             % Calculate liquid pressure using Bernoulli

Pressure = P_liquid + Static_pressure + P_free_surface;

end

