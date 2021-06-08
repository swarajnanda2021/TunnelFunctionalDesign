function [ Ut ] = rise_2( r, params)
%   Rise: Based on tunnel working conditions, displaces the specified bubble
%   population by the terminal velocity estimated from a chosen model. 
   
%   There are 3 models put in so far and their input params are described:
%       1) The terminal velocity correlation for air-bubble in contaminated
%       fluid by D. G. Karamanev.
%           Karamanev, D. G. "Equations for calculation of the terminal 
%           velocity and drag coefficient of solid spheres and gas 
%           bubbles." Chemical engineering communications 147.1 (1996): 75-84.
%       
%       Input params: 
%       bubbles, dt and params (rho_p, rho_f, mu, sigma)
%        


g = 9.81;                                                           % Gravity in m/s^2
if strcmp(params(1),'Karamanev') == 1
        % Initialising parameters from params array
        
        rho_p = str2double(cellstr(params(2)));
        rho_f = str2double(cellstr(params(3)));
        mu    = str2double(cellstr(params(4)));
        sigma = str2double(cellstr(params(5)));
        % Part 1) Calculate Archimedes number first
        Ar = (2*r)^3 * g *abs(rho_p - rho_f) * rho_f / mu^2;
        if Ar > 13000
            Cd = 0.95;
        else
            Cd =  (432/Ar) * (1 + (0.047 * (Ar^(2/3)))) + (0.517/(1 + (154 * (Ar^(-1/3))))) ;
        end
        
        % Part 2) Calculate d_e/d_h
        Eo = g * (rho_f - rho_p) * (2*r)^2 / sigma;
        if Eo > 40
            debydh = 0.62;
        else
            debydh = (1 + 0.163 * (Eo^0.757))^(-1/3);
        end
    
        % Part 3) Calculate the terminal velocity from the expression of D. G.
        % Karamanev
        V = (4/3) * pi * r^3;                                   % Volume of the bubble in consideration
        Ut = 40.3 * (debydh) * (V^(1/3) / Cd)^0.5;
        
end

%
%       
%       2) The stokes law based terminal velocity, simulating the rise of a
%       solid sphere with the idea that since the water is impure,
%       impurities destroy the surface circulation of fluid and make it
%       behave as a solid sphere.
%       
%       Input params: 
%
%       bubbles, dt and params (rho_p, rho_f, mu, sigma)
%        

if strcmp(params(1),'Stokes') == 1
        
        % Initialising parameters from params array
        
        rho_p = str2double(cellstr(params(2)));
        rho_f = str2double(cellstr(params(3)));
        mu    = str2double(cellstr(params(4)));
        
        % Terminal velocity expression by Stokes
        Ut = (g * (2*r)^2 * abs(rho_p - rho_f)) / (18 * mu);
end




%       3) The terminal velocity curve fit for air-bubbles in seawater, by
%       Detsch.
%           Detsch, Richard M. "Small air bubbles in reagent grade water 
%           and seawater: 1. Rise velocities of 20?to 1000??m?diameter 
%           bubbles." Journal of Geophysical Research: Oceans 96.C5 (1991): 8901-8906.
%       
%       Input params: 
%       
%       bubbles and dt
%       

if strcmp(params(1),'Detsch') == 1
        if strcmp(params(2),'PureWater') == 1
            % Initialising parameters from params array
        
            d = 2*r*1e6;
            
            % Terminal velocity expression by Stokes
            Ut = abs(((6.82e-2) + (3.82e-3 * d) + (1.83e-5 * d^2))*0.01);
    
        end
        if strcmp(params(2),'SeaWater') == 1
            % Initialising parameters from params array
        
            d = 2*r*1e6;
            
            % Terminal velocity expression by Stokes
            Ut = abs(((-4.17e-1) + (1.12e-2 * d) + (1.42e-6 * d^2))*0.01);
        
        end

end


end

