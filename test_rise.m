clear all
close all
clc

rho_b = 1.225;
sigma       = 71.99e-3;          % Surface tension for air in water at 298 K [N/m]
rho_f       = 1000;              % Density of water [kg/m^3]
mu_f        = 1e-3;
mu_b        = 1.81e-5;
Ut = [];
RiseParamsStokes    = {'Stokes' num2str(rho_b) num2str(rho_f) num2str(mu_f) num2str(sigma)}';
RiseParamsKaramanev = {'Karamanev' num2str(rho_b) num2str(rho_f) num2str(mu_f) num2str(sigma)}';
RiseParamsDetschPure= {'Detsch' 'PureWater' num2str(rho_f) num2str(mu_b) num2str(sigma)}';
r = 1e-6*2;%linspace(1e-6,1e-3,100);
for i=1%:100
     UtStokes(i) = rise_2( r(i), RiseParamsStokes)
     UtDetsch(i) = rise_2( r(i), RiseParamsDetschPure);
     UtKaramanev(i) = rise_2( r(i), RiseParamsKaramanev);

end

transit_separator = [70.8, 35.4,23.6,17.7,14.2,11.8,10.1,8.9,7.9,7.1,6.4,5.9];

riseheightStokes = UtStokes.*transit_separator;
riseheightDetsch = UtDetsch.*transit_separator;
riseheightKaramanev = UtKaramanev.*transit_separator;