%% Script for initialising calculations for the TUD Multiphase flow tunnel
clear all
close all
clc
%% Table Of Contents
%
% 1. Initialising base geometrical parameters for the various tunnel configurations
% 1.1. Setting base params for the tunnel
% 1.2. Setting params for large test-section configuration
% 1.3. Setting params for small test-section configuration
%
% 2. Initialising parameters for tunnel simulations
% 2.1. Setting params for bubble-dissolution simulation
% 2.2. Setting params for bubble rise simulation
% 2.3. Setting params for turbulence reduction simulation
% 2.4. Setting params for pump power requirement calculations
% 
% 3. Initialisation of calculations
% 3.1. Calculation of derived parameters from the prescribed parameters
% 3.2. For loop for the TunnelMainCode.m code over volumetric flowrate of
% pump
% 
% 4. Compilation of pump and system curves
%
% 5. Calculation and plotting of loading parameters for the tunnel limbs

%% 1 Initialising base geometrical parameters for various tunnel configurations

% 1.1. Setting base params for the tunnel
% Static head extremas, pump size including diffuser, main diffuser params
Hmin            = 0.5;                                                      % Minimum static-head acting on the centre streamline
Dpump           = 0.55;                                                    % Diameter of the pump
ldr             = 5.71;                                                        % Length of downcomer/riser
r               = 0.1;                                                      % Rounding of all bends in the circuit
lsa             = 2.5;                                                        % Lenght of shape-adapter circuit 
l2              = 1.24;                                                        % Length of small riser

% Second pipe vertical has 4.26 

lmd             = 1.5;                                                        % Length of the main diffuser
Wexpansion      = 1;                                                        % Outlet width of the main diffuser
rho             = 1000;                                                     % Fluid density in S.I. units

% Separator internals
nBanks = 3;
nPlates = 8;
nPlatesFunctional = 5; % How many Plates are actually in use? Variable bank separator

% Gas removal
L_KnifeCasing   = 0.2  ; % Length of the casing the holds the exhaust
alpha_opening   = 0  ; % Opening angle of the knife (in degrees)
l_vent          = 0.15  ; % length of the venting knife


% Set configuration of the tunnel 
Config = 2;                                                                 
if Config == 1
% 1.2. Setting params for large test-section configuration
    % Required params
    TotalComponents     = 17;                                               % Total number of components
    NameMain            = {'Riser','Riser-adapter','Bend 1','Separator','Bend 2','Small riser','Bend 3','Honeycomb','Settling chamber','Contraction','Large test-section',...
        'Main diffuser','Exhaust-knife','Bend 4','Downcoming-adapter','Downcomer','180 Bend'};
    
    DiffuserVanesCompNo = [2 12];
    
    Lcontraction2       = 0.75;
    Ldynamometer        = 0.5;
    Lmaindiffuser2      = 5 - 2.13 - Ldynamometer - Lcontraction2;
    

elseif Config == 2
% 1.3. Setting params for small test-section configuration   
    % Required params
    TotalComponents     = 19;                                               % Total number of components
    NameMain            = {'Riser','Riser-adapter','Bend 1','Separator','Bend 2','Small riser','Bend 3','Honeycomb','Settling chamber','Contraction','Dynamometer','Small Contraction',...
        'Cav section','Small main diff.','Main diffuser','Exhaust-knife','Bend 4','Downcoming-adapter','Downcomer'};
    
    DiffuserVanesCompNo = [2 14 15];
    
    Lcontraction2       = 0.75;
    Ldynamometer        = 0.5;
    Lmaindiffuser2      = 5 - 2.13 - Ldynamometer - Lcontraction2;

end

%% 2. Initialising parameters for tunnel simulations

% 2.1. Setting params for bubble-dissolution/rise simulation

% Set Pressure
SetP                    = 0.05*1e5;                                         % Pressure set at Pstat=rho*g*Hmin (EXTREME VALUE)

% 2.2. Setting params for bubble rise simulation (Remember, RiseStart and
% RiseEnd have been set according to Configuration == 1)
% For modification to put separator first 
r_separation            = 0.1*1e-3;                                         % Size of bubble intended for macrobubble separation
RiseStart               = 4;                                                % Component number of macrobubble separator 
RiseEnd                 = RiseStart + 1;                                    % Number following the macrobubble separator

if Config == 1
    Test_Section            =   11;
else
    Test_Section            =   13;
end
    

% 2.3. Setting params for turbulence reduction simulation
Lhoneycomb              = 0.75;                                             % Length of honeycomb
Lsettling               = 0.65;                                             % Length of Settling chamber
MeshWidthHoneycomb      = 0.009;                                            % Width of the honeycomb element
CellLengthHoneycomb     = 0.4;                                              % Length of a honeycomb cell
LInletHoneycomb         = 0.2;                                              % Length scale at the inlet to the honeycomb
Lcontraction            = 1.5;                                             % Length of the circle to square contraction

% 2.4. Setting params for pump power requirement calculations
QPUMP       = 1.08;%linspace(0.108*0.5,1.08,25);                                  % Volumetric flowrate of the pump
RoughMain   = 0.01*1e-3;                                                    % Roughness of the walls, modelled as uniform roughness

%% 3. Initialisation of calculations

l_gap = 1.4 - 2*Dpump; % Increase this gap

% 3.1. Calculation of derived parameters from the prescribed parameters
Ltop                      = 0.5 + r + lmd + 5 + Lcontraction + Lsettling + Lhoneycomb + r + 0.5 + L_KnifeCasing;
Lseparator                = Ltop - Dpump - l_gap - Dpump - r - r - 0.5;

L_limb1                   = 0.5 + r + l2 + r + 0.5 + r + lsa + ldr;

ladj                      = L_limb1 - 0.5 - r - lsa - ldr    ;                                                    % Length of adjuster


% 3.2.      Execute TunnelMainC7.m

for i=1:length(QPUMP)
      [TunnelParams,LossMain,TunnelVolume,NPSH_a(i),TIEvolMain,FinalHeight(i,:)] = ...
        TunnelMainC7(Config,Hmin,Dpump,Lhoneycomb,Lsettling,Lcontraction,lmd,Wexpansion...
        ,r,MeshWidthHoneycomb,CellLengthHoneycomb,LInletHoneycomb,QPUMP(i),RoughMain,SetP,...
        Lseparator,RiseStart,RiseEnd,NameMain,DiffuserVanesCompNo,ldr,lsa,l2,ladj,nBanks,nPlates,r_separation,...
        Lcontraction2,Ldynamometer,Lmaindiffuser2,l_gap,Test_Section,nPlatesFunctional,L_KnifeCasing,alpha_opening,l_vent);
    
    System_Head(i)  = sum(LossMain)/(9.81*1000);
    if Config == 2
        Uts = QPUMP(i)/(0.3^2);
    else
        Uts = QPUMP(i) / (0.5^2);
    end
    Energy_Ratio(i) = (sum(LossMain)/(0.5*1000*Uts^2))^-1;
end


%% 4. Compilation of pump and system curves

% Load tables
load('X_RPM_Y_Flowrate_Z_TotalHead.mat');
load('X_RPM_Y_Flowrate_Z_NPSHR.mat');
load('X_RPM_Y_FlowRate_Z_Power.mat');
load('X_RPM_Y_FlowRate_Z_Eta.mat');     %Use the efficiency data only for plotting purposes

% Construct scattered interpolant for total head produced by pump, required
% suction head, required power and the pump efficiency
f_RPM = scatteredInterpolant(table2array(XRPMYFlowrateZMLC(:,2)),table2array(XRPMYFlowrateZMLC(:,3)),table2array(XRPMYFlowrateZMLC(:,1)));
% f_NPSHR     = scatteredInterpolant(table2array(XRPMYFlowrateZNPSHR(:,1)),table2array(XRPMYFlowrateZNPSHR(:,2)),table2array(XRPMYFlowrateZNPSHR(:,3)));
f_Power     = scatteredInterpolant(table2array(XRPMYFlowRateZPower(:,1)),table2array(XRPMYFlowRateZPower(:,2)),table2array(XRPMYFlowRateZPower(:,3)));
f_Eta2      = scatteredInterpolant(XRPMYFlowrateZEta1(:,2),XRPMYFlowrateZEta1(:,3),XRPMYFlowrateZEta1(:,1));

% Find required power using the RPM
RPM_range = f_RPM(QPUMP.*(3600),System_Head);
Required_Power = f_Power(RPM_range,QPUMP.*(3600));

% Plot System required head curve over pump characteristic curves (efficiency scatter curve and RPM curve)

figure(12)
subplot(3,1,1)
yyaxis left
plot(QPUMP.*(3600),System_Head)
hold on
% xlabel('Capacity (m^3/h)')
ylabel('Total Head (mlc)')
if Config < 3
    scatter3(XRPMYFlowrateZEta1(1:183,2),XRPMYFlowrateZEta1(1:183,3),XRPMYFlowrateZEta1(1:183,1),15,XRPMYFlowrateZEta1(1:183,1))
    xlim([1600 5200])
    ylim([0 3.6])
end

c = colorbar;
c.Label.String = 'Efficiency (%)';
legend('System Curve','Pump efficiency')
grid on

set(gca,'xticklabel',{[]}) 

subplot(3,1,2)
plot(QPUMP.*(3600),NPSH_a)
hold on
% xlabel('Capacity (m^3/h)')
ylabel('NPSH (mlc)')
if Config < 3
    scatter3(table2array(XRPMYFlowrateZNPSHR(:,2)),table2array(XRPMYFlowrateZNPSHR(:,3)),table2array(XRPMYFlowrateZNPSHR(:,1)),15,table2array(XRPMYFlowrateZNPSHR(:,1)))
    xlim([1600 5200])
%     ylim([0 6])
end
legend('Available (mlc)','Required (mlc)')
c = colorbar;
c.Label.String = 'Impeller RPM';
grid on
set(gca,'xticklabel',{[]}) 


subplot(3,1,3)
plot(QPUMP.*(3600),Required_Power)
hold on
xlabel('Capacity (m^3/h)')
ylabel('Power (kW)')
if Config < 3
    scatter3(table2array(XRPMYFlowRateZPower(:,2)),table2array(XRPMYFlowRateZPower(:,3)),table2array(XRPMYFlowRateZPower(:,1)),15,table2array(XRPMYFlowRateZPower(:,1)))
    xlim([1600 5200])
    ylim([0 50])
end
legend('System required Power (kW)','RPM range')
c = colorbar;
c.Label.String = 'Impeller RPM';
grid on

%% 5. Plotting the bubble rise calculations

figure(10)
semilogy(QPUMP.*(3600),FinalHeight(:,1));hold on
semilogy(QPUMP.*(3600),FinalHeight(:,2))
semilogy(QPUMP.*(3600),FinalHeight(:,3))
xlabel('Capacity (m^3/h)')
ylabel('Maximum vertical displacement [m]')
grid on
% set(gca, 'XTick', 1:3, 'XTickLabel', {'Stokes' 'Karamanev' 'Detsch (Pure)'});
legend('Stokes bubble','Karmanevs water','Pure water (Detsch)')


%% 6. Calculation and plotting of loading parameters for the tunnel limbs

% 
% 
% 
% Length_comp     = Lengths(:,1);
% Volume_comp     = TunnelVolume(:,1);
% 
% LoadbyLength    = (rho.*Volume_comp).*1e-3;%./Length_comp;
% 
% figure(13)
% barh(LoadbyLength)
% % barh(Length_comp)
% title(['Component weight (Large test-section) (tonnes)'])
% xlabel('Load per unit length [tonnes]')
% set(gca, 'YTick', 1:length(NameMain), 'YTickLabel', NameMain);
% grid on
% 











