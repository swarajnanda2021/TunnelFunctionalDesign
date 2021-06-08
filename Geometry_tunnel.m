clear all
close all
clc


TunnelParams = importtunnelparams('Tunnel_params_5.csv', 1, 14);
                



r = 1e-4;
resorberdiameter = [1.5];

diffuseroutlet = [0.9];
heightmaximum    = [10];
Dpump   = 0.8;

for aaa = 1:length(heightmaximum)  
for bbb = 1:length(diffuseroutlet)
for ccc = 1:length(resorberdiameter)
%% Evaluate expressions in the loaded data

Hmax = heightmaximum(aaa);
DiffMax = diffuseroutlet(bbb);
Lbottom = resorberdiameter(ccc);
Rdis = DiffMax*2 ;

%%
for i=1:length(TunnelParams(:,1))
    for j=1:length(TunnelParams(1,:))
        if (i>2) && (i<10)
            if isnan(str2double(TunnelParams(i,j))) == 1
                TunnelParams{i,j} = string(eval(char(TunnelParams{i,j})));  
            end
        end
    end
end














% Import main tunnel parameters
NN = 1000;
HTs = 0.3;          % Height of the tunnel test-section      
% TunnelParams = importtunnelparams('Tunnel_params.csv', 1, 14);
% Use the geom function to output a component matrix for both the main tunnel and the Pipe section.
[ComponentsMain,~,~,TunnelVolume, F0byFiMain, MainArInlet] = geom_tunnel(TunnelParams,NN);

% Import pipe section parameters
PipeParams = importpipeparams('Pipe_params.csv',1,14);
% Use the geom function to output a component matrix for both the main tunnel and the Pipe section.
[ComponentsPipe,intake,remerge,PipeVolume, F0byFiPipe, PipeArInlet] = geom_tunnel(PipeParams,NN);


% Plot geometry parameters
figure(1)
subplot(2,1,1)
title('Main tunnel')
xlabel('l [m]')
hold on

yyaxis left
ylabel('A_{cross} [m^2]')
plot(ComponentsMain(1,:),ComponentsMain(2,:))

yyaxis right
ylabel('H_S [m]')
plot(ComponentsMain(1,:),ComponentsMain(5,:))

hold off
grid on

subplot(2,1,2)
title('Pipe section')
xlabel('l [m]')
hold on

yyaxis left
ylabel('A_{cross} [m^2]')
plot(ComponentsPipe(1,:),ComponentsPipe(2,:))

yyaxis right
ylabel('H_S [m]')
plot(ComponentsPipe(1,:),ComponentsPipe(5,:))

hold off
grid on
% Plotting tunnel volumes
figure(2)
subplot(1,2,1)
barh(TunnelVolume(:,1))
title(['Main Tunnel (' num2str(sum(TunnelVolume(:,1))) ' m^3)'])
xlabel('Volume [m^3]')
% ylim([min(TunnelVolume) max(TunnelVolume)])
set(gca, 'YTick', 1:32, 'YTickLabel', TunnelParams(1,:));
grid on
set(gca,'xscale','log')
% xtickangle(90)
subplot(1,2,2)
barh(PipeVolume(:,1))
set(gca, 'YTick', 1:13, 'YTickLabel', PipeParams(1,:));
set(gca,'xscale','log')
%xtickangle(90)
title(['Pipe Section (' num2str(sum(PipeVolume(:,1))) ' m^3)'])
% ylim([min(PipeVolume) max(PipeVolume)])
xlabel('Volume [m^3]')
grid on
%% Applying Continuity to calculate flow-rate in the Main tunnel and the Pipe section
% Remember: Since we are designing a closed loop system, the pump only has to 
% do work to overcome frictional losses in the circuit. Therefore, for a particular 
% discharge, the operating head of the pump must correspond to the frictional 
% losses in that flow-rate. 
% 
% We want to control the flow within the pipe-section by simply introducing 
% a solidity in the main tunnel, thereby running both sections simultaneously. 
% Noting the point of intake of the pipe-section and the point where the pipe-section 
% rejoins the main tunnel, the flow-rate in the tunnel can be estimated for a 
% particular value of solidity in the main tunnel test-section. Therefore, the 
% tunnel is controlled by 2 factors, the solidity at the test-section and the 
% flow rate. When solidity is zero, we assume that the pipe-section is closed 
% off. When the solidity is non-zero, but less than 1, we introduce flow into 
% the pipe-section by continuity. 
% 
% More importantly, when the solidity is non-zero, a new friction factor 
% caused by the points of intake and re-combining of the pipe-section must be 
% accounted.
% 
% Also, when we use the pipe-section, since the pipe section is above the 
% main tunnel, a head difference between the main test-section centreline and 
% the pipe section centreline must be added.

% Calculating pressures from Bernoulli for given flow-rate and main test-section solidity

% Input flow parameters 
ValveDia     = 0.05;                                   % Valve height
theta        = deg2rad(0);                            % Valve opening angle
solidity     = sin(theta) * ValveDia / HTs;            % Test-section solidity
QPump        = 1.08;                                   % Pump output flowrate [m^2/hr]   
rho          = 1000;                                   % Fluid density [kg/m^3]
SetP         = 101325*0.05;                                   % Tunnel set pressure [Pascals] 
nu           = 1e-6;                                   % Fluid kinematic viscosity [m^2/sec]


% Dividing the flow-rates as per solidity
QTunnel      = zeros(size(ComponentsMain(2,:)));
QPipe        = ones(size(ComponentsPipe(2,:))) * (QPump * solidity); 
for i = 1:length(QTunnel)
   if ComponentsMain(1,i) < intake
       QTunnel(i) = QPump;
   elseif (ComponentsMain(1,i) >= intake && ComponentsMain(1,i) <= remerge)
       QTunnel(i) = QPump*(1-solidity);
   elseif ComponentsMain(1,i) > intake
       QTunnel(i) = QPump;
   end
end

%Plotting flow-rates
if solidity > 0
    figure(3)
    subplot(2,1,1)
    hold on
    plot(ComponentsMain(1,:),(QTunnel./ComponentsMain(2,:)))
    title('Fluid speeds and flowrates (Main tunnel)')
    xlabel('l [m]')
    ylabel('U [m/sec]')
    yyaxis right
    ylabel('Q [m^3/hr]')
    plot(ComponentsMain(1,:),QTunnel)
    grid on
    
    
    subplot(2,1,2)
    hold on
    plot(ComponentsPipe(1,:),(QPipe./ComponentsPipe(2,:)))
    title('Pipe section')
    xlabel('l [m]')
    ylabel('U [m/sec]')
    yyaxis right
    ylabel('Q [m^3/hr]')
    plot(ComponentsPipe(1,:),QPipe)
    grid on
    hold off
elseif solidity ==0
    figure(3)
    hold on
    plot(ComponentsMain(1,:),(QTunnel./ComponentsMain(2,:)))
    title('Fluid speeds and flowrates (Main tunnel)')
    xlabel('l [m]')
    ylabel('U [m/sec]')
    yyaxis right
    ylabel('Q [m^3/hr]')
    plot(ComponentsMain(1,:),QTunnel)
    grid on
end

% Calculating the Static pressure [Pa] for a particular operating fluid density
PStatMain = ComponentsMain(5,:) .* rho * 9.81;

% Calculating the Dynamic pressure [Pa] of the flow in both sections
PDynMain  = (0.5*rho).* (QTunnel./(ComponentsMain(2,:))).^2;

% Plotting the Static and Dynamic pressures
if solidity > 0
    PDynPipe  = (0.5*rho).* QPipe./(ComponentsPipe(2,:));
    PStatPipe = ComponentsPipe(5,:) .* rho * 9.81;

    figure(4)
    subplot(2,1,1)
    title('Main tunnel')
    xlabel('l [m]')
    hold on
    
    yyaxis left
    ylabel('P_{static} [Pa]')
    plot(ComponentsMain(1,:),PStatMain)
    
    yyaxis right
    ylabel('P_{Dynamic} [Pa]')
    plot(ComponentsMain(1,:),PDynMain)
    grid on
    hold off
    
    subplot(2,1,2)
    title('Pipe section')
    xlabel('l [m]')
    hold on
    
    yyaxis left
    ylabel('P_{static} [Pa]')
    plot(ComponentsPipe(1,:),PStatPipe)
    
    yyaxis right
    ylabel('P_{dynamic} [Pa]')
    plot(ComponentsPipe(1,:),PDynPipe)
    grid on
    hold off
elseif solidity ==0
    figure(4)
    title('Main tunnel')
    xlabel('l [m]')
    hold on
    
    yyaxis left
    ylabel('P_{static} [Pa]')
    plot(ComponentsMain(1,:),PStatMain)
    
    yyaxis right
    ylabel('P_{Dynamic} [Pa]')
    plot(ComponentsMain(1,:),PDynMain)
    grid on
    hold off
end

%% Calculating frictional losses for each component of the Main Tunnel as well as the Pipe section
% Theoretical frictional losses for standard components such as bends, expansions, 
% contractions and straight sections are calculated from the "Handbook of Hydraulic 
% Resistance" by I. E. Idelchik. Theoretical frictional losses for honeycombs 
% are obtained from the iso-reduction curves by John L. Lumley "Reducing water-tunnel 
% turbulence by means of a honeycomb (1967)." Theoretical frictional losses for 
% meshes are obtained from the studies of Kurian and Fransson in "Grid-generated 
% turbulence revisited."
% 
% These theoretical frictional losses provided an introductory basis for 
% choosing appropriate components to minimize losses and improve the energy-ratio 
% of the wind tunnel. Indeed, on further refinement of the design, the bends, 
% contractions, and diffusers are expected to have smaller friction factors. It 
% is important to note that the friction factors are calculated on the bases of 
% the flow-rate at that component, influenced by the pump discharge and the test-section 
% solidity.
% 
% Another noteworthy point is the flow-control device for the pipe-section. 
% Here, we assume a butterfly valve for throttling the flow in the main test-section 
% which is taken from Diagram 9-17 of Idelchik. It is important to note that the 
% throttling valve must be used in the range of 33% to 100% open. This allows 
% the use of Idelchik's expressions for a square throttling valve, since observed 
% friction factors appear to be independant of the valve shape and thickness.
% 
% 

% Applying the function friction_factor.m to estimate theoretical frictional losses

if solidity == 0     % Only the Main Tunnel is in use
    [ FricMain, QIMain ] = friction_factors( QTunnel, TunnelParams, ComponentsMain, NN, nu, QPump, max(QPipe), solidity,ValveDia);
elseif solidity > 0
    [ FricMain, QIMain ] = friction_factors( QTunnel, TunnelParams, ComponentsMain, NN, nu, QPump, max(QPipe), solidity,ValveDia);
    [ FricPipe, QIPipe ] = friction_factors( QPipe  , PipeParams  , ComponentsMain, NN, nu, QPump, max(QPipe), solidity,ValveDia);
    FricPipe(end) = []; %Remove last component of FricPipe, not needed
end

% Modify diffusers of the Main Tunnel to account for vanes (Reduction in \zeta to 65%). Comment out if not needed
% VaneCompListMain = [1,5,14,21];    
% for i=1:length(VaneCompListMain)
%     FricMain(VaneCompListMain(i)) = 0.65*FricMain(VaneCompListMain(i));    
% end

% Calculating individual loss components
if solidity > 0
    for i=1:length(FricMain)-1
        LossMain(i) = FricMain(i)*(0.5 * rho * (QIMain(i)/MainArInlet(i))^2);
    end
    LossMain = [LossMain (FricMain(end) * (0.5*rho*((QPump-QPipe(1))/(HTs^2)))) ];
    for i=1:length(FricPipe)
        LossPipe(i) = FricPipe(i)*(0.5 * rho * (QIPipe(i)/PipeArInlet(i))^2);
    end
    
elseif solidity == 0
    for i=1:length(FricMain)
        LossMain(i) = FricMain(i)*(0.5 * rho * (QIMain(i)/MainArInlet(i))^2) ;
    end
end

% % % Calculating system loss coefficients
% if solidity == 0 % Normalize w.r.t. main test-section jet KE
%     TestSectionJetKE = 0.5 * rho * (min(QTunnel)/min(ComponentsMain(2,:)))^2;
%     for i=1:length(FricMain)
%         
%         zeta0Main(i) = FricMain(i) * (HTs^4 /    MainArInlet(i)^2) * 1.5; % Correlction factor of 1.5
%     end
%     LossMain = sum(zeta0Main) * TestSectionJetKE;
% elseif solidity > 0 % Normalize w.r.t. pipe section jet KE
%     PipeSectionJetKE = 0.5 * rho * (max(QPipe)/min(ComponentsPipe(2,:)))^2;
%     for i=1:length(FricMain)
%         zeta0Main(i) = FricMain(i) * (min(ComponentsPipe(2,:))^2 / MainArInlet(i)^2);
%     end
%     for i=1:length(FricPipe)
%         zeta0Pipe(i) = FricPipe(i) * (min(ComponentsPipe(2,:))^2 / MainArInlet(i)^2);
%     end
%     LossMain = sum(zeta0Main) * PipeSectionJetKE;
%     LossPipe = sum(zeta0Pipe) * PipeSectionJetKE;
% end

% Plotting friction factors
if solidity > 0
    figure(6)
    grid on
    subplot(1,2,1)
    barh(FricMain)
    title(['Main Tunnel Head loss: ' num2str(LossMain/(9.81*rho)) ' m.'])
    xlabel('Global loss coeff')
    set(gca, 'YTick', 1:21, 'YTickLabel', [TunnelParams(1,:) 'Throttle']);
    grid on
    set(gca,'xscale','log')
    
    subplot(1,2,2)  
    barh(FricPipe)
    set(gca, 'YTick', 1:13, 'YTickLabel', PipeParams(1,:));
    grid on
    set(gca,'xscale','log')
    title(['Pipe Tunnel Head Loss: ' num2str(LossPipe) ' m.'])
    xlabel('Global loss coeff')
elseif solidity == 0
    figure(6)
    grid on
    barh(LossMain)
    title(['Main Tunnel Head loss: ' num2str(sum(LossMain)/(9.81*1000)) ' m.'])
    xlabel('HYKAT: Total pressure loss (Pascals)')
    set(gca, 'YTick', 1:24, 'YTickLabel', TunnelParams(1,:));
    grid on
    set(gca,'xscale','log')
end
%% Susceptibility to cavitation of the tunnel: Using bernoulli's theorem to calculate pressure along a streamline
% We do not know what the fluid pressure at the outlet to the pump is. Given 
% this, and taking this pressure to be zero, we can calculate $C_p$ across the 
% tunnel length to help design a tunnel without dangerously low $C_p$ values. 
% Once $C_p$ goes below zero, we expect to be concerned since the pressure in 
% the fluid will be lower than that at the pump outlet. Similarly for the pipe 
% section. In both cases, we non-dimensionalize with respect to the main tunnel 
% test-section jet kinetic energy.


% Using Bernoulli to find the constant energy at the test-section
BernoulliConstant = max(PDynMain) + SetP ; % Bernoulli constant is determined with respect to set pressure at the test-section, elevation is not considered



TestSectionJetKE  = (0.5 * rho * ((QPump - max(QPipe))/(HTs^2))^2);
% Calculating $C_p$

%% Estimating Required Net Positive Suction Head, Pump Power and investigating the need of Heat Exchangers
% For the net positive suction head, we require the calculation of loss head, 
% absolute head by tunnel atmosphere, static head and the velocity head. When 
% solidity is non-zero, the local friction factors are made global by using the 
% jet kinetic energy of the pipe-section, while when the solidity is 0 (the throttle 
% valve is fully open), the non-dimensionalization is with respect to the main 
% test-section jet kinetic energy.

if solidity == 0
    PLiquidMain = BernoulliConstant - PDynMain;
    
    UTs     = QPump / (HTs^2);
    H_abs   = SetP/(rho * 9.81);                                % Absolute head above the free-surface
    H_z     = (PStatMain(end))/(rho * 9.81);                      % Static pressure above the pump, the pump is working on positive static head conditions
    H_vp    = 3173.0724/(rho * 9.81);                           % Water Vapor Pressure at 298K.
    H_v     = ((QTunnel(end)/ComponentsMain(2,end))^2)/(2*9.81);% Velocity head at pump section
    H_f     = (sum(LossMain))/(rho * 9.81);     % Total frictional head
    NPSH_a = H_abs + H_z - (H_f + H_v) - H_vp;                  % Available Net Positive suction head
    % Display warning for low NPSH at the pump
    if NPSH_a < 0
        disp(['NPSH available at the pump is below zero. Provide atleast ' num2str(abs(NPSH_a)) ' m of static head more, upstream of the pump.'])
    end
elseif solidity > 0
    DiaPipe = min(ComponentsPipe(3,:));
    UTs     = QPipe(1) / (DiaPipe^2);
    H_abs   = (SetP - min(PStatPipe))/(rho * 9.81);                                % Absolute head above the free-surface
    H_z     = PStatPipe(end)/(rho * 9.81);                      % Static pressure above the pump, the pump is working on positive static head conditions
    H_vp    = 3173.0724/(rho * 9.81);                           % Water Vapor Pressure at 298K.
    H_v     = ((QTunnel(end)/ComponentsMain(2,end))^2)/(2*9.81);% Velocity head at pump section
    H_f     = (sum(LossMain) + sum(LossPipe))*(UTs^2)/(rho * 9.81);     % Total frictional head
    NPSH_a  = H_abs + H_z - (H_f + H_v) - H_vp;                  % Available Net Positive suction head
    % Display warning for low NPSH at the pump
    if NPSH_a < 0
        disp(['NPSH available at the pump is below zero. Provide atleast ' num2str(abs(NPSH_a)) ' m of static head more, upstream of the pump.'])
    end 
end
%% 
% With the required NPSH estimated, and with a knowledge of the theoretical 
% frictional losses, it is possible to calculate the required pump power for a 
% given pump efficiency. 

if solidity == 0
    EtaPump             = 0.7;
    RequiredPumpPower   = (sum(FricMain)*(rho/2) * UTs^2) * QPump / EtaPump;                 % Power in Watts   
elseif solidity >0
    EtaPump             = 0.7;
    RequiredPumpPower   = (sum(FricMain) + sum(FricPipe))*(rho/2) * UTs^2 * QPump / EtaPump; % Power in Watts
end
%% 
% Now, since we have assumed a pump efficiency, the lost power is expectedly 
% converted to dissipative heat. By using the specific heat of the operating fluid 
% and the tunnel capacity, one can estimate the heat capacity of the tunnel and 
% from the dissipated heat in KW, one can find the rise in temperature of the 
% fluid caused by the pumps.

CpFluid = 4.186; % Specific head of operating fluid in kW/kg-sec oC
% Estimating rise in temperature every second
if solidity == 0
    TempRise   =  (sum(TunnelVolume(:,1))*rho) * CpFluid/RequiredPumpPower;                 % oC per second
elseif solidity >0
    TempRise   =  ((sum(TunnelVolume(:,1)) + sum(PipeVolume(:,1)))*rho) * CpFluid/RequiredPumpPower;                 % oC per second
end

%% Estimating turbulence reduction factors for turbulence management devices in both sections
% Turbulence reduction factors are empirically calculated for each turbulence 
% management device. Due to the high speed in the tunnel and the susceptibility 
% of noise caused by screens, turbulence management devices are restricted to 
% the mesh, the honeycomb and the contraction (otherwise the hydroacoustic integrity 
% of the tunnel is severely compromised). The expressions for turbulence intensity 
% reduction factors have been derived in theory and validated through experiments 
% for each of the selected components. The theory is predominantly from Batchelor's 
% rapid distortion analysis for the Contraction, Lumley's approach to estimating 
% the honeycomb turbulence reduction factor and recent works on grid generated 
% turbulence by Kurian and Fransson.
% 
% 

% Using the function TurbRed to estimate turbulence reduction factor for main and pipe section. Does not require flow-rate values
TI_Inlet = linspace(0.01,0.30,5);    % Inlet turbulence intensity array
if solidity == 0
    TIEvolMain = [];
    for i=1:length(TI_Inlet)
        [ TIEvolMain(:,i) ] = TurbRed( QIMain, TunnelParams, TI_Inlet(i), NN, nu );
    end
elseif solidity >0
    TIEvolMain = [];
    TIEvolPipe = [];
    for i=1:length(TI_Inlet)
        [ TIEvolMain(:,i) ] = TurbRed( QIMain, TunnelParams, TI_Inlet(i), NN, nu );
        [ TIEvolPipe(:,i) ] = TurbRed( QIPipe, PipeParams  , TI_Inlet(i), NN, nu );
    end
    
end
    
% Plotting turbulence reduction performance

if solidity > 0 
    figure(8)
    subplot(2,1,1)
    for j = 1:length(TI_Inlet)
            plot(TIEvolMain(:,j))
            hold on
    end
    set(gca,'yscale','log')
    grid on
    title('Turbulence intensity (Main Tunnel)')
    xlabel('Component no.')
    ylabel('TI evolution')
    subplot(2,1,2)
    for j = 1:length(TI_Inlet)
            plot(TIEvolPipe(:,j))
            hold on
    end
    set(gca,'yscale','log')
    grid on
    title('(Pipe section)')
    xlabel('Component no.')
    ylabel('TI evolution')
    
elseif solidity == 0
    figure(8)
    for j = 1:length(TI_Inlet)
            hold all
            plot(TIEvolMain(:,j))
    end
    set(gca,'yscale','log')
    grid on
    title('Turbulence intensity (Main Tunnel)')
    xlabel('Component no.')
    ylabel('TI evolution')
    
end



% eta(aaa,bbb,ddd,eee) = TIEvolMain(3,3)/TIEvolMain(1,3); 
% LCell(aaa,bbb,ddd,eee) = l_cell(aaa);
% MH(aaa,bbb,ddd,eee) = M_h(bbb);
% LS(aaa,bbb,ddd,eee) = l_s(ddd);
% LIn(aaa,bbb,ddd,eee) = L_in(eee);



%% Evaluating the dissolution and separation quality in the Main Tunnel
% Any cavitation tunnel must be designed to ensure that any free-gas introduced 
% in the test-section vanishes along the tunnel length. Traditionally, this process 
% is achieved through gravity separation of large bubbles (> 100 microns) and 
% the dissolution of small bubbles (< 100 microns) principally because switching 
% these processes between the 2 ranges of bubbles described increases the separation 
% and dissolution times by orders of magnitudes. 
% 
% For the dissolution of microbubbles, it shall be assumed that there is 
% no slip between the bulk fluid and a bubble of certain size (which holds true 
% for bubbles smaller than a millimeter), otherwise convective mass transfer would 
% drastically enhance the dissolution processes. In such a case, the dissolution 
% model proposed by Epstein and Plesset shall be utilized. While the model is 
% applicable to pure-systems (pure water), workers have shown that this model 
% could be applicable for contaminated water as well, althought this would be 
% highly uncertain in reality. The issue of contaminated water has to do with 
% the choice of seawater as the operating fluid. Therefore, 2 models, i.e. the 
% Epstein Plesset model and the model by Levich shall be used here. The objective 
% is to ensure that the tunnel is long enough to dissolve 100 micron air-bubbles.
% 
% For the separation of macrobubbles, it suffices to indicate the maximum 
% vertical distance travelled by the bubble over a specified length. We assume 
% the bubble starts at the bottom of the test-section and rises by gravity as 
% the flow continues along the upper limb of the tunnel. We simply calculate whether 
% the bubble has crossed the vertical dimension of the test-section along the 
% upper limb. It muse be noted, however, that when there are horizontal space 
% constraints, the bubble may not travell sufficiently high before encountering 
% a bend. In such a case, it may be natural to then provide a separation length 
% along an intermediary limb, just like the EPFL tunnel. In such a case it should 
% be noted that the bend will negate the progress in bubble-rise made in the upper 
% limb as the bend would turn by 180 degrees. The separation length in the intermediary 
% limb, then, must be designed for the entire rise range. For seawater, the terminal 
% velocity curve-fit of Detsch will be used, while that of Karamanev will be used 
% for relatively pure water.



% Modify geometry definition to start from separator
Separator    = 21; % Separator component number
for i=1:length(ComponentsMain(:,1))
    ComponentsMain(i,:) = circshift(ComponentsMain(i,:), [0 -(NN*(Separator-1))]);
end
changing=size(TunnelParams(1,:),2)-Separator + 1;
ComponentsMain(1,changing*NN+1:end) = ComponentsMain(1,changing*NN+1:end) + (ComponentsMain(1,changing*NN));


PStatMain      = circshift(PStatMain, [0 -(NN*(Separator-1))]);
PLiquidMain    = circshift(PLiquidMain, [0 -(NN*(Separator-1))]);
QTunnel        = circshift(QTunnel, [0 -(NN*(Separator-1))]);
%% Calculating time of flow along tunnel length
TimeTunnel = zeros(1,length(ComponentsMain(2,:)));
for ii=1:(length(ComponentsMain(2,:))-1)
    
    dx = abs(ComponentsMain(1,ii) - ComponentsMain(1,ii+1));
    U(ii)  = QTunnel(ii)/ComponentsMain(2,ii);
    TimeTunnel(ii+1) = TimeTunnel(ii)+ (dx/U(ii));
    
end
%% Dissolution
% delete(findall(0,'Type','figure'))

%% Dissolution

% Tunnel Parameters
% r       = 1*1e-4;               % Initial bubble radius [m]
f       = 0.2;                  % Degassing factor [C_i/C_S]
PInfty  = PStatMain +  PLiquidMain; 

% Test-section component number
TestSection = 21;

% Fluid and Model Parameters
D_AB    = 2.82e-9;              % Diffusion constant (Air-Water)[m^2/sec]
rho_b   = 1.225;                 % Density of bubble (Air) [kg/m^3]    
mu_b    = 1.81e-5;              % Dynamic viscosity of air
fluid   = {'Pure'; 'SeaWater'};               % Working fluid (Pure/SeaWater)             
for ii=1:length(fluid)
        
    if strcmp(fluid(ii),'SeaWater') == 1
        
        % Calculating solubility of air in sea water at P_infty
        K_H_O2      = 1.2e-5;    % Henry's constant for oxygen in air [mol/m^3Pa] % Check Zotero for link to paper
        K_H_N2      = 6.4e-6;    % Henry's constant for nitrogen in air [mol/m^3Pa]
        P_O2_air    = 0.21;      % Partial fraction for O2 and N2 in air
        P_N2_air    = 0.79;
        M_O2        = 31.99e-3;  % Molar mass of Oxygen and Nitrogen [kg/mol]
        M_N2        = 28.0134e-3;

        C_S = ((M_O2 * P_O2_air * K_H_O2) .* PInfty) + (M_N2 * P_N2_air * K_H_N2) .* PInfty;     
        Salinity    = 35 ;                                                                      % Salinity in grams/kilogram of water
        temp        = 25;                                                                       % Temperature in degrees celcius
    
%       sigma       = SW_SurfaceTension(temp,'C',Salinity,'ppt') * 1e-3 * 0.3;                        % Correlation taken from the paper 'surface tension of seawater by L.G. Nayar (MIT)
        sigma = 0.022; %Lower limit of Lozano's measurements
        rho_f       = SW_Density(temp,'C',Salinity,'ppt',1,'bar');   
        mu_f        = SW_Viscosity(temp,'C',Salinity,'ppt');
        M_a         = 0.29e-3;                                              % Molar mass of air [kg/mol]

        % Terms of the Epstein Plesset model
        d           = C_S*(rho_b^(-1));
        b           = (d./(2*pi)).^0.5;
    
        % Non-dimensionalising time
    
        for j = 1:length(d)
            XprimeBernoulli(j) = sqrt((2 * D_AB * d(j)/r^2) * TimeTunnel(j));
        end
    
    
        
        % Solving Epstein-Plesset's model
        epsilonini = 1;
        RSaltyEvolution = [];
        RSaltyEvolution = [RSaltyEvolution r];
        D_calc = [];
        RNonDimNext = 1;
        TNonDimEvolution = [];
        RNonDimEvolution = [];
        for i=1:length(TunnelParams(1,:))
            if i < ((TestSection-1)) % Terminate before test-section
            RLoop = RSaltyEvolution(i);
 %             delta       = (2 * M_a * sigma/(8.314*298))/(RLoop * rho_f);
            rho_infty = (M_a/(8.314*298)) .* PInfty(((i-1)*NN + 1):(i*NN));
            delta       = ((2 * M_a * sigma/(8.314*298))/(RLoop)) ./ rho_infty;
            opts                = odeset('NonNegative',1);     % No negative terms
            [TNonDim,RNonDim]   = ode45(@(t_nd,r_nd) epstein_plesset_nondim_2(t_nd,r_nd, XprimeBernoulli(((i-1)*NN + 1):(i*NN)), f , b(((i-1)*NN + 1):(i*NN)), delta), [XprimeBernoulli(((i-1)*NN + 1)) XprimeBernoulli((i*NN))], RNonDimNext,opts);
            RNonDimNext = RNonDim(end);
            Rnew = RNonDim(end)*RLoop;
            RSaltyEvolution = [RSaltyEvolution Rnew];
            
            TNonDimEvolution = [TNonDimEvolution TNonDim'];
            RNonDimEvolution = [RNonDimEvolution RNonDim'];
            % Redimensionalizing Time
            D_calc = [D_calc (interp1(XprimeBernoulli(((i-1)*NN + 1):(i*NN)),d(((i-1)*NN + 1):(i*NN)),TNonDim))'];
            
            end
        end
    
        RSaltyDimEvol = RNonDimEvolution.*r;
        TSaltyDimEvol = (TNonDimEvolution).^2 ./ ((2 .* D_AB .*D_calc ./r.^2));

    end
    
    if strcmp(fluid(ii),'Pure') == 1
        % Calculating solubility of air in water at P_infty
        K_H_O2      = 5.51e-4 / 101.630;    % Henry's constant for oxygen in air [mol/m^3Pa]
        K_H_N2      = 1.1e-3  / 101.630;    % Henry's constant for nitrogen in air [mol/m^3Pa]
        P_O2_air    = 0.21;                 % Partial fraction for O2 and N2 in air
        P_N2_air    = 0.79;
        M_O2        = 31.99e-3;             % Molar mass of Oxygen and Nitrogen [kg/mol]
        M_N2        = 28.0134e-3;
    
        C_S = ((M_O2 * P_O2_air * K_H_O2) .* PInfty) + (M_N2 * P_N2_air * K_H_N2) .* PInfty;     
        sigma       = 71.99e-3;          % Surface tension for air in water at 298 K [N/m]
        rho_f       = 1000;              % Density of water [kg/m^3]
        mu_f        = 1e-3;              % Dynamic viscosity of pure water
    
        M_a         = 0.29e-3;                                              % Molar mass of air [kg/mol]
    
        % Terms of the Epstein Plesset model
        d           = C_S*(rho_b^(-1));
        b           = (d./(2*pi)).^0.5;
    
            % Non-dimensionalising time
    
        for j = 1:length(d)
            XprimeBernoulli(j) = sqrt((2 * D_AB * d(j)/r^2) * TimeTunnel(j));
        end
    
    
    
        % Solving Epstein-Plesset's model
        epsilonini = 1;
        REvolution = [];
        REvolution = [REvolution r];
        D_calc = [];
        RNonDimNext = 1;
        TNonDimEvolution = [];
        RNonDimEvolution = [];
        for i=1:length(TunnelParams(1,:))
            
            if i < ((TestSection-1)) % Terminate before test-section
            RLoop = REvolution(i);
%             delta       = (2 * M_a * sigma/(8.314*298))/(RLoop * rho_f);
            rho_infty = (M_a/(8.314*298)) .* PInfty(((i-1)*NN + 1):(i*NN));
            delta       = ((2 * M_a * sigma/(8.314*298))/(RLoop)) ./ rho_infty;
            opts                = odeset('NonNegative',1);     % No negative terms
            [TNonDim,RNonDim]   = ode45(@(t_nd,r_nd) epstein_plesset_nondim_2(t_nd,r_nd, XprimeBernoulli(((i-1)*NN + 1):(i*NN)), f , b(((i-1)*NN + 1):(i*NN)), delta), [XprimeBernoulli(((i-1)*NN + 1)) XprimeBernoulli((i*NN))], RNonDimNext,opts);
%             figure(i)
%             plot(TNonDim,RNonDim)
%             drawnow
            
            
            RNonDimNext = RNonDim(end);
            Rnew = RNonDim(end)*RLoop;
            REvolution = [REvolution Rnew];
    
    
    
            TNonDimEvolution = [TNonDimEvolution TNonDim'];
            RNonDimEvolution = [RNonDimEvolution RNonDim'];
    
            % Redimensionalizing Time
            D_calc = [D_calc (interp1(XprimeBernoulli(((i-1)*NN + 1):(i*NN)),d(((i-1)*NN + 1):(i*NN)),TNonDim))'];
            end
        end
        RDimEvol = RNonDimEvolution.*r;
        TDimEvol = (TNonDimEvolution).^2 ./ ((2 .* D_AB .*D_calc ./r.^2));
    
    end
       
end
    
            
% Plotting
figure(9)
plot(TSaltyDimEvol,RSaltyDimEvol)
hold on
plot(TDimEvol,RDimEvol)
hold off
% plot(1:length(TunnelParams(1,:))+1,REvolution')    
% hold on
% plot(1:length(TunnelParams(1,:))+1,RNewLevich')
% hold off
title('Behaviour of bubble along the tunnel length')
xlabel('Time (seconds)')
ylabel('Bubble radius (m)')
legend('Natural sea-water','Pure')
% Rise
% We avoid the notion of bubble dissolution in the calculations. We want to
% design the loop for a stationary bubble rise only. In this case, the
% simulation is conducted for finding if the bubble rises above the
% cross-section along the length.
RiseStart = 1;
RiseEnd   = 2;

Ut = [];
RiseParamsStokes    = {'Stokes' num2str(rho_b) num2str(rho_f) num2str(mu_b) num2str(sigma)}';
RiseParamsKaramanev = {'Karamanev' num2str(rho_b) num2str(rho_f) num2str(mu_f) num2str(sigma)}';
RiseParamsDetschPure= {'Detsch' 'PureWater' num2str(rho_f) num2str(mu_b) num2str(sigma)}';
RiseParamsDetschSW  = {'Detsch' 'SeaWater' num2str(rho_f) num2str(mu_b) num2str(sigma)}';

Ut = [Ut  rise_2(r,RiseParamsStokes) rise_2(r,RiseParamsKaramanev) rise_2(r,RiseParamsDetschSW) rise_2(r,RiseParamsDetschPure)];
BubTime ={};
BubRise ={};
for i=1:length(Ut)
    opts = odeset('NonNegative',1);     % No negative terms
    [t_r,Y]=ode45(@(t,y) dydt_terminal(t,y,Ut(i)), [TimeTunnel((RiseStart-1)*NN + 1) TimeTunnel(RiseEnd*NN)], 0,opts);
    
    FinalHeight(i) = Y(end);
end
figure(10)
bar(FinalHeight)
ylabel('Maximum vertical displacement [m]')
set(gca, 'XTick', 1:4, 'XTickLabel', {'Stokes' 'Karamanev' 'Detsch (Sea)' 'Detsch (Pure)'});










% %%
% % 
% delete(findall(0,'Type','figure'))


% save(['R'   num2str(r*1e6) 'Lv'  num2str(Hmax) 'Lb' num2str(Lbottom) 'Ds' num2str(Rdis) 'Dr' num2str(Lbottom) 'Dp' num2str(Dpump) 'Dm' num2str(DiffMax) '.mat'])




% save(['TLc' num2str(l_cell(aaa)*100) 'Mh' num2str(M_h(bbb)*100) 'Ls' num2str(l_s(ddd)*10) 'LhdbyM' num2str(40) 'Lin' num2str(L_in(eee)*100) '.mat'])

end
end
end


% save('TurbParametric_xbyD=40')