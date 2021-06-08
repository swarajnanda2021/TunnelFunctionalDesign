% function [TunnelParams,LossMain,TunnelVolume,EnergyRatio,NPSH_a,TIEvolMain,SepVel,Finaldiabubble,ContractionWidth_stage] = TunnelMainCode( Stage, HTs,r,Hmax,Hmin,SepIn,SepOut,HeightSec,WidthSec,HeightSep,WidthSep,Dpump,Lpump,Lhoneycomb,Lsettling,Lcontraction,Lts,Lmaindiffuser,MainDiffOut,Lpumpleft,RT1,RT2,RT3,RT4,RT5,MeshWidthHoneycomb,CellLengthHoneycomb,LinletHoneycomb,QPUMP,Dpipe,Buffer,lunchboxwidth,CRpipe,theta ,RoughMain,thetaangle,betaangle,CR1x,CR1y,CR2,CStage1Length,CStage2Length,SetP,Lseparator,RiseStart,RiseEnd,Separator,TestSection,VaneCompListMain,Insertlength,NameMain,NamePipe)
% function [TunnelParams,LossMain,TunnelVolume,EnergyRatio,NPSH_a,TIEvolMain] = TunnelMainCode( Stage, HTs,r,Hmax,Hmin,SepIn,SepOut,HeightSec,WidthSec,HeightSep,WidthSep,Dpump,Lpump,Lhoneycomb,Lsettling,Lcontraction,Lts,Lmaindiffuser,MainDiffOut,Lpumpleft,RT1,RT2,RT3,RT4,RT5,MeshWidthHoneycomb,CellLengthHoneycomb,LinletHoneycomb,QPUMP,Dpipe,Buffer,lunchboxwidth,CRpipe,theta ,RoughMain,thetaangle,betaangle,CR1x,CR1y,CR2,CStage1Length,CStage2Length,SetP,Lseparator,RiseStart,RiseEnd,Separator,TestSection,VaneCompListMain,Insertlength,NameMain,NamePipe,LinletCollector,LhoseInlet,LinletExpander,LoutletCollector,LhoseOutlet,LoutletExpander,DiaHose,Rtpipe,Ldynamometer,Lmaindiffuser2,PipeIntercept,PipeRejoin)
function [TunnelParams,LossMain,TunnelVolume,NPSH_a,TIEvolMain,FinalHeight] = ...
    TunnelMainC7(Config,Hmin,Dpump,Lhoneycomb,Lsettling,Lcontraction,lmd,Wexpansion...
        ,r,MeshWidthHoneycomb,CellLengthHoneycomb,LinletHoneycomb,QPUMP,RoughMain,SetP,...
        Lseparator,RiseStart,RiseEnd,NameMain,VaneCompListMain,ldr,lsa,l2,ladj,nBanks,nPlates,r_separation,...
        Lcontraction2,Ldynamometer,Lmaindiffuser2,lgap,Test_Section,nPlatesFunctional,L_KnifeCasing,alpha_opening,l_vent)
            


if Config == 1
    TunnelParams = importtunnelparams('Tunnel_params_large.csv', 1, 14,Config);
    HTs          = 0.5;
elseif Config == 2
    TunnelParams = importtunnelparams('Tunnel_params_cav.csv', 1, 14,Config);
    HTs          = 0.3;
end

%% Configurations 1,2 and 3


    
%% Evaluate expressions in the loaded data 

% Evaluate expressions in the parametrized format from the Excel files with
% exceptional components treated specially
for i=1:length(TunnelParams(:,1))
    for j=1:length(TunnelParams(1,:))
        % For all general parts
        if (i>2) && (i<10)
            TunnelParams{i,j} = string(eval(char(TunnelParams{i,j})));  
        end
        
        % For honeycomb
        if (strcmp(TunnelParams(1,j),'Honeycomb') == 1) && (i>2)
            if isnan(str2double(TunnelParams(i,j))) == 1
                TunnelParams{i,j} = string(eval(char(TunnelParams{i,j})));  
            end
        end
        
        % For honeycomb
        if (strcmp(TunnelParams(1,j),'Separator') == 1) && (i>2)
            if isnan(str2double(TunnelParams(i,j))) == 1
                TunnelParams{i,j} = string(eval(char(TunnelParams{i,j})));  
            end
        end
        
        % For exhaust knife
        if (strcmp(TunnelParams(2,j),'ductedProtuberance') == 1) && (i>2)
            if isnan(str2double(TunnelParams(i,j))) == 1
                TunnelParams{i,j} = string(eval(char(TunnelParams{i,j})));  
            end
        end
    end
end


% Import main tunnel parameters
NN = 1000;
% TunnelParams = importtunnelparams('Tunnel_params.csv', 1, 14);
% Use the geom function to output a component matrix for both the main tunnel and the Pipe section.
[ComponentsMain,~,~,TunnelVolume, F0byFiMain, MainArInlet,ContractionWidth_stage,lengths] = geom_tunnel(TunnelParams,NN,HTs);


% Plot geometry parameters
figure(1)
% subplot(2,1,1)
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

% Plotting tunnel volumes
figure(2)
% subplot(1,2,1)
barh(TunnelVolume(:,1))
title(['Main Tunnel (' num2str(sum(TunnelVolume(:,1))) ' m^3)'])
xlabel('Volume [m^3]')
% ylim([min(TunnelVolume) max(TunnelVolume)])
set(gca, 'YTick', 1:length(NameMain), 'YTickLabel', NameMain);
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

    % ValveDia     = 0.05;                                   % Valve height
    % solidity     = sin(theta) * ValveDia / HTs;            % Test-section solidity

    QPump        = QPUMP;                                   % Pump output flowrate [m^2/hr]   
    rho          = 1000;                                   % Fluid density [kg/m^3]
    nu           = 1e-6;                                   % Fluid kinematic viscosity [m^2/sec]

    QTunnel = ones(size(ComponentsMain(2,:))) * QPump;
   



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



    [ FricMain, QIMain ] = friction_factors( QTunnel, TunnelParams, ComponentsMain, NN, nu, QPump,0,RoughMain  );    



    for i=1:length(VaneCompListMain)
        FricMain(VaneCompListMain(i)) = 0.65*FricMain(VaneCompListMain(i));    
    end

    % end

    % Calculating individual loss components

   
    for i=1:length(FricMain)
        LossMain(i) = FricMain(i)*(0.5 * rho * (QIMain(i)/MainArInlet(i))^2) ;
    end
    
    LossMain(Test_Section) = 2*LossMain(Test_Section);

    % Plotting friction factors

    % elseif solidity == 0
    figure(6)
    grid on
    barh(LossMain)
    title(['Main Tunnel Head loss: ' num2str(sum(LossMain)/(9.81*rho)) ' m.'])
    xlabel('Pressure loss (Pascals)')
    set(gca, 'YTick', 1:length(NameMain), 'YTickLabel', NameMain);
    grid on
    %     set(gca,'xscale','log')
    % end


    %% Plotting of Static and Dynamic Pressures, Velocities, Areas and Static Heads
    %Plotting flow-rates

    % elseif solidity ==0
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
    % end

    %Calculating the Static pressure [Pa] for a particular operating fluid density
    PStatMain = ComponentsMain(5,:) .* rho * 9.81;

    %Calculating the Dynamic pressure [Pa] of the flow in both sections
    PDynMain  = (0.5*rho).* (QTunnel./(ComponentsMain(2,:))).^2;

    % Plotting the Static and Dynamic pressures

    % elseif solidity ==0
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
    % end



    %% Using bernoulli's theorem to calculate pressure along a streamline


    % Using Bernoulli to find the constant energy at the test-section
    BernoulliConstant = max(PDynMain) + SetP ; % Bernoulli constant is determined with respect to set pressure at the test-section, elevation is not considered

%     TestSectionJetKE  = (0.5 * rho * ((QPump - max(QPipe))/(HTs^2))^2);
    % Calculating $C_p$

    %% Estimating Required Net Positive Suction Head, Pump Power and investigating the need of Heat Exchangers
    % For the net positive suction head, we require the calculation of loss head, 
    % absolute head by tunnel atmosphere, static head and the velocity head. When 
    % solidity is non-zero, the local friction factors are made global by using the 
    % jet kinetic energy of the pipe-section, while when the solidity is 0 (the throttle 
    % valve is fully open), the non-dimensionalization is with respect to the main 
    % test-section jet kinetic energy.

    % if solidity == 0
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

    %% 
    % With the required NPSH estimated, and with a knowledge of the theoretical 
    % frictional losses, it is possible to calculate the required pump power for a 
    % given pump efficiency. 

    % if solidity == 0
        EnergyRatio = sum(LossMain)/((rho/2) * UTs^2);
        EtaPump             = linspace(0.4,0.8,100);
        for pumpeff=1:length(EtaPump)
            RequiredPumpPower(pumpeff)   = (sum(FricMain)*(rho/2) * UTs^2) * QPump / EtaPump(pumpeff);                 % Power in Watts   

        end

    %% 
    % Now, since we have assumed a pump efficiency, the lost power is expectedly 
    % converted to dissipative heat. By using the specific heat of the operating fluid 
    % and the tunnel capacity, one can estimate the heat capacity of the tunnel and 
    % from the dissipated heat in KW, one can find the rise in temperature of the 
    % fluid caused by the pumps.

    CpFluid = 4.186; % Specific head of operating fluid in kW/kg-sec oC
    % Estimating rise in temperature every second
    % if solidity == 0
        for temprise=1:length(RequiredPumpPower)
            TempRise(temprise)   =  (sum(TunnelVolume(:,1))*rho) * CpFluid*1000/RequiredPumpPower(temprise);                 % oC per second
        end
    % elseif solidity >0
    %     for temprise=1:length(RequiredPumpPower)
    %         TempRise(temprise)   =  ((sum(TunnelVolume(:,1)) + sum(PipeVolume(:,1)))*rho) * CpFluid*1000/RequiredPumpPower(temprise);                 % oC per second
    %     end
    % end



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
    for i=1:length(TI_Inlet)
        [ TIEvolMain(:,i) ] = TurbRed( QIMain, TunnelParams, TI_Inlet(i), NN, nu );
    end

    % Plotting turbulence reduction performance

    % elseif solidity == 0
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

    
    %% Calculating time of flow along tunnel length
    TimeTunnel = zeros(1,length(ComponentsMain(2,:)));
    for ii=1:(length(ComponentsMain(2,:))-1)

        dx(ii) = abs(ComponentsMain(1,ii) - ComponentsMain(1,ii+1));
        U(ii)  = QTunnel(ii)/ComponentsMain(2,ii);
        TimeTunnel(ii+1) = TimeTunnel(ii)+ (dx(ii)/U(ii));

    end


    %% Rise
    % We avoid the notion of bubble dissolution in the calculations. We want to
    % design the loop for a stationary bubble rise only. In this case, the
    % simulation is conducted for finding if the bubble rises above the
    % cross-section along the length.
    % Change r to r_separation
    r = r_separation;
    rho_b = 1.225;
    sigma       = 71.99e-3;          % Surface tension for air in water at 298 K [N/m]
    rho_f       = 1000;              % Density of water [kg/m^3]
    mu_f        = 1e-3;
    mu_b        = 1.81e-5;
    Ut = [];
    RiseParamsStokes    = {'Stokes' num2str(rho_b) num2str(rho_f) num2str(mu_f) num2str(sigma)}';
    RiseParamsKaramanev = {'Karamanev' num2str(rho_b) num2str(rho_f) num2str(mu_f) num2str(sigma)}';
    RiseParamsDetschPure= {'Detsch' 'PureWater' num2str(rho_f) num2str(mu_b) num2str(sigma)}';
%     RiseParamsDetschSW  = {'Detsch' 'SeaWater' num2str(rho_f) num2str(mu_b) num2str(sigma)}';

%     Ut = [Ut  rise_2(r,RiseParamsStokes) rise_2(r,RiseParamsKaramanev) rise_2(r,RiseParamsDetschSW) rise_2(r,RiseParamsDetschPure)];
    Ut = [Ut  rise_2(r,RiseParamsStokes) rise_2(r,RiseParamsKaramanev) rise_2(r,RiseParamsDetschPure)];
    BubTime ={};
    BubRise ={};
    
    % Recalculating time of flow along tunnel length
    TimeTunnel = zeros(1,length(ComponentsMain(2,:)));
    for ii=1:(length(ComponentsMain(2,:))-1)

        dx(ii) = abs(ComponentsMain(1,ii) - ComponentsMain(1,ii+1));
        U(ii)  = QTunnel(ii)/ComponentsMain(2,ii);
        TimeTunnel(ii+1) = TimeTunnel(ii)+ (dx(ii)/U(ii));

    end
    
    TimeTunnel((RiseStart - 1)*NN + 1)
    TimeTunnel((RiseEnd - 1)*NN)
    
%     figure(15)
%     hold all
%     plot(ComponentsMain(1,(RiseStart - 1)*NN + 1:(RiseEnd - 1)*NN),ComponentsMain(2,(RiseStart - 1)*NN + 1:(RiseEnd - 1)*NN))
%     plot(ComponentsMain(1,(RiseStart - 1)*NN + 1:(RiseEnd - 1)*NN),ComponentsMain(5,(RiseStart - 1)*NN + 1:(RiseEnd - 1)*NN))    
%     pause
    
    for i=1:length(Ut)
        opts = odeset('NonNegative',1);     % No negative terms
        [t_r,Y]=ode45(@(t,y) dydt_terminal(t,y,Ut(i)), [TimeTunnel((RiseStart - 1)*NN + 1) TimeTunnel((RiseEnd - 1)*NN)], 0,opts);

        FinalHeight(i) = Y(end);
    end
%     
%     
%     figure(10)
%     bar(FinalHeight)
%     ylabel('Maximum vertical displacement [m]')
% %     set(gca,'yscale','log')
%     set(gca, 'XTick', 1:3, 'XTickLabel', {'Stokes' 'Karamanev' 'Detsch (Pure)'});
% %     set(gca, 'XTick', 1:4, 'XTickLabel', {'Stokes' 'Karamanev' 'Detsch (Sea)' 'Detsch (Pure)'});

    

end



