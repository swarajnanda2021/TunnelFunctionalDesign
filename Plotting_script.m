clear all
close all
clc



%% Turning radius and height



% Plotting script for parametric calculations
load('Test1_Diffmax_RT1_Hmax_1013p25Pa.mat')

figure(1)
[C1,h1] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),TunnelVolumeTest1(:,:,1),5,'Linecolor',[0 1 0]);   
hold on; 
[C2,h2] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),TunnelVolumeTest1(:,:,2),5,'Linecolor',[0 0 1]);
[C8,h8] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),TunnelVolumeTest1(:,:,end),5,'Linecolor',[0 0.5 1]);
[C3,h3] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),(LossMainTest1(:,:,1)./(9.81*1000)).^(1),3,'Linecolor',[0 1 0]);
[C4,h4] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),(LossMainTest1(:,:,2)./(9.81*1000)).^(1),3,'Linecolor',[0 0 1]);
[C7,h7] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),(LossMainTest1(:,:,end)./(9.81*1000)).^(1),3,'Linecolor',[0 0.5 1]);
[C5,h5] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),LTop(:,:,2),[10.2 10.5 10.9],'Linecolor',[0 0 0]);
[C6,h6] = contour(HeightMaxTest1(:,:,1),TurnRadiusTest1(:,:,1),NPSH_aTest1(:,:,2),[2 2.5 3 3.5 4],'Linecolor',[0 0 0]);
hold off
clabel(C1,h1)
clabel(C2,h2)
clabel(C3,h3)
clabel(C4,h4)
clabel(C5,h5)
clabel(C6,h6)
clabel(C7,h7)
clabel(C8,h8)
set(h1,'linestyle','-');
set(h2,'linestyle','-');
set(h3,'linestyle','-.');
set(h4,'linestyle','-.');
set(h5,'linestyle','-.');
set(h6,'linestyle','-');
set(h7,'linestyle','-.');
set(h8,'linestyle','-');
xlabel('Tunnel centreline maximum height (m)')
ylabel('Main-bend turning radius (m)')
legend('TunnelVolume (wo/main diffuser)','TunnelVolume (w/main diffuser D_{out} = 0.7m)','TunnelVolume (w/main diffuser D_{out} = 0.8m)','Loss [m] (wo/main diffuser)','Loss [m] (w/main diffuser D_{out} = 0.7m)','Loss [m] (w/main diffuser D_{out} = 0.8m)','Top length','NPSH available (w/main diffuser D_{out} = 0.7m)')
grid on

%% Turbulence parameters: Settling chamber length and contraction length

% Part 1: Decision regarding settling chamber length and contraction length
% Plotting script for parametric calculations
load('Test2_Lsettling_Lcontraction_TurbRed.mat')

figure(2)
[C1,h1] = contour(ContractionLengthTest1,SettlingLengthTest1,TunnelVolumeTest1,5,'Linecolor',[0 1 0]);   
hold on; 
[C2,h2] = contour(ContractionLengthTest1,SettlingLengthTest1,TIEvolMainTest1,5,'Linecolor',[0 0 1]);
[C4,h4] = contour(ContractionLengthTest1,SettlingLengthTest1,LTop,5,'Linecolor',[0 0 1]);
hold off
clabel(C1,h1)
clabel(C2,h2)
clabel(C4,h4)
set(h1,'linestyle','-');
set(h2,'linestyle','-');
set(h4,'linestyle','-.');
xlabel('Contraction length (m)')
ylabel('Settling chamber length (m)')
legend('Tunnel Volume','Turbulence Reduction ratio','Top Length')
grid on

%%
% Part 2: Decision regarding Honeycomb parameters

load('Test3_HoneycombCellWidth_HoneycombCellWidth_TurbRed.mat')

figure(3)
[C1,h1] = contour(MeshWidthHoneycombTest1,CellLengthHoneycombTest1,LossMainTest1./(9.81*1000),5,'Linecolor',[0 1 0]);   
hold on; 
[C2,h2] = contour(MeshWidthHoneycombTest1,CellLengthHoneycombTest1,TIEvolMainTest1,5,'Linecolor',[0 0 1]);
[C3,h3] = contour(MeshWidthHoneycombTest1,CellLengthHoneycombTest1,CellLengthHoneycombTest1./MeshWidthHoneycombTest1,[25 35 40],'Linecolor',[0 0 1]);
hold off
clabel(C1,h1)
clabel(C2,h2)
clabel(C3,h3)
set(h1,'linestyle','-');
set(h2,'linestyle','-');
set(h3,'linestyle','-.');
xlabel('Cell Width (m)')
ylabel('Cell length (m)')
legend('Tunnel Loss (m)','Turbulence reduction ratio','Honeycomb L/D')
grid on




%% Separator considerations
% 
% % Length and volume consideration
% load('Test4_SepVolume_length_widthsec.mat')
% 
% figure(4)
% [C1,h1] = contour(LengthSepSecTest1,WidthSepSecTest1,SepSecVelTest1,[0.8 0.85 0.9 0.95 1],'Linecolor',[0 1 0]);   
% hold on; 
% [C2,h2] = contour(LengthSepSecTest1,WidthSepSecTest1,TunnelVolumeTest1,[13 13.6 14.2 14.8 15.3 16.1 16.7 17.9 18.4],'Linecolor',[0 0 1]);
% [C3,h3] = contour(LengthSepSecTest1,WidthSepSecTest1,LossMainTest1./(9.81*1000),[1.3 1.31 1.32 1.33 1.34],'Linecolor',[1 0 1]);
% hold off
% clabel(C1,h1)
% clabel(C2,h2)
% clabel(C3,h3)
% set(h1,'linestyle','-');
% set(h2,'linestyle','-');
% set(h3,'linestyle','-.');
% xlabel('Separator length (m)')
% ylabel('Section width (m)')
% legend('Separator velocity (m/sec)','Tunnel volume (m^3)','Tunnel loss (m)')
% grid on


%% Necessity of diffuser after pump
load('Test5_SepVolume_diain_widthsec.mat')

figure(5)
[C1,h1] = contour(SepInTest1,WidthSepSecTest1,SepSecVelTest1,5,'Linecolor',[0 1 0]);   
hold on; 
[C2,h2] = contour(SepInTest1,WidthSepSecTest1,TunnelVolumeTest1,5,'Linecolor',[0 0 1]);
[C3,h3] = contour(SepInTest1,WidthSepSecTest1,LossMainTest1./(9.81*1000),5,'Linecolor',[1 0 1]);
hold off
clabel(C1,h1)
clabel(C2,h2)
clabel(C3,h3)
set(h1,'linestyle','-');
set(h2,'linestyle','-');
set(h3,'linestyle','-.');
xlabel('Separator inlet dia (m)')
ylabel('Section width (m)')
legend('Separator velocity (m/sec)','Tunnel volume (m^3)','Tunnel loss (m)')
grid on




%%

%% Necessity of diffuser after pump
load('Test7_Tunnel_params_vs_Qpump_LargeTest.mat')

Losshead_large = LossMainTest1(:,1)/(1000*9.81);
NPSHA_large = NPSH_aTest1(:,1);
ReqPower_large = Losshead_large'.*QPUMP.*(2*(9.81));

load('Test6_Tunnel_params_vs_Qpump_CavTunnel.mat')

Losshead_small = LossMainTest1(:,1)/(1000*9.81);
NPSHA_small = NPSH_aTest1(:,1);
ReqPower_small = Losshead_small'.*QPUMP.*(2*(9.81));

figure(6)
grid on
hold on
yyaxis left
plot(QPUMP,Losshead_large,'Color','b')
plot(QPUMP,NPSHA_large,'Color','b')
plot(QPUMP,Losshead_small,'Color','r')
plot(QPUMP,NPSHA_small,'Color','r')
ylabel('Head (m)')
refline(0,5.5)
yyaxis right
plot(QPUMP,ReqPower_large,'Color','k')
plot(QPUMP,ReqPower_small,'Color','k')
xlabel('Pump discharge (m^3/sec)')
ylabel('Required power (kW)')
legend('Loss head (large test-section mode)','NPSH available','Loss head (cavitation test-section mode)','NPSH available','Pumping line datum','Required power (large test-section) (kW)','Required power (cavitation) (kW)')
title('Pump curves')
hold off














