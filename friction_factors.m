% function [ fric_facs, Q_inlet_comp ] = friction_factors( Q, Params, ComponentsMain, NN, nu, Q_pump, Q_pipe, TunnelSolidity, ValveDia, roughness )
function [ fric_facs, Q_inlet_comp ] = friction_factors( Q, Params, ComponentsMain, NN, nu, Q_pump,Q_pipe,roughness )
%friction_factors Calculates friction factor for each component in Params
%and outputs an array of friction factors. Also estimates the turbulence
%reduction factor.


n_comp      = size(Params,2);
fric_facs   = zeros(size(n_comp));
counter     = 1;

% Assign lengths array
lengths = zeros(size(Params,2));
for k = 1:size(Params,2) 
    lengths(k) = str2double(cellstr(Params{3,k}));
end
k_friction = roughness;%0.01e-3;


for j=1:n_comp
    
    % Input component categorical information
    type_comp = cellstr(Params{1,j});
    type_cross= cellstr(Params{2,j});
    % Find out the inlet flow-rate
    Q_comp    = Q((counter-1)*NN + 1);
    Q_inlet_comp(j) = Q_comp;
    
    
    if strcmp(type_comp,'Separator') == 1 
        if strcmp(type_cross,'ducted') == 1
            disp('Separator bank')
            
            % Input data
            channel_height  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            channel_width   = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
            nBanks      = str2double(cellstr(Params{10,j}));
            nPlates     = str2double(cellstr(Params{11,j}));
            
            total_channels = nBanks*nPlates;
            
            h = channel_height;
            b = channel_width;
            D_eff = ((2 * b*h) / (b + h)) * 0.64;            

            w = (Q_comp/total_channels)/((h/nPlates)*(b/nBanks));
            Re = w*D_eff/nu;

            % Calculation of friction facgtor from Haaland's explicit
            % formula

            f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);

            fric_facs(j) = f*lengths(j)/D_eff;
            
        elseif strcmp(type_cross,'Top') == 1
            disp('Separator top section')
            
            % Input data
            channel_height  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_separator   = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
            % Area
            R               = dia_separator/2;
            alpha_subtended = 2*pi - (2*(acos((channel_height - R)/R)));
            Area_channel    = R^2/2 * abs((alpha_subtended)- sin(alpha_subtended));
            
            % Hydraulic Dia
            D_h = sqrt(Area_channel*4/pi);
            
            w = Q_comp/(pi/4 * D_h^2);
            Re = w*D_h / nu;
                 
            f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
            
            fric_facs(j) = f*lengths(j)/D_h;
            
        elseif strcmp(type_cross,'Bottom') == 1
            
            disp('Separator bottom section')
            
            % Input data
            channel_height  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_separator   = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
            % Area
            R               = dia_separator/2;
            alpha_subtended = 2*(acos((R-channel_height)/R));
            Area_channel    = R^2/2 * abs((alpha_subtended)- sin(alpha_subtended));
            
            % Hydraulic Dia
            D_h = sqrt(Area_channel*4/pi);
                
            w = Q_comp/(pi/4 * D_h^2);
            Re = w*D_h / nu;
                 
            f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
            
            fric_facs(j) = f*lengths(j)/D_h;
            
        elseif strcmp(type_cross,'ExpandingBend') == 1
            
            disp('Separator expanding bend')
            
            % Input data
            height_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_separator   = str2double(cellstr(Params{5,j}));                                      % Cross-section based information
            height_out = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
            % Area at inlet
            % Area
            R               = dia_separator/2;
            alpha_subtended1 = 2*(acos((R - height_in)/R));
            Area_channel1    = R^2/2 * abs((alpha_subtended1)- sin(alpha_subtended1));
            
            % Hydraulic Dia
            D_h1 = sqrt(Area_channel1*4/pi);

            % Area at outlet
            % Area
            R               = dia_separator/2;
            alpha_subtended2 = 2*pi - (2*(acos((height_out - R)/R)));
            Area_channel2    = R^2/2 * abs((alpha_subtended2)- sin(alpha_subtended2));
            
            % Hydraulic Dia
            D_h2 = sqrt(Area_channel2*4/pi);
            
            % Turning parameters
            D_h = 0.5*(D_h1 + D_h2); % Average of top and bottom hydraulic dia is taken
            R0  = 1.1*(0.5*D_h);     % 1.1 time the average hydraulic radius is taken as the turning dia
            w_in = Q_comp/(pi/4 * D_h^2); % Reynolds number is also calculated using these average params
            Re_in = w_in*D_h / nu;
            
            % Loss due to 180 degree bend, incorporated within the
            % separator
            A1 = 0.7 + (0.35*180/90);
            B1 = 0.21*((R0/D_h)^(-2.5));
            C1 = 1;
            zeta_loc_bend = A1*B1*C1;
            % Calculation of f
            f = (-1.8 * log10(6.9/Re_in + ((k_friction/D_h)/3.7)^1.11))^(-2);
            zeta_fr_bend  = 0.0175 * 180 * f * (R0/D_h);
            
            % Loss due to the expansion
            expansion = (1 - ((pi/4 * D_h1^2)/(pi/4 *D_h2^2)))^2;
            
            % Summary of frictions
            fric_facs(j) = zeta_loc_bend + zeta_fr_bend + expansion;
                
            
        elseif strcmp(type_cross,'ContractingBend') == 1
            
            disp('Separator contracting bend')
            
            
            % Input data
            height_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_separator   = str2double(cellstr(Params{5,j}));                                      % Cross-section based information
            height_out = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
           % Area at inlet
            % Area
            R               = dia_separator/2;
            alpha_subtended1 = 2*pi - (2*(acos((height_in - R)/R)));
            Area_channel1    = R^2/2 * abs((alpha_subtended1)- sin(alpha_subtended1));
            
            % Hydraulic Dia
            D_h1 = sqrt(Area_channel1*4/pi);
            
            % Area at outlet
            % Area
            R               = dia_separator/2;
            alpha_subtended2 = 2*(acos((R - height_out)/R));
            Area_channel2    = R^2/2 * abs((alpha_subtended2)- sin(alpha_subtended2));
            
            % Hydraulic Dia
            D_h2 = sqrt(Area_channel2*4/pi);
            
            % Turning parameters
            D_h = 0.5*(D_h1 + D_h2); % Average of top and bottom hydraulic dia is taken
            R0  = 1.1*(0.5*D_h);     % 1.1 time the average hydraulic radius is taken as the turning dia
            w_in = Q_comp/(pi/4 * D_h^2); % Reynolds number is also calculated using these average params
            Re_in = w_in*D_h / nu;
            
            % Loss due to 180 degree bend, incorporated within the
            % separator
            A1 = 0.7 + (0.35*180/90);
            B1 = 0.21*((R0/D_h)^(-2.5));
            C1 = 1;
            zeta_loc_bend = A1*B1*C1;
            % Calculation of f
            f = (-1.8 * log10(6.9/Re_in + ((k_friction/D_h)/3.7)^1.11))^(-2);
            zeta_fr_bend  = 0.0175 * 180 * f * (R0/D_h);
            
            % Loss due to the contraction
            contraction = 0.5*(1 - ((pi/4 * D_h2^2)/(pi/4 *D_h1^2)))^0.75;
            
            % Summary of frictions
            fric_facs(j) = zeta_loc_bend + zeta_fr_bend + contraction;
            
                
            
        elseif strcmp(type_cross,'SuddenExpansion') == 1
            % Remember, you always have an inlet that is sudden expanding
            % and an outlet that is sudden contracting. This means that the
            % separator always has the largest dia in the system.
            
            % For this component, the geometry params only show a bend
            % pathway, i.e., constant cross section of length r*theta. The
            % expansion and contraction is incorporated into the friction
            % calculations
            
            disp('Separator sudden expansion')
            
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            
            dia_separator   = str2double(cellstr(Params{5,j}));                                   % Cross-section based information
            height_out = str2double(cellstr(Params{6,j}));                                        % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
            
            % Area at outlet
            % Area
            % If the separator radius is greater than the channel height
            if height_out > (dia_separator/2)
                R               = dia_separator/2;
                alpha_subtended2 = 2*pi - (2*(acos((height_out - R)/R)));
                Area_channel2    = R^2/2 * ((alpha_subtended2)- sin(alpha_subtended2));
            else % If the separator radius is smaller than the channel height
                R               = dia_separator/2;
                alpha_subtended2 = 2*(acos((R - height_out)/R));
                Area_channel2    = R^2/2 * ((alpha_subtended2)- sin(alpha_subtended2));
            end
            
            % Hydraulic Dia
            D_h = sqrt(Area_channel2*4/pi);
            
            % Friction modelling of this component is simply that of
            % expansion and rounded 90 degree bend always
            
            %Diagram 6-1 of Idelchik
            %A1
            A1 = 1; % due to 90 degree bend
            %B1
            if ((R0/dia_in) <= 1 && (R0/dia_in) >= 0.5) 
                B1 = 0.21*((R0/dia)^(-2.5));
            elseif (R0/dia) > 1
                B1 = 0.21 * sqrt((R0/dia)^(-0.5));
            end
            C1 = 1;
            w = Q_comp/(pi/4 * dia_in^2);
            Re = w * dia_in/nu;
            lambda = (-1.8 * log10(6.9/Re + ((k_friction/dia_in)/3.7)^1.11))^(-2);
            zeta_fr  = 0.0175 * 90 * lambda * R0/dia_in;
            zeta_loc = A1*B1*C1;
            
            % Loss due to the expansion
            expansion = (1 - ((pi/4 * dia_in^2)/(pi/4 *D_h^2)))^2;
            
            
            zeta = zeta_loc + zeta_fr+ expansion;
            fric_facs(j) = zeta;
            
            
        elseif strcmp(type_cross,'SuddenContraction') == 1
            % Remember, you always have an inlet that is sudden expanding
            % and an outlet that is sudden contracting. This means that the
            % separator always has the largest dia in the system. 
            
            % For this component, the geometry params only show a bend
            % pathway, i.e., constant cross section of length r*theta. The
            % expansion and contraction is incorporated into the friction
            % calculations
            
            disp('Separator sudden contraction')
            
            % Input data
            height_in = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out  = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            dia_separator   = str2double(cellstr(Params{5,j}));                                      % Cross-section based information
            
            
            if height_out > (dia_separator/2)
                R               = dia_separator/2;
                alpha_subtended2 = 2*pi - (2*(acos((height_in - R)/R)));
                Area_channel2    = R^2/2 * ((alpha_subtended2)- sin(alpha_subtended2));
            else % If the separator radius is smaller than the channel height
                R               = dia_separator/2;
                alpha_subtended2 = 2*(acos((R - height_in)/R));
                Area_channel2    = R^2/2 * ((alpha_subtended2)- sin(alpha_subtended2));
            end
            
            % Hydraulic Dia
            D_h = sqrt(Area_channel2*4/pi);
                       
            
            % Friction modelling of this component is simply that of
            % expansion and rounded 90 degree bend always
            
            %Diagram 6-1 of Idelchik
            %A1
            A1 = 1; % due to 90 degree bend
            %B1
            if ((R0/dia_in) <= 1 && (R0/dia_in) >= 0.5) 
                B1 = 0.21*((R0/dia)^(-2.5));
            elseif (R0/dia) > 1
                B1 = 0.21 * sqrt((R0/dia)^(-0.5));
            end
            C1 = 1;
            w = Q_comp/(pi/4 * dia_out^2);
            Re = w * dia_out/nu;
            lambda = (-1.8 * log10(6.9/Re + ((k_friction/dia_out)/3.7)^1.11))^(-2);
            zeta_fr  = 0.0175 * 90 * lambda * R0/dia_out;
            zeta_loc = A1*B1*C1;
            
            % Loss due to the contraction
            contraction = 0.5*(1 - ((pi/4 * dia_out^2)/(pi/4 *D_h^2)))^0.75;
            
            zeta = zeta_loc + zeta_fr+ contraction;
            fric_facs(j) = zeta;
            
        end
    end
    
    
    
    
    if strcmp(type_comp,'Straight') == 1 
        
        if strcmp(type_cross,'pipe') == 1
            disp('Straight pipe')
    
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            if dia_in==dia_out
                
                % Parameters for the calculation
            
                D_h = dia_in;
                w = Q_comp/(pi/4 * D_h^2);
                Re = w*D_h / nu;
                 
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                fric_facs(j) = f*lengths(j)/D_h;
                
            elseif dia_in < dia_out
                % Formula from 4th edition of Idelchik
                D_h = dia_in;
                w = Q_comp/(pi/4 * D_h^2);
                Re = w*D_h / nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
                n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
                friction = (f/(8*sin(Alphaangleby2))) * (1 - (n_ar1)^-2); 
                
                expansion = 3.2 * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                
                fric_facs(j) = friction + expansion;
                
            elseif dia_in > dia_out
                D_h = dia_in;
                w = Q_comp/(pi/4 * D_h^2);
                Re = w*D_h / nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
                n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
                friction = (f/(8*sin(Alphaangleby2))) * (1 - (n_ar1)^-2); 
                
                n0 = n_ar1^-1;
                Alphaanglep = 0.01745 * 2 * Alphaangleby2;
                contraction = ( (-0.0125 * n0^4) + (0.0224 * n0^3) - (0.00723 * n0^2) + (0.00444 * n0) - 0.00745) * (Alphaanglep^3 - (2*pi*Alphaanglep^2) - (10*Alphaanglep));
                fric_facs(j) = contraction + friction;
            end
                
                
%                 D_h = dia_in;
%                 w = Q_comp/(pi/4 * D_h^2);
%                 Re = w*D_h/nu;
%                 n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
%                 Alphaangle = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
%                 friction = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2) * (lengths(j)/D_h);
%                 expansion= 3.2 * tan(2*Alphaangle) * tan(2*Alphaangle)^(1/4) * (1 - (1/n_ar1))^2;
%                 fric_facs(j) = (friction + expansion)*1.2; % 1.2 is a non-uniformity correction
%                  
%             end
%                 % A generalised formula found in Idelchik is used here
%                 % (para 38 39 of Idelchik, ch-5, ed-3)
%                 D_h = dia_in;
%                 Alphaangle = 2* rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
%                 x_bar = lengths(j)/D_h;
%                 x_tilde = log(1 + 2*x_bar * tan(deg2rad(Alphaangle/2)))/(2*tan(deg2rad(Alphaangle/2)));
%                 n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
%                 
%                 l0_bar = str2double(cellstr(Params{14,j}));
%                 
%                 l0_bar = l0_bar/dia_in;
%                 Re_0 = dia_in * (Q_comp/(pi * dia_in^2/ 4))/nu;
%                
%                 % Using look-up table for Phi value
%                 if (Re_0 < 600000 && Re_0 > 50000)
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0);
%                 elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 600000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0_new);
%                 elseif Re_0 < 0.5e5            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 50000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0_new);
%                 end
%                 zeta_un = Phi * (1 - (1/n_ar1))^1.92;
%                  
%                 a = 0.924 / (1 + (1.3 * 1e-5 * Alphaangle^pi));
%                 b = (0.3 + (1.55  * 1.1^(-Alphaangle)))/((1 + (1.03 * 1e-8 * l0_bar^7.5)));
%                 c = 1.05 / (1 + (2.3 * 1e-62 * Re_0^11));
%                 zeta_non = 0.044 * (0.345*Alphaangle)^a * (1 - (0.2*n_ar1 + 0.8)^(-3.82)) * (0.154 * l0_bar)^b * ((2.31 * 1e-6 * Re_0) + (1.1 * 10e-30 * Re_0^5.62))^c;
%                 
%                 lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                 zeta_fr = (lambda / (8 * sin(deg2rad(Alphaangle/ 2) ) ) ) * (1 - (1/n_ar1^2) ); 
%                 zeta_frprime = (1 + (0.5/(1.5 ^ x_tilde)))*zeta_fr;
%                 
%                 % Calculation of friction factor
%                 zeta = zeta_frprime + zeta_non + zeta_un ;
%                 % Assigning to fric_facs matrix
%                 fric_facs(j) = zeta;
%             end
%Method2
% Calculating the diffuser angle (theta/2)
%                 Alphaangle = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
%                 n0      = ((min(dia_in, dia_out)/max(dia_in, dia_out)) ^ 2);
%                 D_h = dia_in;
%                 w = Q_comp/(pi/4 * D_h^2);
%                 Re = w*D_h/nu;
% 
%                 if dia_in < dia_out
%                     zeta  = (1 - (n0))^2 * 2.6 * sin(Alphaangle);
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h ;%* 1.5; % 1.5 is a correction factor to account for a non-uniform profile
%                 
%                 elseif dia_in > dia_out
%                     % If the pipe-section is converging (can be expected),
%                     % we apply the formula in diagram 5-23 of Idelchik.
%                     
% %                     Alphaangle   = rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
% %                     Alphaangle_r = 0.01745*Alphaangle;
% %                     zeta_converging    = (((-0.0125*n0^4) + (0.0224*n0^3) + (-0.00723*n0^2) + (0.00444*n0) + (-0.00745)) * ((Alphaangle_r^3) - (10*Alphaangle_r) - (2*pi * Alphaangle_r^2)))  + zeta;
%                     if Alphaangle <= 45
%                         zeta = 0.5*(1 - (n0))^0.75 * 1.6 * sin(Alphaangle);
%                     elseif Alphaangle > 45
%                         zeta = 0.5*(1 - (n0))^0.75 * sqrt( sin(Alphaangle));
%                     end
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h ;%* 1.5; % 1.5 is a correction factor to account for a non-uniform profile
%                 end
%                 
%             end
 
        elseif strcmp(type_cross,'flexible') == 1
            disp('Flexible hose')
    
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            % Parameters for the calculation
            
            D_h = dia_in;
            w = Q_comp/(pi/4 * D_h^2);
            Re = w*D_h / nu;
                 
            f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
            fric_facs(j) = 3*f*lengths(j)/D_h; % Flexible hoses are assume to have 3 times the losses of normal hoses
        

            
                        
        elseif strcmp(type_cross,'ducted') == 1
            
            disp('Straight duct')
    
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            if (height_in==height_out && width_in==width_out)
                % Friction calculations (fully turbulent case with smooth
                % walls, Idelchik)
            
                % Parameters for the calculation
                h = height_in;
                b = width_in;
                D_eff = ((2 * b*h) / (b + h)) * 0.64;            
            
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;
                        
                % Calculation of friction facgtor from Haaland's explicit
                % formula
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                
                fric_facs(j) = f*lengths(j)/D_eff;
                
            elseif (height_in < height_out) && (width_in < width_out)     % Axisymmetric diffuser
                
                h = height_in;
                b = width_in;
                D_eff = (2 * b*h) / (b + h) * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
                betaby2 = atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j));
                n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(betaby2)^(-1) + sin(Alphaangleby2)^(-1)); 
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                if (2*rad2deg(betaby2) > 4 && 2*rad2deg(betaby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(betaby2)); 
                    expansionbeta = 3.2 * K * (tan(betaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(betaby2) > 12 && 2*rad2deg(betaby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(betaby2));
                    expansionbeta = 3.2 * K * (tan(betaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionbeta = 0;
                end
                
                expansion = (expansionAlphaangle + expansionbeta) / 2;
                fric_facs(j) = friction + expansion;
            
            elseif (height_in == height_out && width_in < width_out) || (height_in < height_out && width_in == width_out)   % Plane diffusers
                
                h = height_in;
                b = width_in;
                D_eff = (2 * b*h) / (b + h) * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
                betaby2 = atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j));
                if rad2deg(Alphaangleby2) > rad2deg(betaby2)
                    n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                    friction = (f/(4*sin(Alphaangleby2))) * ( (b/h * (1 - (n_ar1)^-1)) + (0.5*(1 - (n_ar1)^-2)) ); 
                else
                    n_ar1 = (max(width_out,width_in) / min(width_out,width_in))^2;
                    friction = (f/(4*sin(betaby2))) * ((b/h * (1 - (n_ar1)^-1)) + (0.5*(1 - (n_ar1)^-2))); 
                end
                
                
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 2 - (0.03*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 20)
                    K = 2 - (0.04*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                if (2*rad2deg(betaby2) > 4 && 2*rad2deg(betaby2) < 12)
                    K = 2 - (0.03*2*rad2deg(betaby2)); 
                    expansionbeta = 3.2 * K * (tan(betaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(betaby2) > 12 && 2*rad2deg(betaby2) < 20)
                    K = 2 - (0.04*2*rad2deg(betaby2));
                    expansionbeta = 3.2 * K * (tan(betaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else % Assume negligible expansion losses
                    expansionbeta = 0;
                end
                
                
                
              
                
                expansion = (expansionAlphaangle + expansionbeta) / 2;
                
                fric_facs(j) = friction + expansion;
            
            elseif (height_in == height_out && width_in > width_out) || (height_in > height_out && width_in == width_out)   % Plane contractions
                % Treated as normal contractions
                
                disp('Ducted plane-contraction')
                % Input data
                height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
                height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
                width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
                width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information

                % Friction calculations, based on axisymmetric contraction with
                % the same expansion angle in both directions, smooth walls and
                % fully turbulent case. This is unlike the flat top desired.

                % Parameters for the calculation
                h = height_in;
                b = width_in;
                D_eff = (2 * b*h) / (b + h) *0.64;            
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;

                Re_sh = (Re/4) * (1 + (b/h));
                Re_long = (Re/4) * ((1 + (b/h))/(b/h));
                lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
                lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
                lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
                Alphaangle = atan(((width_out - width_in) / 2) / lengths(j));
                % Calculation of friction factor
                f = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((height_in/height_out) * (width_in/width_out)));

                fric_facs(j) = f;

                
            end
            

% Method 3) Idelchik's own formula
% 
%                 D_h = (2 * b*h) / (b + h);            
%                 w = Q_comp/(height_in*width_in);
%                 Re = w*D_h/nu;
%                 Alphaangle = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
%                 friction = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2) * (lengths(j)/D_h);
%                 expansion= 3.2 * tan(2*Alphaangle) * tan(2*Alphaangle)^(1/4) * (1 - (1/n_ar1))^2;
%                 fric_facs(j) = (friction + expansion)*1.2; % 1.2 is a non-uniformity correction
                
                
% Method 1) A generalised formula found in Idelchik is used here
% (para 40 41 of Idelchik, ch-5, ed-3)
%                 D_h = (2 * height_in * width_in)/(height_in + width_in);
%                 
%                 Alphaangle = 2*rad2deg(atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j)));
%                 beta = 2*rad2deg(atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j)));
%                 
%                 a0 = height_in/D_h;
%                 b0 = width_in/D_h;
%                 x_bar = lengths(j)/D_h;
%                 x_tilde = (1/(4*tan(deg2rad(Alphaangle/2))))*log(((4 * x_bar^2 * tan(deg2rad(Alphaangle/2))^2) + (2 * x_bar * (a0 + b0) * tan(deg2rad(Alphaangle/2))) + (a0*b0))/(a0*b0));
%                 
%                 n_ar1 = max(height_out*width_out,height_in*width_in) / min(height_in*width_in,height_out*width_out);
%                 
%                 l0_bar = str2double(cellstr(Params{14,j})) / D_h;
%                                 
%                 Re_0 = D_h * (Q_comp/(pi/4 * D_h^2))/nu;
%                 
%                 
%                 % Using look-up table for Phi value
%                 if (Re_0 < 600000  && Re_0 > 50000)
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0);
%                 elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 600000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0_new);
%                 elseif Re_0 < 0.5e5            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 0.5e5;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(Alphaangle,Re_0_new);
%                 end
%                 zeta_un = Phi * (1 - (1/n_ar1))^1.76;
%                  
%                 s = 1.06 / (1 + (2.82 * 1e-3 * Alphaangle^2.24));
%                 t = 0.73 / (1 + (4.31 * 1e-6 * l0_bar^7.31));
%                 u = 1.06 / (1 + (1.1 * 1e-30 * Re_0^5.62));
%                 zeta_non = 0.024 * (0.625*Alphaangle)^s * (1 - (2.81*n_ar1 - 1.81)^(-1.04)) * (0.303 * l0_bar)^t * ((4.8 * 1e-7 * Re_0) + 1.8)^u;
%                 
%                 lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                 zeta_fr = (lambda / 16 ) * (1 - (1 / n_ar1^2) ) * ((1/sin(deg2rad(Alphaangle/2))) + (1/sin(deg2rad(beta/2))) );
%                 zeta_frprime = (1 + (0.5/(1.5 ^ x_tilde)))*zeta_fr;
%                 
%                 % Calculation of friction factor
%                 zeta = zeta_frprime + zeta_non + zeta_un ;
%                 % Assigning to fric_facs matrix
%                 fric_facs(j) = zeta;
% 
% %                 % Method 2
%                 h = height_in;
%                 b = width_in;
%                 D_h = (2 * b*h) / (b + h);            
%             
%                 w = Q_comp/(pi/4 * D_h^2);
%                 Re = w*D_h/nu;
%                 if height_in < height_out
%                     zeta  = (1 - (n0))^2 * 2.6 * sin(Alphaangle);    
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h;
%                 
%                 elseif height_in > height_out
%                     % If the pipe-section is converging (can be expected),
%                     % we apply the formula in diagram 5-23 of Idelchik.
%                     
% %                     Alphaangle   = rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
% %                     Alphaangle_r = 0.01745*Alphaangle;
% %                     zeta_converging    = (((-0.0125*n0^4) + (0.0224*n0^3) + (-0.00723*n0^2) + (0.00444*n0) + (-0.00745)) * ((Alphaangle_r^3) - (10*Alphaangle_r) - (2*pi * Alphaangle_r^2)))  + zeta;
%                     if Alphaangle <= 45
%                         zeta = 0.5*(1 - (n0))^0.75 * 1.6 * sin(Alphaangle);
%                     elseif Alphaangle > 45
%                         zeta = 0.5*(1 - (n0))^0.75 * sqrt( sin(Alphaangle));
%                     end
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h;
%             end
%             end

%             
%         end
        

        elseif strcmp(type_cross,'ductedProtuberance') == 1
            disp('Ducted section with Protuberance')
            
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            l_vent      = str2double(cellstr(Params{11,j}));     % length of the knife
            alpha_opening = str2double(cellstr(Params{10,j}));  % Opening angle (degrees)
            
            protuberance = l_vent*sin(deg2rad(alpha_opening)); % The opening angle
            
            % Friction factor of the case of the straight section
            h = height_in;
            b = width_in;
            D_eff = ((2 * b*h) / (b + h)) * 0.64;            

            w = Q_comp/(height_in*width_in);
            Re = w*D_eff/nu;

            % Calculation of friction facgtor from Haaland's explicit
            % formula

            f1 = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
            
            % Friction of the protuberance (Diagram 10-6 of Idelchik)
            c_x = 2.4;
            tau = 1; % approximately
            S_mbyF0 = (protuberance / str2double(cellstr(Params{4,j})));
            if S_mbyF0 > 0.3
               disp('Reduce opening angle, else formula is invalid') 
            end
            
            w_loc = w/(1 - tau*S_mbyF0);
            
            f2 = c_x * (S_mbyF0) * (w_loc/w)^3;
            
            
            fric_facs(j) = f1*lengths(j)/D_eff    + f2;
            
            
            

        elseif strcmp(type_cross,'dissolutionPipe')
                % Loss due to straight section
                dia_sec  = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
                dia_in   = str2double(cellstr(Params{4,j}));
                dia_out  = str2double(cellstr(Params{5,j}));
                      
                XX = (dia_sec*0.5) - (dia_out + 0.015 + 0.03);
                theta = acos(XX/(dia_sec*0.5));
                Abottom = (dia_sec*0.5)^2/2 * (theta - sin(theta));
                Atop = (pi/4 * dia_sec^2) - Abottom;
                Pbottom = theta*dia_sec*0.5;
                Ptop   = pi*dia_sec - Pbottom;
            
                D_h_in = 4*Atop/Ptop;
                D_h_out = 4*Abottom/Pbottom;
            
            
                w_in = Q_comp/(pi/4 * D_h_in^2); 
                Re_in = w_in*D_h_in/nu;
                
                w_out = Q_comp/(pi/4 * D_h_out^2); 
                Re_out = w_out*D_h_out/nu;
                
                
                % Calculation of friction factor from Haaland's explicit
                % formula
                
                f_top = (-1.8 * log10(6.9/Re_in + ((k_friction/D_h_in)/3.7)^1.11))^(-2);
                f_bottom = (-1.8 * log10(6.9/Re_out + ((k_friction/D_h_out)/3.7)^1.11))^(-2);
                
                
                friction_top    = f_top*lengths(j)/D_h_in;
                friction_bottom = f_bottom*lengths(j)/D_h_out;
                
                % Loss due to 180 degree bend, incorporated within the
                % component
                A1 = 0.7 + (0.35*180/90);
                B1 = 0.21*((dia_sec*0.25/dia_sec)^(-2.5));
                C1 = 1;
                zeta_loc_bend = A1*B1*C1;
                zeta_fr_bend  = 0.0175 * 180 * f * (dia_sec*0.5)/dia_sec;
                
                
                % Add friction due to sudden expansion and contraction at
                % the exits
                expansion = (1 - ((pi/4 * dia_in^2)/(pi/4 *D_h_in^2)))^2;
                if D_h_out > dia_out
                    contraction = 0.5*(1 - ((pi/4 * dia_out^2)/(pi/4 *D_h_out^2)))^0.75;
                else
                    contraction = (1 - ((pi/4 * dia_out^2)/(pi/4 *D_h_out^2)))^2;
                end
                
                % Account for differing cross-sectional areas in the
                % calculation of the pressure loss
                friction_top    = friction_top * ((pi/4 * dia_in^2)/Atop)^2;
                friction_bottom = friction_bottom * ((pi/4 * dia_in^2)/Abottom)^2;
                contraction     = contraction* ((pi/4 * dia_in^2)/Atop)^2;
                zeta_fr_bend    = zeta_fr_bend*((pi/4 * dia_in^2)/Abottom)^2;
                expansion       = expansion * ((pi/4 * dia_in^2)/Abottom)^2;
                
                
                fric_facs(j) = friction_top + friction_bottom + zeta_loc_bend + zeta_fr_bend + expansion + contraction;
        
        elseif strcmp(type_cross,'dissolutionBox')
        
                % Loss due to straight section
                height_sep  = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
                width_sep  = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
                height_sec = str2double(cellstr(Params{5,j}));
                width_sec  = str2double(cellstr(Params{10,j}));
                dia_in   = str2double(cellstr(Params{4,j}));
                dia_out  = str2double(cellstr(Params{5,j}));        % Here, dia_out is actually taken as square cross-section, :P. So is dia_in
                      
                Abottom = (height_sep*width_sep) - (dia_in*width_sec);
                Atop = (height_sep*width_sep) - Abottom;
                Pbottom = 2*(width_sep + dia_in);
                Ptop   = 2*(width_sep + dia_out);
            
                D_h_in = 4*Atop/Ptop;
                D_h_out = 4*Abottom/Pbottom;
            
                SepVel = Q_comp/(height_sec*width_sec);
                
                w_in = Q_comp/(pi/4 * D_h_in^2); 
                Re_in = w_in*D_h_in/nu;
                
                w_out = Q_comp/(pi/4 * D_h_out^2); 
                Re_out = w_out*D_h_out/nu;
                
                
                % Calculation of friction factor from Haaland's explicit
                % formula
                
                f_top = (-1.8 * log10(6.9/Re_in + ((k_friction/D_h_in)/3.7)^1.11))^(-2);
                f_bottom = (-1.8 * log10(6.9/Re_out + ((k_friction/D_h_out)/3.7)^1.11))^(-2);
                
                
                friction_top    = f_top*lengths(j)/D_h_in;
                friction_bottom = f_bottom*lengths(j)/D_h_out;
                
                % Loss due to 180 degree bend, incorporated within the
                % component
                % Here we refer to Diagram 6.15 of Idelchik assuming a
                % diverging bend with l_el=0 and bel/b0 = dia_out/dia_in of
                % 2, coming out to be 1.2
                
                C1 = 0.97 - 0.13*log(dia_in/width_sep);%From diagram 6.7 of Idelchik, tentatively
                zeta_loc_bend = C1*1.2; % 1.2 comes form item 4 of diagram 6.15 of Idelchik
                zeta_fr_bend  = (-1.8 * log10(6.9/Re_out + ((k_friction/D_h_out)/3.7)^1.11))^(-2);
                
                
                % Add friction due to sudden expansion and contraction at
                % the exits
                expansion = (1 - ((dia_in*width_sec)/(pi/4 *D_h_in^2)))^2;
                if D_h_out > dia_out
                    contraction = 0.5*(1 - ((dia_out^2)/(pi/4 *D_h_out^2)))^0.75;
                else
                    contraction = (1 - ((dia_out^2)/(pi/4 *D_h_out^2)))^2;
                end
                
                % Account for differing cross-sectional areas in the
                % calculation of the pressure loss
                friction_top    = friction_top * ((pi/4 * dia_in^2)/Atop)^2;
                friction_bottom = friction_bottom * ((pi/4 * dia_in^2)/Abottom)^2;
                contraction     = contraction* ((pi/4 * dia_in^2)/Atop)^2;
                zeta_fr_bend    = zeta_fr_bend*((pi/4 * dia_in^2)/Abottom)^2;
                expansion       = expansion * ((pi/4 * dia_in^2)/Abottom)^2;
                
                fric_facs(j) = friction_top + friction_bottom + zeta_loc_bend + zeta_fr_bend + expansion + contraction;
                
        end        
    end
    
    if strcmp(type_comp,'ContractionSingle') == 1
        % Noting that we use 5th order polynomials for the contraction, we
        % expect no separation such that the only loss component is by
        % friction alone. Chapters 2 and 5 of Idelchik are used here.
        
        if strcmp(type_cross,'ducted') == 1
            disp('Ducted contraction')
            % Input data
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            % Friction calculations, based on axisymmetric contraction with
            % the same expansion angle in both directions, smooth walls and
            % fully turbulent case. This is unlike the flat top desired.
            if height_in > height_out
                % Parameters for the calculation
                h = height_in;
                b = width_in;
                D_eff = (2 * b*h) / (b + h) *0.64;            
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;

                Re_sh = (Re/4) * (1 + (b/h));
                Re_long = (Re/4) * ((1 + (b/h))/(b/h));
                lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
                lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
                lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
                Alphaangle = atan(((height_out - height_in) / 2) / lengths(j));
                % Calculation of friction factor
                f = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((height_in/height_out) * (width_in/width_out)));

                fric_facs(j) = f;
            elseif height_in == height_out
               % Parameters for the calculation
                h = height_in;
                b = width_in;
                D_eff = (2 * b*h) / (b + h) *0.64;            
                w = Q_comp/(height_in*width_in);
                Re = w*D_eff/nu;

                Re_sh = (Re/4) * (1 + (b/h));
                Re_long = (Re/4) * ((1 + (b/h))/(b/h));
                lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
                lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
                lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
                Alphaangle = atan(((width_out - width_in) / 2) / lengths(j));
                % Calculation of friction factor
                f = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((height_in/height_out) * (width_in/width_out)));

                fric_facs(j) = f; 
                
                
                
            end
                        
        elseif strcmp(type_cross,'pipe') == 1    
            
            disp('Conical contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j}));
    
            % A generalised formula found in Idelchik is used here
            % (para 38 39 of Idelchik, ch-5, ed-3)
            D_h = dia_in;
            Alphaangle = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)); 
            Re_0 = dia_in * (Q_comp/(pi * dia_in^2/ 4))/nu;
               
            % Using look-up table for Phi value
            if Re_0 < 600000
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0);
            elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
                Re_0_new = 600000;
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0_new);
            end
                
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            zeta_fr = (lambda / (8 * sin(Alphaangle / 2) ) ) * (1 - (1 / ( (max(dia_in, dia_out)/min(dia_in, dia_out)) ^ 2) ) ); 
            
            % Calculation of friction factor
            zeta = zeta_fr ;
            % Assigning to fric_facs matrix
            fric_facs(j) = zeta;
            
        elseif strcmp(type_cross,'C2R') == 1
            disp('Circle to Square Contraction')
        
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
             
            % Parameters for the calculation
            h = height_out;
            b = width_out;
            
            D_h = (2 * b*h) / (b + h) ;            
            w = Q_comp/(height_out*width_out);
            
            
            
            if ((pi/4*D_h^2 - pi/4*dia_in^2) < 1e-1) % Converging and straight pieces. Ignore contraction loss due to it being very small.
                D_eff = 0.64*D_h;
                % Calculation of friction factor
                Re_0 = D_eff*w/nu;               
                lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                % Calculation of friction factor
                zeta = lambda * lengths(j)/(D_eff);
                % Assigning to fric_facs matrix
                fric_facs(j) = zeta;

            elseif (pi/4*D_h^2 - pi/4*dia_in^2) > 1e-1 % Diverging piece
            

                D_eff = D_h * 0.64;            
                
                w = Q_comp/(height_out*width_out);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan((dia_in - 2*(sqrt(height_out*width_out/pi)))/2*lengths(j));
               
                n_ar1 = (max(height_out,height_out) / min(height_out,height_out))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(Alphaangleby2)^(-1) + sin(Alphaangleby2)^(-1));
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                
                expansion = (expansionAlphaangle) / 2;
                fric_facs(j) = friction + expansion;
                
            end
                      
            
            
        end
            
    
    end
    
    if strcmp(type_comp,'ContractionDouble') == 1
        % Noting that we use 5th order polynomials for the contraction, we
        % expect no separation such that the only loss component is by
        % friction alone. Chapters 2 and 5 of Idelchik are used here.
        
        if strcmp(type_cross,'ducted') == 1
            disp('Ducted contraction')
            % Input data
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            CR1x         = str2double(cellstr(Params{10,j}));
            CR1y         = str2double(cellstr(Params{11,j}));
            CR1          = (CR1x + CR1y)/2;
            CR2          = str2double(cellstr(Params{12,j}));
            % Length of individual stages
            length1         = str2double(cellstr(Params{13,j}));
            length2         = str2double(cellstr(Params{14,j}));
            
            % Friction calculations, based on axisymmetric contraction with
            % the same expansion angle in both directions, smooth walls and
            % fully turbulent case. This is unlike the flat top desired.
            
            % Parameters for the calculation
            % Stage 1
            h = height_in;
            b = width_in;
            D_eff = (2 * b*h) / (b + h) *0.64;            
            w = Q_comp/(h*b);
            Re = w*D_eff/nu;
                    
            Re_sh = (Re/4) * (1 + (b/h));
            Re_long = (Re/4) * ((1 + (b/h))/(b/h));
            lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
            lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
            lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
            Alphaangle = atan(((sqrt(height_out^2*CR2) - h) / 2) / length1);
            % Calculation of friction factor
            f1 = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((h/sqrt(height_out^2*CR2)) * (b/sqrt(height_out^2*CR2))));
            
            % Stage 2
            h = sqrt(height_out^2*CR2);
            b = sqrt(width_out^2*CR2);
            D_eff = (2 * b*h) / (b + h) *0.64;            
            w = Q_comp/(h*b);
            Re = w*D_eff/nu;
                    
            Re_sh = (Re/4) * (1 + (b/h));
            Re_long = (Re/4) * ((1 + (b/h))/(b/h));
            lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
            lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
            lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
            Alphaangle = atan(((height_out - h) / 2) / length2);
            % Calculation of friction factor
            f2 = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((h/height_out) * (b/width_out)));
            
            fric_facs(j) = f1 + f2;
            
                        
        elseif strcmp(type_cross,'pipe') == 1    
            
            disp('Conical contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j}));
            
            
            CR1x         = str2double(cellstr(Params{10,j}));
            CR1y         = str2double(cellstr(Params{11,j}));
            CR1          = (CR1x + CR1y)/2;
            CR2          = str2double(cellstr(Params{12,j}));
            % Length of individual stages
            length1         = str2double(cellstr(Params{13,j}));
            length2         = str2double(cellstr(Params{14,j}));
            
            % Stage 1
            % A generalised formula found in Idelchik is used here
            % (para 38 39 of Idelchik, ch-5, ed-3)
            CR = (dia_in/dia_out)^2;
            d1 = dia_in;
            d2 = sqrt(dia_in^2/CR);
            D_h = d1;
            
            Alphaangle = atan(((max(d1,d2) - min(d1,d2))/2) / length1); 
            Re_0 = d1 * (Q_comp/(pi * d1^2/ 4))/nu;
               
            % Using look-up table for Phi value
            if Re_0 < 600000
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0);
            elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
                Re_0_new = 600000;
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0_new);
            end
            % Haaland formula    
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            % Idelchik
            zeta_fr1 = (lambda / (8 * sin(Alphaangle / 2) ) ) * (1 - (1 / ( (max(d1, d2)/min(d1, d2)) ^ 2) ) ); 
            
            % Stage 2
            % A generalised formula found in Idelchik is used here
            % (para 38 39 of Idelchik, ch-5, ed-3)
            d1 = sqrt(dia_in^2/CR);
            d2 = dia_out;
            
            D_h = d1;
            
            Alphaangle = atan(((max(d1,d2) - min(d1,d2))/2) / length2); 
            Re_0 = d1 * (Q_comp/(pi * d1^2/ 4))/nu;
               
            % Using look-up table for Phi value
            if Re_0 < 600000
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0);
            elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
                Re_0_new = 600000;
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(Alphaangle),Re_0_new);
            end
            % Haaland formula    
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            % Idelchik
            zeta_fr2 = (lambda / (8 * sin(Alphaangle / 2) ) ) * (1 - (1 / ( (max(d1, d2)/min(d1, d2)) ^ 2) ) );
            
            % Calculation of friction factor
            zeta = zeta_fr1 + zeta_fr2 ;
            % Assigning to fric_facs matrix
            fric_facs(j) = zeta;
            
        elseif strcmp(type_cross,'C2R') == 1
           disp('Circle to Square adapter')
        
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            CR1x         = str2double(cellstr(Params{10,j}));
            CR1y         = str2double(cellstr(Params{11,j}));
            CR1          = (CR1x + CR1y)/2;
            CR2          = str2double(cellstr(Params{12,j}));
            % Length of individual stages
            length1         = str2double(cellstr(Params{13,j}));
            length2         = str2double(cellstr(Params{14,j}));
            
            % Stage 1 : Circle to square contraction
            % Parameters for the calculation
            h = sqrt(height_out^2 *CR2);
            b = sqrt(width_out^2 *CR2);
            
            D_h = (2 * b*h) / (b + h) ;            
            w = Q_comp/(h*b);
            
            D_eff = 0.64*D_h;
            % Calculation of friction factor
            Re_0 = D_eff*w/nu;               
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
            % Calculation of friction factor
            zeta = lambda * lengths1/(D_eff);
            % Assigning to fric_facs matrix
            fric_facs1 = zeta;

            
            
            % Stage 2 : Fully ducted contraction
            % Parameters for the calculation
            h = sqrt(height_out^2*CR2);
            b = sqrt(width_out^2*CR2);
            D_eff = (2 * b*h) / (b + h) *0.64;            
            w = Q_comp/(h*b);
            Re = w*D_eff/nu;
                    
            Re_sh = (Re/4) * (1 + (b/h));
            Re_long = (Re/4) * ((1 + (b/h))/(b/h));
            lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
            lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
            lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
            Alphaangle = atan(((height_out - h) / 2) / length2);
            % Calculation of friction factor
            f2 = (lambda / (8 * abs(sin(Alphaangle/2)))) * (1 - 1/((h/height_out) * (b/width_out)));
            
            fric_facs(j) = fric_facs1 + f2;
        end
            
    
    end
    
    if strcmp(type_comp,'Bend') == 1
        % Cross-section based information, Chapter-6 of Idelchik.
        
        if strcmp(type_cross,'pipe') == 1
            disp('Pipe bend')
            % Input data
            dia         = str2double(cellstr(Params{4,j}));               % Cross-section based information
             
            % Modify lengths(j)
            lengths(j) = pi*dia/4;
            
            corner_mod = Params{12,j};
            
            
            D_h = dia;
            w = Q_comp/(pi/4 * dia^2);
            Re = w * D_h/nu;
            k_re = 0.8 + (4.02e4/Re);
                
            
            if strcmp(corner_mod,'thin_rounded') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.23 * k_re) + (1.28 * lambda);
            elseif strcmp(corner_mod,'thin_bevelled') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.3 * k_re) + (1.28 * lambda);
            elseif strcmp(corner_mod,'profiled_rounded_normal') == 1
                %rounding is r=0.18*/D0
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.23 * k_re) + (1.28 * lambda);
            elseif strcmp(corner_mod,'profiled_rounded_reduced') == 1
                %rounding is r=0.18*/D0
                % Vanes installed according to arithmetic progression
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.15 * k_re) + (1.28 * lambda);
                
            elseif strcmp(corner_mod,'profiled_bevelled_normal') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.30 * k_re) + (1.28 * lambda);
            
            elseif strcmp(corner_mod,'profiled_bevelled_reduced') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.23 * k_re) + (1.28 * lambda);
            
            end
           
            fric_facs(j) = zeta;       
            
        elseif strcmp(type_cross,'ducted') == 1
            disp('Ducted bend')
            % Input data
                        
            height          = str2double(cellstr(Params{4,j}));               % Cross-section based information
            width           = str2double(cellstr(Params{6,j}));    
            
            % Modify lengths(j)
            lengths(j) = height;
            
            corner_mod = Params{12,j};
            
            D_eff = 2*height*width/(height+width) * 0.64;
            w = Q_comp/(height*width);
            Re = w * D_eff/nu;
            k_re = 0.8 + (4.02e4/Re);
            
            if strcmp(corner_mod,'thin_inner_sharp_normal') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                zeta = (0.45 * k_re) +  lambda;
                
            elseif strcmp(corner_mod,'thin_inner_sharp_reduced') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                zeta = (0.45 * k_re) +  lambda;
                            
            elseif strcmp(corner_mod,'thin_inner_bevelled') == 1
                % Bevelling is t = 0.25*width_inlet
                
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                zeta = (0.36 * k_re) + (1.28 * lambda);
                
            elseif strcmp(corner_mod,'thin_rounded') == 1
                % Rounding is r = 0.2*width_inlet
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                zeta_fr = (1 + (1.57*0.2)) * lambda;
                zeta = (k_re * 0.14) + zeta_fr;         
            end

            fric_facs(j) = zeta;
            
            
        elseif strcmp(type_cross,'arbitrary') == 1
            disp('Arbitrary bend')
            dia         = str2double(cellstr(Params{4,j}));
            Alphaangle  = str2double(cellstr(Params{5,j}));                                       % The bending angle in degrees
            R0          = str2double(cellstr(Params{6,j}));
            
            
            %Diagram 6-1 of Idelchik
            %A1
            if Alphaangle <= 70
                A1 = 0.9*sin(deg2rad(Alphaangle));
            elseif Alphaangle == 90
                A1 = 1;
            elseif Alphaangle >= 100
                A1 = 0.7 + (0.35*Alphaangle/90);
            end
            %B1
            
            
            if ((R0/dia) <= 1 && (R0/dia) >= 0.5) 
                B1 = 0.21*((R0/dia)^(-2.5));
            elseif (R0/dia) > 1
                B1 = 0.21 * sqrt((R0/dia)^(-0.5));
            end
            C1 = 1;
            w = Q_comp/(pi/4 * dia^2);
            Re = w * dia/nu;
            lambda = (-1.8 * log10(6.9/Re + ((k_friction/dia)/3.7)^1.11))^(-2);
            zeta_fr  = 0.0175 * Alphaangle * lambda * R0/dia;
            zeta_loc = A1*B1*C1;
            
            zeta = zeta_loc + zeta_fr;
            fric_facs(j) = zeta;          
            
        end
        
    end
    
    if strcmp(type_comp,'Adapter') == 1
        % IMPORTANT NOTICE:
        % For cross-sectional adapters, we simply want to change the
        % cross-section. Therefore,we assume only lengthwise friction
        % factor, similar to a pipe-section or a duct, for simplicity.
        % Transition from a circle to a rectangle
        if strcmp(type_cross,'C2R') == 1
            disp('Circle to Square adapter')
        
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
             
            % Parameters for the calculation
            h = height_out;
            b = width_out;
            
            D_h = (2 * b*h) / (b + h) ;            
            w = Q_comp/(height_out*width_out);
            
            
            
            if ((pi/4*D_h^2 - pi/4*dia_in^2) < 1e-1) % Converging and straight pieces
                D_eff = 0.64*D_h;
                % Calculation of friction factor
                Re_0 = D_eff*w/nu;               
                lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                % Calculation of friction factor
                zeta = lambda * lengths(j)/(D_eff);
                % Assigning to fric_facs matrix
                fric_facs(j) = zeta;

            elseif (pi/4*D_h^2 - pi/4*dia_in^2) > 1e-1 % Diverging piece
            

                D_eff = D_h * 0.64;            
                
                w = Q_comp/(height_out*width_out);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = abs(atan((dia_in - 2*(sqrt(height_out*width_out/pi)))/2*lengths(j)));
               
                n_ar1 = (max(h,b) / min(h,b))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(Alphaangleby2)^(-1) + sin(Alphaangleby2)^(-1));
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                expansion = (expansionAlphaangle) / 2;
                fric_facs(j) = friction + expansion;
                
            end
                      
            
                       
        end
        
        % Transition from a circle to a rectangle
        if strcmp(type_cross,'R2C') == 1
            disp('Square to circle adapter')
        
            % Input data
            dia_out    = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            height_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            width_in   = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            % Parameters for the calculation
            h = height_in;
            b = width_in;
            
            D_h = (2 * b*h) / (b + h) * 0.64;            
            w = Q_comp/(height_in*width_in);
            
            if ((pi/4*D_h^2 - pi/4*dia_in^2) < 1e-1) % Converging and straight pieces
                D_eff = 0.64*D_h;
                % Calculation of friction factor
                Re_0 = D_eff*w/nu;               
                lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                % Calculation of friction factor
                zeta = lambda * lengths(j)/(D_eff);
                % Assigning to fric_facs matrix
                fric_facs(j) = zeta;

            elseif (pi/4*D_h^2 - pi/d*dia_in^2) > 1e-1 % Diverging piece
            

                D_eff = D_h * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan((2*(sqrt(h*b/pi)) - dia_out)/2*lengths(j));
               
                n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(betaby2)^(-1) + sin(Alphaangleby2)^(-1)); 
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                expansion = (expansionAlphaangle) / 2;
                fric_facs(j) = friction + expansion;
                
            end           
        end
        
        
        
    end
    
    
    if strcmp(type_comp,'Mesh') == 1
        
        if strcmp(type_cross,'ducted') == 1
            disp('Mesh in a duct')
            % Input params
            height      = str2double(cellstr(Params{4,j}));
            width       = str2double(cellstr(Params{6,j}));
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
                        
            % Friction calculations
            porosity = (1 - (d/M))^2;
            solidity = 1 - porosity;
            Re_d = d * (Q_comp/(height*width))/nu;
            Re_M = M * (Q_comp/(height*width))/nu;            
            % We restrict to typical values of solidity for which investigations on grid generated turbulence have happened
            % We keep porosity above 0.55 to ensure that merger of jets
            % donot cause in increase in turbulence length scale
            % downstream.
                
            if (porosity > 0.55) 
                disp(['Porosity of ' num2str(porosity) ' greater than 0.55. Increase d/M.'])
            end
            
            
            % Friction factor calculations are taken from the paper "Grid
            % -generated turbulence revisited" by Kurian and Fransson
        
            if (35 < Re_d || Re_d < 800)
                f = 0.5 + (26 * (Re_d^(-1))); 
            elseif Re_d > 800
                f = 0.52 + (66 * (Re_d^(-4/3)));
            elseif Re_d < 35
                f = Re_d^(-1);
            end
        
            fric_facs(j) = f * ((1 - porosity^2)/ porosity^2);
            
            
        elseif strcmp(type_cross,'pipe') == 1
            disp('Mesh in a pipe')
            % Input params
            dia         = str2double(cellstr(Params{4,j}));
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
            
            % Friction calculations
            porosity = (1 - (d/M))^2;
            solidity = 1 - porosity;
            Re_d = d * (Q_comp/(pi/4 * dia^2))/nu;
            Re_M = M * (Q_comp/(pi/4 * dia^2))/nu;            
            if (porosity > 0.55) 
                disp(['Porosity of ' num2str(porosity) ' greater than 0.55. Increase d/M.'])
            end
            
            
            % Friction factor calculations are taken from the paper "Grid
            % -generated turbulence revisited" by Kurian and Fransson
        
            if (35 < Re_d || Re_d < 800)
                f = 0.5 + (26 * (Re_d^(-1))); 
            elseif Re_d > 800
                f = 0.52 + (66 * (Re_d^(-4/3)));
            elseif Re_d < 35
                f = Re_d^(-1);
            end
        
            fric_facs(j) = f * ((1 - porosity^2)/ porosity^2); 
            
    
    
        end
        
        
    end
    
    
    if strcmp(type_comp,'Honeycomb') == 1
        disp('Lumleys honeycomb')
        if strcmp(type_cross,'ducted') == 1
            
            height = str2double(cellstr(Params{4,j}));
            width  = str2double(cellstr(Params{6,j}));
            D_h    = (2 * height*width/(height+width));
        
        
            l = str2double(cellstr(Params{12,j}));
            D = str2double(cellstr(Params{11,j}));
            
            lbyD = l/D;
            Re_d = (Q_comp/(height*width)) * D / nu;
            
            if lbyD < 10
                disp(['l/D of ' num2str(lbyD) ' is lower than 10. Increase.'])
            elseif lbyD > 100
                disp(['l/D of ' num2str(lbyD) ' is greater than 100. Reduce.'])
            end
            if Re_d > 30000
                disp(['Mesh Reynolds number of ' num2str(Re_d) ' is greater than 30000. Reduce.'])
            elseif Re_d < 1000
                disp(['Mesh Reynolds number of ' num2str(Re_d) ' is lower than 1000. Increase.'])
            end
            % Friction calculations (look up table from the iso-reduction
            % curves of Lumley in "Reducing water tunnel turbulence by means of
            % a honeycomb"
            
            data = csvread('Lumley_honeycomb_pressuredrop.txt');
            F = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
            
            if Re_d > 29000
                fric_facs(j) = F(lbyD,29000);
            else
                fric_facs(j) = F(lbyD,Re_d);
            end
            
            
        
        elseif strcmp(type_cross,'pipe')
            % For a honeycomb, we use the data of Lumley's iso-reduction curves
            dia         = str2double(cellstr(Params{4,j}));    
                  
            l = str2double(cellstr(Params{12,j}));
            D = str2double(cellstr(Params{11,j}));
            
            lbyD = l/D;
            Re_d = (Q_comp/(pi/4 * dia^2)) * D / nu;
            
            
            if lbyD < 10
                disp(['l/D of ' num2str(lbyD) ' is lower than 10. Increase.'])
            elseif lbyD > 100
                disp(['l/D of ' num2str(lbyD) ' is greater than 100. Reduce.'])
            end
            if Re_d > 30000
                disp(['Mesh Reynolds number of ' num2str(Re_d) ' is greater than 30000. Reduce.'])
            elseif Re_d < 1000
                disp(['Mesh Reynolds number of ' num2str(Re_d) ' is lower than 1000. Increase.'])
            end
            % Friction calculations (look up table from the iso-reduction
            % curves of Lumley in "Reducing water tunnel turbulence by means of
            % a honeycomb"
        
            data = csvread('Lumley_honeycomb_pressuredrop.txt');
            F = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
            if Re_d > 29000
                fric_facs(j) = F(lbyD,29000);
            else
                fric_facs(j) = F(lbyD,Re_d);
            end
        
            
        end
        
    
    end


    if strcmp(type_comp,'Collector') == 1
        if strcmp(type_cross,'inlet') == 1
           
            length_hose = str2double(cellstr(Params{8,j}));
            dia_hose =  str2double(cellstr(Params{4,j}));
            hose_outlet = str2double(cellstr(Params{5,j}));
            length_contraction = str2double(cellstr(Params{9,j}));
            length_expansion = str2double(cellstr(Params{10,j}));
            % Friction due to shape adapting contraction
        
            % Input data
            dia_out    = dia_hose;                                           % Cross-section based information
            height_in  = 0.3;                                           % Cross-section based information
            width_in   = 0.3;                                           % Cross-section based information
            
            % Parameters for the calculation
            h = height_in;
            b = width_in;
            
            D_h = (2 * b*h) / (b + h) * 0.64;            
            w = Q_comp/(height_in*width_in);
            
            if ((pi/4*D_h^2 - pi/4*dia_out^2) < 1e-1) % Converging and straight pieces
                D_eff = 0.64*D_h;
                % Calculation of friction factor
                Re_0 = D_eff*w/nu;               
                lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                % Calculation of friction factor
                zeta = lambda * length_expansion/(D_eff);
                % Assigning to fric_facs matrix
                zeta_contraction = zeta;

            elseif (pi/4*D_h^2 - pi/d*dia_out^2) > 1e-1 % Diverging piece
            

                D_eff = D_h * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan((2*(sqrt(h*b/pi)) - dia_out)/2*length_expansion);
               
                n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(betaby2)^(-1) + sin(Alphaangleby2)^(-1)); 
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                zeta_contraction = ((expansionAlphaangle) / 2) + friction;
                
            end           
         
            % Friction due to flexible hose (3 times that in smooth pipe)
            
            Re_0 = dia_hose*(Q_comp/(pi/4 * dia_hose^2))/nu;               
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/dia_hose)/3.7)^1.11))^(-2);
            % Calculation of friction factor
            zeta_hose = lambda * length_hose/(dia_hose) * 3; % multiply with 3 to get conservative estimate for a flexible hose
                
            % Friction due to conical expansion
            D_h = dia_hose;
            w = Q_comp/(pi/4 * D_h^2);
            Re = w*D_h / nu;
            
            dia_in = dia_hose;
            dia_out = hose_outlet;
            
            f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
            Alphaangleby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / length_contraction);
            n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
            friction = (f/(8*sin(Alphaangleby2))) * (1 - (n_ar1)^-2); 
                
            expansion = 3.2 * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                
            zeta_expansion = friction + expansion;
            
            % Re-scale the zeta values to 0.3^2
            zeta_contraction = zeta_contraction;
            zeta_hose = zeta_hose*(0.3^2/(pi/4 * dia_hose^2));
            zeta_expansion = zeta_expansion*(0.3^2/(pi/4 * dia_hose^2));
            
            fric_facs(j) = zeta_contraction + zeta_hose + zeta_expansion;
            
        elseif strcmp(type_cross,'outlet') == 1
        
            length_hose = str2double(cellstr(Params{8,j}));
            dia_hose =  str2double(cellstr(Params{4,j}));
            hose_inlet = str2double(cellstr(Params{5,j}));
            length_contraction = str2double(cellstr(Params{9,j}));
            length_expansion = str2double(cellstr(Params{10,j}));
            % Friction due to shape adapting contraction
        
            % Input data
            dia_in    = dia_hose;                                           % Cross-section based information
            height_out  = 0.3;                                           % Cross-section based information
            width_out   = 0.3;                                           % Cross-section based information
            
            % Parameters for the calculation
            h = height_out;
            b = width_out;
            
            D_h = (2 * b*h) / (b + h) * 0.64;            
            w = Q_comp/(height_out*width_out);
            
            if ((pi/4*D_h^2 - pi/4*dia_in^2) < 1e-1) % Converging and straight pieces
                D_eff = 0.64*D_h;
                % Calculation of friction factor
                Re_0 = D_eff*w/nu;               
                lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                % Calculation of friction factor
                zeta = lambda * length_expansion/(D_eff);
                % Assigning to fric_facs matrix
                zeta_expansion = zeta;

            elseif (pi/4*D_h^2 - pi/d*dia_in^2) > 1e-1 % Diverging piece
            

                D_eff = D_h * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_eff)/3.7)^1.11))^(-2);
                Alphaangleby2 = atan((2*(sqrt(h*b/pi)) - dia_out)/2*length_expansion);
               
                n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(betaby2)^(-1) + sin(Alphaangleby2)^(-1)); 
                
                if (2*rad2deg(Alphaangleby2) > 4 && 2*rad2deg(Alphaangleby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(Alphaangleby2)); 
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(Alphaangleby2) > 12 && 2*rad2deg(Alphaangleby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(Alphaangleby2));
                    expansionAlphaangle = 3.2 * K * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionAlphaangle = 0;
                end
                
                
                zeta_expansion = ((expansionAlphaangle) / 2) + friction;
                
            end           
         
            % Friction due to flexible hose (3 times that in smooth pipe)
            
            Re_0 = dia_hose*(Q_comp/(pi/4 * dia_hose^2))/nu;               
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/dia_hose)/3.7)^1.11))^(-2);
            % Calculation of friction factor
            zeta_hose = lambda * length_hose/(dia_hose) * 3; % multiply with 3 to get conservative estimate for a flexible hose
                
            % Friction due to conical contraction
            D_h = dia_hose;
            w = Q_comp/(pi/4 * D_h^2);
            Re = w*D_h / nu;
            
            dia_out = dia_hose;
            dia_in = hose_outlet;
            
            f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
            Alphaangleby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / length_contraction);
            n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
            friction = (f/(8*sin(Alphaangleby2))) * (1 - (n_ar1)^-2); 
                
            contraction = 3.2 * (tan(Alphaangleby2))^1.25 * (1 - (n_ar1)^-1)^2;
                
            zeta_contraction = friction + contraction;
            
            % Re-scale the zeta values to 0.3^2
            zeta_contraction = zeta_contraction;
            zeta_hose = zeta_hose*(hose_inlet^2/dia_hose^2);
            zeta_expansion = zeta_expansion*(hose_inlet^2/dia_hose^2);
            
            fric_facs(j) = zeta_contraction + zeta_hose + zeta_expansion;
            
        end
    end
    
    if strcmp(type_comp,'Lunchbox') == 1
        if strcmp(type_cross,'inlet') == 1
            disp('Intake pipe for pipe section')
            % NOTE: For the intake, the length input is, infact, the angular
            % intake from the main pipe or riser
            % Input data
            intercept       = str2double(cellstr(Params{12,j}));
            dia_inlet       = str2double(cellstr(Params{4,j}));                        
            Alphaangle      = str2double(cellstr(Params{7,j}));               % The divergence angle in degrees
            length_inlet    = str2double(cellstr(Params{3,j}));
            width_lunchbox  = str2double(cellstr(Params{5,j}));
            dia_out         = str2double(cellstr(Params{6,j}));
            d               = str2double(cellstr(Params{10,j}));
            M               = str2double(cellstr(Params{11,j}));
            
            
            % Friction calculations
            A_pipe      =   pi/4 * dia_inlet^2 ;
            [x, index] = unique(ComponentsMain(1,:));
            A_tunnel    =   interp1(x,ComponentsMain(2,index),intercept);    % Finding the cross-section of the main tunnel at the location of the intake
        
            %tau_st (Table 7-20) from Idelchik
            if (A_pipe/A_tunnel) <= 0.4
                tau_st = 0.4;            
            elseif (A_pipe/A_tunnel) > 0.4
                if (Q_pipe/Q_pump) <= 0.5        
                    tau_st = 2*(2*(Q_pipe/Q_pump) - 1); 
                elseif (Q_pipe/Q_pump) > 0.6                   
                    tau_st = 0.3*(2*(Q_pipe/Q_pump) - 1);
                end        
            end
            
            %A_prime (Table 7-4) from Idelchik
            if (A_pipe/A_tunnel) <= 0.35
                if (Q_pipe/Q_pump) <= 0.4
                    A_prime = 1.1 - 0.7*(Q_pipe/Q_pump);                
                elseif (Q_pipe/Q_pump) > 0.4                       
                    A_prime = 0.85;
                end            
            elseif (A_pipe/A_tunnel) > 0.35
                if (Q_pipe/Q_pump) <= 0.6        
                    A_prime = 1.0 - 0.6*(Q_pipe/Q_pump); 
                elseif (Q_pipe/Q_pump) > 0.6                   
                    A_prime = 0.6;
                end        
            end
        
            zeta_cs     = A_prime * (1 + ((Q_pipe / (Q_pump - Q_pipe))*(A_pipe / A_tunnel))^2 - (2 * ((Q_pipe / (Q_pump - Q_pipe))*(A_pipe / A_tunnel)) * cos(deg2rad(Alphaangle))));
            zeta_cst    = tau_st * (Q_pipe/(Q_pump - Q_pipe))^2;
            
            
            zeta_s  = zeta_cs / ((Q_pipe/Q_pump)*(A_tunnel/A_pipe))^2;
            zeta_st = zeta_cst / ((1 - (Q_pipe/Q_pump))^2  * (A_tunnel/A_pipe));
        
            zeta_divergence = ((Q_pipe/Q_pump) * zeta_s) + (((Q_pump - Q_pipe)/Q_pump) * zeta_st);
                    
            % Part 2) Friction due to horizontal pipe section
            Re_0   = (Q_pipe/(pi/4 * dia_inlet^2) * dia_inlet)/nu;
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/dia_inlet)/3.7)^1.11))^(-2);
            % multiply length by 3 to account for flexible piping
            zeta_hor = lambda * 3* length_inlet/dia_inlet;
            
            
            % Part 3) Friction due to main lunchbox (Contraction
            % and expansion losses
            expansion = (1 - ((pi/4*dia_inlet^2)/(width_lunchbox^2)) )^2;
            contraction = 0.5*(1 - ((pi/4*dia_out^2)/(width_lunchbox^2)) )^0.75;
            
            % Part 4) Friction due to screen placed at the inlet section

            porosity = (1 - (d/M))^2;
            Re_d = d * (Q_pipe/(pi/4 * dia_out^2))/nu;           
            if (porosity > 0.55) 
                disp(['Porosity of ' num2str(porosity) ' greater than 0.55. Increase d/M.'])
            end
            
            
            % Friction factor calculations are taken from the paper "Grid
            % -generated turbulence revisited" by Kurian and Fransson
        
            if (35 < Re_d || Re_d < 800)
                f = 0.5 + (26 * (Re_d^(-1))); 
            elseif Re_d > 800
                f = 0.52 + (66 * (Re_d^(-4/3)));
            elseif Re_d < 35
                f = Re_d^(-1);
            end
        
            screenloss = f * ((1 - porosity^2)/ porosity^2);
            
            % Modify friction factors to account for variation of
            % cross-sections within the unit
            contraction = contraction * (pi/4 * dia_inlet^2 / (width_lunchbox^2))^2;
            screenloss = screenloss * (pi/4 * dia_inlet^2 / (width_lunchbox^2))^2;
 
            % Calculation of friction factor
    
            fric_facs(j) = zeta_divergence + zeta_hor + expansion + contraction + screenloss;
        
        
        elseif strcmp(type_cross,'outlet') == 1
    
    
            disp('Merging connection for pipe')
        
            % NOTE: For the intake, the length input is, infact, the angular
            % intake from the main pipe or riser
            % Input data
            remerge         = str2double(cellstr(Params{12,j}));
            dia_out         = str2double(cellstr(Params{4,j}));                        
            alfa            = deg2rad(str2double(cellstr(Params{7,j})));     
            length_outlet   = str2double(cellstr(Params{3,j}));
            width_lunchbox  = str2double(cellstr(Params{5,j}));
            dia_in          = str2double(cellstr(Params{6,j}));
            
            
            % Friction calculations
            A_pipe      =   pi/4 * dia^2; 
            [x, index]  =   unique(ComponentsMain(1,:));
            A_tunnel    =   interp1(x,ComponentsMain(2,index),remerge);    % Finding the cross-section of the main tunnel at the location of the convergence
        
      
            % Part 1) Friction due to the convergence
            %A (Table 7-1) from Idelchik
            if (A_pipe/A_tunnel) <= 0.35
                if (Q_pipe/Q_pump) <= 1
                    A = 1;                
                end            
            elseif (A_pipe/A_tunnel) > 0.35
                if (Q_pipe/Q_pump) <= 0.4        
                    A = 0.9*(1.0 - (Q_pipe/Q_pump)); 
                elseif (Q_pipe/Q_pump) > 0.4                   
                    A = 0.55;
                end        
            end
        
            %Kstprime (Table 7-1) from Idelchik
            if (A_pipe/A_tunnel) <= 0.35
                Kstprime = 0.8 * (Q_pipe/Q_pump);                      
            elseif (A_pipe/A_tunnel) > 0.35
                if (Q_pipe/Q_pump) <= 0.6       
                    Kstprime = 0.5; 
                elseif (Q_pipe/Q_pump) > 0.6                  
                    Kstprime = 0.8 * (Q_pipe/Q_pump);
                end        
            end
        
            zeta_cs     = A * (1 + ((Q_pipe / (Q_pump))*(A_tunnel / A_pipe))^2 - (2 * (1 - (Q_pipe / (Q_pump))))    - (2 * ((Q_pipe / (Q_pump)) * (A_pipe / A_tunnel)^2 * cos(alfa))) );
            zeta_cst    = 1 - (1 - (Q_pipe / (Q_pump)) )^2    - ((1.4 - (Q_pipe / (Q_pump)) )   * (Q_pipe / (Q_pump))^2 * sin(alfa)) - (2*Kstprime* cos(alfa) * A_tunnel/Q_pump);
        
            zeta_s  = zeta_cs / ((Q_pipe/Q_pump)*(A_tunnel/A_pipe))^2;
            zeta_st = zeta_cst / (1 - (Q_pipe/Q_pump)^2); 
        
            zeta_convergence = ((Q_pipe/Q_pump) * zeta_s) + (((Q_pump - Q_pipe)/Q_pump) * zeta_st);
                  
            % Part 2) Friction due to flexible hose connection
            Re_0   = (Q_pipe/(pi/4 * dia_out^2) * dia_out)/nu;
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/dia_out)/3.7)^1.11))^(-2);
            % multiply length by 3 to account for flexible piping
            zeta_hor = lambda * 3*length_outlet/dia_out;
            
            % Part 3) Friction due to main lunchbox (Contraction, bending
            % and expansion losses
            expansion = (1 - ((pi/4*dia_in^2)/(width_lunchbox^2)) )^2;
            contraction = 0.5*(1 - ((pi/4*dia_out^2)/(width_lunchbox^2)) )^0.75;
            
            % Calculation of friction factor
            
            fric_facs(j) = zeta_convergence + zeta_hor + expansion + contraction;
                
        end
        
    end
    
    if strcmp(type_comp,'Dynamometer') == 1
        disp('Dynamometer channel')
        % NOTE: We represent the dynamometer with a sphere, for
        % conservative design. We add the channel length it is located in
        % here as well.
        % Input data
        sphere_dia  = str2double(cellstr(Params{10,j}));         % Sphere diameter in the tunnel
        if strcmp(type_cross,'pipe') == 1
            dia      = str2double(cellstr(Params{4,j}));         % Duct height                       
            w        = Q_comp/(pi/4 * dia^2);
            D_h      = dia;
        elseif strcmp(type_cross,'ducted') == 1
            height   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            width    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            D_h      = 0.64 * 4*height*width/(2*(height+width));       % 64% of the hydraulic diameter, or simply the effective diameter
            w        = Q_comp/(pi/4 * D_h^2);
        end
        
        % Friction by the channel
        Re  = w*D_h/nu;
        lambda      = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
        f_channel   = lambda*lengths(j)/D_h;  
        
        % Friction by the sphere
        Re_sph   = sphere_dia * w/nu;
        
        % From the correlation of Morris:
        Cd       = (24/Re_sph) + ((2.6*Re_sph/5)/(1 + (Re_sph/5)^1.52)) + ((0.411 * (Re_sph/ 2.63e5)^(-7.94) )/(1 + ((Re_sph/ 2.63e5))^-8)) + ((0.25*Re_sph/1e6)/(1 + (Re_sph/1e6)));
        
        % Adding individual contributions
        zeta     = Cd + f_channel;
        fric_facs(j) = zeta;
            
        
    end    
    
    counter = counter+1;
end



% if TunnelSolidity > 0
%     theta    = asin(TunnelSolidity * 0.3 / ValveDia); % Opening angle of the valve
%     % Display warning for lightly opened throttling valve
%     if rad2deg(theta) < 25
%         disp(['The valve opening angle ' num2str(theta) ' is less than 25 degrees. Friction factor will not be accurately calculated.'])
%     end
%     Re      = ((Q_pump - Q_pipe)/(0.3^2)) * 0.3 / nu;
%     Dd_bar  = ValveDia/0.3; 
%     A       = abs(60 * ((1 + (0.5 * Dd_bar * (1 + sin(theta)))) / (1 - (Dd_bar^2 * sin(theta))^2 )));
%     zeta_qu = (((1.56) / (1 - (Dd_bar*sin(theta)))) - 1)^2;
%     
%     zeta_throttling = abs((A/Re) + ((1- (50/Re)) * zeta_qu));
%     
%     fric_facs = [fric_facs zeta_throttling]; 
%         
% end


end

