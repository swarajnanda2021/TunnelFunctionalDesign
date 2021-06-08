function [ TI ] = TurbRed( Q_comp, Params, TI_Inlet, NN, nu )
%TurbRed: Calculates the turbulence reduction factor of various turbulence
%control devices after matchin with them in Params.


n_comp      = size(Params,2);
turbcounter = 1;
% counter = 0;

% Setting up the TI Matrix
TI = [];
TI = [TI TI_Inlet];

% Assign lengths array
lengths = zeros(size(Params,2));
for k = 1:size(Params,2) 
    lengths(k) = str2double(cellstr(Params{3,k}));
end


for j = 1:n_comp
     type_comp = cellstr(Params{1,j});
     type_cross= cellstr(Params{2,j});   
     if strcmp(type_comp,'Mesh') == 1
        
        if strcmp(type_cross,'ducted') == 1
%             disp('Mesh in a duct')
            % Input params
            height      = str2double(cellstr(Params{4,j}));
            width       = str2double(cellstr(Params{6,j}));
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
                        
            % Friction calculations
            porosity = (1 - (d/M))^2;
            solidity = 1 - porosity;
            Re_d = d * (Q_comp(j)/(height*width))/nu;
            Re_M = M * (Q_comp(j)/(height*width))/nu;            
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
        
            fric_facs = f * ((1 - porosity^2)/ porosity^2);
            
            
        elseif strcmp(type_cross,'pipe') == 1
%             disp('Mesh in a pipe')
            % Input params
            dia         = str2double(cellstr(Params{4,j}));
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
            
            % Friction calculations
            porosity = (1 - (d/M))^2;
            solidity = 1 - porosity;
            Re_d = d * (Q_comp(j)/(pi/4 * dia^2))/nu;
            Re_M = M * (Q_comp(j)/(pi/4 * dia^2))/nu;            
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
        
            fric_facs = f * ((1 - porosity^2)/ porosity^2); 
            
    
    
        end
        RedFac = (1+fric_facs)^-1 + (2*(1+fric_facs)^-2);
        
        TI = [TI (TI(end) * RedFac)];
        
        
        turbcounter = turbcounter+1;
     end

     if strcmp(type_comp,'ContractionSingle') == 1
         if strcmp(type_cross,'ducted') == 1
%             disp('Contraction in a duct')
            % Input data
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information 
            CR = (height_in*width_in)/(height_out*width_out); 
            Mu = (3/4) * (CR^-2) * (log((4 * CR^3) - 1));
            Nu =  (3/4) * CR;
            RedFac = (Mu + 2*Nu) / (3 * CR^2);
            
            
            
            
         end
         
         if strcmp(type_cross,'pipe') == 1
%             disp('Contraction in a pipe')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j})); 
            CR = (dia_in^2)/(dia_out^2); 
            Mu = (3/4) * (CR^-2) * (log((4 * CR^3) - 1));
            Nu =  (3/4) * CR;
            RedFac = (Mu + 2*Nu) / (3*CR^2);
            
            
     
         end
         
         if strcmp(type_cross,'C2R') == 1
%             disp('Contraction in a pipe')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            height_out  = str2double(cellstr(Params{5,j})); 
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information 
            CR = (pi/4 * dia_in^2)/(height_out*width_out); 
            Mu = (3/4) * (CR^-2) * (log((4 * CR^3) - 1));
            Nu =  (3/4) * CR;
            RedFac = (Mu + 2*Nu) / (3*CR^2);
            
         end
         
         % Assuming that a settling chamber is located before hand, we also
         % take into account the decay of turbulence along the settling
         % length, i.e., from the plenum (end of honeycomb section) to the
         % contraction. (Wetzels paper)
         
         if strcmp(cellstr(Params{1,j-2}),'Honeycomb')
            LPlenum = lengths(j-1);
            xPlenum  = lengths(j-2) - str2double(cellstr(Params{12,j-2})); % Decay length of the honeycomb before it
            xTestSec = (LPlenum + (1/3)*lengths(j));
            RedFacPlenum = xPlenum / xTestSec;
            TI = [TI (TI(end) * (RedFac * RedFacPlenum)^0.5)];
         else
            TI = [TI TI(end)* (RedFac^0.5)];
         end
         turbcounter = turbcounter+1;
     end
     
     
     if strcmp(type_comp,'ContractionDouble') == 1
         if strcmp(type_cross,'ducted') == 1
%             disp('Contraction in a duct')
            % Input data
            CR1x         = str2double(cellstr(Params{10,j}));
            CR1y         = str2double(cellstr(Params{11,j}));
            CR2          = str2double(cellstr(Params{12,j}));
            CR1 = (CR1x + CR1y)/2;
            Mu1 = (3/4) * (CR1^-2) * (log((4 * CR1^3) - 1));
            Nu1 =  (3/4) * CR1;
            RedFac1 = (Mu1 + 2*Nu1) / (3 * CR1^2);
            
            Mu2 = (3/4) * (CR2^-2) * (log((4 * CR2^3) - 1));
            Nu2 =  (3/4) * CR2;
            RedFac2 = (Mu2 + 2*Nu2) / (3 * CR2^2);
            
            RedFac = RedFac1*RedFac2;
            
            
         end
         
         if strcmp(type_cross,'pipe') == 1
%             disp('Contraction in a pipe')
            % Input data
            % Input data
            CR1         = str2double(cellstr(Params{10,j}));
            CR2         = str2double(cellstr(Params{11,j}));
            
            Mu1 = (3/4) * (CR1^-2) * (log((4 * CR1^3) - 1));
            Nu1 =  (3/4) * CR1;
            RedFac1 = (Mu1 + 2*Nu1) / (3 * CR1^2);
            
            Mu2 = (3/4) * (CR2^-2) * (log((4 * CR2^3) - 1));
            Nu2 =  (3/4) * CR2;
            RedFac2 = (Mu2 + 2*Nu2) / (3 * CR2^2);
            
            RedFac = RedFac1*RedFac2;
            
     
         end
         
         if strcmp(type_cross,'C2R') == 1
%             disp('Contraction in a pipe')
            % Input data
            CR1         = str2double(cellstr(Params{10,j}));
            CR2         = str2double(cellstr(Params{11,j}));
            
            Mu1 = (3/4) * (CR1^-2) * (log((4 * CR1^3) - 1));
            Nu1 =  (3/4) * CR1;
            RedFac1 = (Mu1 + 2*Nu1) / (3 * CR1^2);
            
            Mu2 = (3/4) * (CR2^-2) * (log((4 * CR2^3) - 1));
            Nu2 =  (3/4) * CR2;
            RedFac2 = (Mu2 + 2*Nu2) / (3 * CR2^2);
            
            RedFac = RedFac1*RedFac2;
            
         end
         
         % Assuming that a settling chamber is located before hand, we also
         % take into account the decay of turbulence along the settling
         % length, i.e., from the plenum (end of honeycomb section) to the
         % contraction. (Wetzels paper)
         
         if strcmp(cellstr(Params{1,j-2}),'Honeycomb')
            LPlenum = lengths(j-1);
            xPlenum  = lengths(j-2) - str2double(cellstr(Params{12,j-2})); % Decay length of the honeycomb before it
            xTestSec = (LPlenum + (1/3)*lengths(j));
            RedFacPlenum = xPlenum / xTestSec;
            TI = [TI (TI(end) * (RedFac * RedFacPlenum)^0.5)];
         else
            TI = [TI TI(end)* (RedFac^0.5)];
         end
         turbcounter = turbcounter+1;
     end
     
     
     if strcmp(type_comp,'Honeycomb') == 1
        
        % Remove the turbulence reduction caused by the mesh, if upstream of the honeycomb
        
%         if strcmp(Params{1,j-1},'Mesh')==1
%             TI = [TI_Inlet];
%         end   
        
        TI = [TI_Inlet];
         
        if strcmp(type_cross,'ducted') == 1
%             disp('Honeycomb in a duct')
            height = str2double(cellstr(Params{4,j}));
            width  = str2double(cellstr(Params{6,j}));
            D_h    = (2 * height*width/(height+width));
        
        
            l = str2double(cellstr(Params{12,j}));
            D = str2double(cellstr(Params{11,j}));
    
            L = str2double(cellstr(Params{13,j}));     % Length scale of inlet turbulence, set as constant due to the mesh upstream
    
            lbyD = l/D;
            Lbyl = L/l;
    
            Re_d = (Q_comp(j)/(pi*D^2/4)) * D / nu;
    
%             data = csvread('Lumley_honeycomb_pressuredrop.txt');
%             F = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
%             fric_facs = F(lbyD,Re_d);
%             
            
            % New approach: Friction factor computed from the Haaland's
            % explicit format for a single cell flow, reduction factor
            % directly calculated from the formula
            
            % Friction factor is taken to be the pressure drop in a single
            % honeycomb cell, assuming it to be a circular mesh.
            k_friction = 0.01e-3;
            f = (-1.8 * log10(6.9/Re_d + ((k_friction/D)/3.7)^1.11))^(-2)* l/D;

            
            % The turbulence reduction factor is calculated directly from
            % the formula given by Lumley's analysis.
            k = f;
            fun = @(y) 8*((k/(3*pi))*(Lbyl))^2.*((4*k*L*(k+1)/(3*pi*l))^2 + y.^2.*(1 + ( (8*k/(3*pi)) .*(Lbyl) .* (1 - y.^2).^(-0.5) ) ).^2).^(-1) ;
            
            red_fac = integral(fun,0,1);
            
            
            
%             data_tur_red = csvread('Lumley_honeycomb_turbulence_reduction_isoreduction.txt');
%             G = scatteredInterpolant(data_tur_red(:,1),data_tur_red(:,2),data_tur_red(:,3));
%             sqrt_eta = G(fric_facs,Lbyl);
%             red_fac = sqrt_eta^2; 
            
            
            
            
            
            % Decay of honeycomb turbulence
            lengths_tur  = lengths(j) - l;
    
    
            if lengths_tur > 0
                if lengths_tur < 20*D
                    disp(['Decay length too small since lbyD is ' num2str(lengths_tur/D) '.']);
                end
                
                Decay_length_array = linspace(0,lengths_tur,1000);
                
                Re_l = (Q_comp(j)/(pi*D^2/4)) * l/nu; % Calculating the reynolds number of the flow in each cell
                delta = (l*0.37) * (Re_l^(-1/3));
    
                for m=1:length(Decay_length_array)
                    if (Decay_length_array(m)/D) < ((Re_d/16)^(5/3)) % Initial period of decay
%                         disp('Initial period of honeycomb tur decay')
                        if delta >= D/2  % Turbulent cell flow
                            TI_decay(m) = sqrt(0.0072/(Decay_length_array(m)/D));
                        else % Laminar cell flow
                            TI_decay(m) = sqrt(0.03/(Decay_length_array(m)/D));
                        end
    
    
                    else % Final period of decay
%                         disp('Final period of honeycomb tur decay')
                        TI_decay(m) = sqrt((Decay_length_array(m)/D)^(-5/3));
                    end
                end
    
            else
                Decay_length_array = linspace(0,lengths_tur,1000);
            end   
            
            if exist('RedFac') == 0
                TI_inflow = TI_Inlet;
            else
                TI_inflow = TI_Inlet * prod(RedFac);
            end
            
            Final_TI = TI_decay(length(Decay_length_array)) +  red_fac*TI_inflow ;
            % Effective reduction factor
            TI = [TI Final_TI];
            
            turbcounter = turbcounter+1;
            
            
            
        end
            
        if strcmp(type_cross,'pipe') == 1
%           disp('Honeycomb in a pipe')
            dia = str2double(cellstr(Params{4,j}));
            
        
            l = str2double(cellstr(Params{12,j}));
            D = str2double(cellstr(Params{11,j}));
    
            L = str2double(cellstr(Params{13,j}));     % Length scale of inlet turbulence, set as constant due to the mesh upstream
    
            lbyD = l/D;
            Lbyl = L/l;
            
            
            Re_d = (Q_comp(j)*(D^2/dia^2)/(pi/4 * D^2)) * D / nu;   % From Reynolds number of the flow in the pipe, avoiding Lumley's graph
            
            
            % Older approach: Use the diagrams of Lumleu
%             data = csvread('Lumley_honeycomb_pressuredrop.txt');
%             F = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
%             fric_facs = F(lbyD,Re_d);
    
            % Turbulence calculations
%             data_tur_red = csvread('Lumley_honeycomb_turbulence_reduction_isoreduction.txt');
%             G = scatteredInterpolant(data_tur_red(:,1),data_tur_red(:,2),data_tur_red(:,3));
%             sqrt_eta = G(fric_facs,Lbyl);
%             red_fac = sqrt_eta^2; 
%     
            
            % New approach: Friction factor computed from the Haaland's
            % explicit format for a single cell flow, reduction factor
            % directly calculated from the formula
            
            % Friction factor is taken to be the pressure drop in a single
            % honeycomb cell, assuming it to be a circular mesh.
            k_friction = 0.01e-3;
            f = (-1.8 * log10(6.9/Re_d + ((k_friction/D)/3.7)^1.11))^(-2)* l/D;

            
            % The turbulence reduction factor is calculated directly from
            % the formula given by Lumley's analysis.
            k = f;
            fun = @(y) 8*((k/(3*pi))*(Lbyl))^2.*((4*k*L*(k+1)/(3*pi*l))^2 + y.^2.*(1 + ( (8*k/(3*pi)) .*(Lbyl) .* (1 - y.^2).^(-0.5) ) ).^2).^(-1) ;
            
            red_fac = integral(fun,0,1);
            

            % Decay of honeycomb turbulence
            lengths_tur  = lengths(j) - L;
    
    
            if lengths_tur > 0
                if lengths_tur < 20*D
                    disp(['Decay length too small since lbyD is ' num2str(lengths_tur/D) '.']);
                end
                
                Decay_length_array = linspace(0,lengths_tur,1000);
                
    
                Re_l = (Q_comp(j)/(pi*D^2/4)) * l/nu; % Calculating the reynolds number of the flow in each cell
                delta = (l*0.37) * (Re_l^(-1/3));
    
                for m=1:length(Decay_length_array)
                    if (Decay_length_array(m)/D) < ((Re_d/16)^(5/3)) % Initial period of decay
%                         disp('Initial period of honeycomb tur decay')
                        if delta >= D/2  % Turbulent cell flow
                            TI_decay(m) = sqrt(0.0072/(Decay_length_array(m)/D));
                        else % Laminar cell flow
                            TI_decay(m) = sqrt(0.03/(Decay_length_array(m)/D));
                        end
    
    
                    else % Final period of decay
%                         disp('Final period of honeycomb tur decay')
                        TI_decay(m) = sqrt((Decay_length_array(m)/D)^(-5/3));
                    end
                end
    
            else
                Decay_length_array = linspace(0,lengths_tur,1000);
            end    
            
            if exist('RedFac') == 0
                TI_inflow = TI_Inlet;
            else
                TI_inflow = TI_Inlet * prod(RedFac(1:turbcounter-1,1));
            end
            Final_TI = TI_decay(end) +  red_fac*TI_inflow ;
            % Effective reduction factor
            TI = [TI Final_TI];
             
        
            turbcounter = turbcounter+1;
            
        end 
            
     end
     
     
end




end

