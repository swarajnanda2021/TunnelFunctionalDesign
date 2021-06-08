function [components, intercept, remerge, volume, f0byfi, A_comp,ContractionWidth_stage2,lengths] = geom_tunnel( Params , NN , HTs)
% SUMMARY: Code for constructing the geometry matrix for both the main
% tunnel as well as the pipe section, taking information from the cell
% matrix Tunnel_params.

ContractionWidth_stage2 = 0;
A_comp = [];
counter=0;
intercept = 0;
remerge   = 0;
n_comp = size(Params,2);


% Assign lengths array
lengths = zeros(size(Params,2));
for k = 1:size(Params,2) 
    lengths(k,1) = str2double(cellstr(Params{3,k}));
end

% Pre-assigning the volume array
volume = zeros(size(lengths));
f0byfi = zeros(size(lengths));
indexcounter    = 0;
lengthcounter   = 0;

for j=1:n_comp
    
    % Input component categorical information
    type_comp = cellstr(Params{1,j});
    type_cross= cellstr(Params{2,j});
    
    % Lunchbox element with a simple hose connection
    if strcmp(type_comp,'Lunchbox') == 1
        
        if strcmp(type_cross,'inlet') == 1
            % Input data
            length_inlet    = str2double(cellstr(Params{3,j}));
            dia_inlet       = str2double(cellstr(Params{4,j}));                        
            width_lunchbox  = str2double(cellstr(Params{5,j}));
            dia_out         = str2double(cellstr(Params{6,j}));
            Alphaangle      = str2double(cellstr(Params{7,j}));
            datum_in        = str2double(cellstr(Params{8,j}));                                       
            datum_out       = str2double(cellstr(Params{9,j}));
            d               = str2double(cellstr(Params{10,j}));
            M               = str2double(cellstr(Params{11,j}));
            intercept       = str2double(cellstr(Params{12,j}));
                                      
            % Modify length
            lengths(j,1) = length_inlet + width_lunchbox;
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            for i=1:NN
                if components(1,1+(indexcounter*NN)+(i-1)) < length_inlet
                    components(2,1+(indexcounter*NN)+(i-1))   = pi/4 * (dia_inlet)^2;
                else
                    components(2,1+(indexcounter*NN)+(i-1))   = width_lunchbox^2;
                end
            end
                
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_inlet, width_lunchbox, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_inlet, width_lunchbox, NN);
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        elseif strcmp(type_cross,'outlet') == 1
            
            remerge         = str2double(cellstr(Params{12,j}));
            length_outlet   = str2double(cellstr(Params{3,j}));
            dia_outlet      = str2double(cellstr(Params{4,j}));                        
            width_lunchbox  = str2double(cellstr(Params{5,j}));
            dia_in          = str2double(cellstr(Params{6,j}));
            alfa            = deg2rad(str2double(cellstr(Params{7,j})));
            
            % Modify length
            lengths(j,1) = length_outlet + width_lunchbox;
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
           
            % Assigning cross-sectional area array of component matrix
            for i=1:NN
                if components(1,1+(indexcounter*NN)+(i-1)) - lengthcounter < width_lunchbox
                    components(2,1+(indexcounter*NN)+(i-1))   = width_lunchbox^2;
                else
                    components(2,1+(indexcounter*NN)+(i-1))   = pi/4 * (dia_outlet)^2;
                end
            end

            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_outlet, width_lunchbox, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_outlet, width_lunchbox, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        end
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
    end
    
    
    
    if strcmp(type_comp,'Separator') == 1 
        if strcmp(type_cross,'Top') == 1
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
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = Area_channel;
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = D_h;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = D_h;
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
            
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
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = Area_channel;
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = D_h;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = D_h;
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);            
            
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
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4 * (D_h1)^2, pi/4 * (D_h2)^2, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(D_h1, D_h2, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(D_h1, D_h2, NN);
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
%             volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            volume(j) = 0.5*((4/3) * pi * (0.5*dia_separator)^3);
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN));
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
            
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
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4 * (D_h1)^2, pi/4 * (D_h2)^2, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(D_h1, D_h2, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(D_h1, D_h2, NN);
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
%             volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            volume(j) = 0.5*((4/3) * pi * (0.5*dia_separator)^3);
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN));
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);  
            
            
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
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4 * (dia_in)^2,pi/4 * (D_h)^2,NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = dia_in;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = dia_in;
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
            
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
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
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
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(pi/4 * (D_h)^2,pi/4 * (dia_out)^2,NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = dia_out;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = dia_out;
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
        elseif strcmp(type_cross,'ducted') == 1 
            
            nPlatesFunctional   = str2double(cellstr(Params{12,j}));
            nPlates             = str2double(cellstr(Params{11,j}));
            
            
            height_in   = str2double(cellstr(Params{4,j})) * nPlatesFunctional/nPlates;                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j})) * nPlatesFunctional/nPlates;                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(height_in*width_in, height_out*width_out, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(height_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(width_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);
            
        end
    end   
    
    
    if strcmp(type_comp,'Straight') == 1 
        
        if strcmp(type_cross,'pipe') == 1
            disp('Straight pipe')
            
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4 * (dia_in)^2, pi/4 * (dia_out)^2, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_in, dia_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_in, dia_out, NN);
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
        
        elseif strcmp(type_cross,'flexible') == 1
            disp('Flexible hose')
    
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4 * (dia_in)^2, pi/4 * (dia_out)^2, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_in, dia_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(dia_in, dia_out, NN);
            
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);           
            
                        
        elseif strcmp(type_cross,'ducted') == 1
            
            disp('Straight duct')
    
            
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(height_in*width_in, height_out*width_out, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(height_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(width_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);
        
            
        elseif strcmp(type_cross,'ductedProtuberance') == 1
            
            disp('Straight duct with protuberance loss')
    
            l_vent      = str2double(cellstr(Params{11,j}));     % length of the knife
            alpha_opening = str2double(cellstr(Params{10,j}));  % Opening angle (degrees)
            
            protuberance = l_vent*sin(deg2rad(alpha_opening)); % The opening angle
            
            height_in   = str2double(cellstr(Params{4,j}))  - protuberance;                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}))  - protuberance;                                          % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));     % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(height_in*width_in, height_out*width_out, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(height_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(width_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter +lengths(j,1);
        
            
            
            
            
        elseif strcmp(type_cross,'dissolutionPipe')
            
            disp('Dissolution or Separation section-Pipe')
    
            % Input data
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
            
            datum_in = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out= str2double(cellstr(Params{9,j}));                                       % Datum related information
               
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + 2*lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            
            components(2,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = Atop;
            components(2,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = Abottom;
            
            
            
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(3,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            
            components(4,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(4,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            % Assigning vertical of the centreline
            components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = datum_in;%components(3,:); 
            components(5,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = datum_out;
            
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter + 2*lengths(j,1);
        
        elseif strcmp(type_cross,'dissolutionBox')
            
            disp('Dissolution or Separation section-Box')
    
            % Input data
            height_sec  = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            width_sec   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            dia_in     = str2double(cellstr(Params{4,j}));
            height_out = str2double(cellstr(Params{5,j}));
            width_out  = str2double(cellstr(Params{10,j}));
            
            Abottom = dia_in*width_sec;
            Atop = (height_sec*width_sec) - Abottom;
            Pbottom = 2*(width_sec + dia_in);
            Ptop   = 2*(width_sec + height_out);
            
            D_h_in = 4*Atop/Ptop;
            D_h_out = 4*Abottom/Pbottom;
            
            
            datum_in = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out= str2double(cellstr(Params{9,j}));                                       % Datum related information
            
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + 2*lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            
            components(2,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = Abottom;
            components(2,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = Atop;
            
            
            
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(3,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            
            components(4,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(4,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            % Assigning vertical of the centreline
            components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = datum_in;%components(3,:); 
            components(5,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = datum_out;
            
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) + (pi/4*(dia_in+height_out)^2*0.5*width_out);    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter + 2*lengths(j,1);
        
        elseif strcmp(type_cross,'Separator')
            
            disp('Separator chamber')
    
            % Input data
            lseparator      = str2double(cellstr(Params{3,j}));                                           % Cross-section based information
            dia_in          = str2double(cellstr(Params{4,j}));
            dia_out         = str2double(cellstr(Params{5,j}));
            plate_allowance = str2double(cellstr(Params{6,j}));
            dia_sec         = str2double(cellstr(Params{7,j}));
            datum_in        = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out       = str2double(cellstr(Params{9,j}));                                       % Datum related information
               
            % Calculate separator properties from the above params
            
            height_top      = max(dia_in,dia_out);
            height_bottom   = dia_sec - (height_top + plate_allowance);
            
            %Angle subtended by the segment from the smaller channel height
            %side:
            theta       = acos((0.5*dia_sec - height_bottom)/(0.5*dia_sec));
            A_bottom    = ((dia_sec/2)^2/2) * (theta*pi/180 - sin(theta));
            A_top       = (pi/4 * dia_sec^2) - A_bottom - (plate_allowance*2*(0.5*dia_sec)*cos(theta));

            Pbottom = theta*dia_sec*0.5;
            Ptop   = pi*dia_sec - Pbottom;
            
            D_h_top = 4*Atop/Ptop;
            D_h_bottom = 4*Abottom/Pbottom;
            
            % Modify length of the separator, since you only specify
            % channel length in row 3.
            lengths(j,1) = lseparator*2 + (0.5*dia_sec*pi/2); % Considering only the 2 channels and r*theta from the bend
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            
            % Assigning cross-sectional area array of component matrix
            
            for iiii = 1+(indexcounter*NN):NN*(indexcounter+1)
               if (components(1,(indexcounter*NN)) - components(1,iiii)) <= lseparator % For first channel 
                  if dia_in > dia_out % Macrobubble separator
                      components(2,iiii) = A_top;
                      components(5,iiii) = datum_in;
                  else % Microbubble separator
                      components(2,iiii) = A_bottom;
                      components(5,iiii) = datum_in;
                  end
                   
               end
                
            end
            
            components(2,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = Atop;
            components(2,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = Abottom;
            
            % Assigning vertical of the centreline
            components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = datum_in;%components(3,:); 
            components(5,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = datum_out;
            
            
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(3,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            
            components(4,1+(indexcounter*NN):((indexcounter*NN) + NN/2))        = D_h_in;
            components(4,((indexcounter*NN) + NN/2 + 1) :NN*(indexcounter+1))   = D_h_out;
            
            
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            
            indexcounter = indexcounter+1;
            lengthcounter = lengthcounter + 2*lengths(j,1);    
            
        end
        
        
    end
    
    if strcmp(type_comp,'ContractionSingle') == 1
        
        
        if strcmp(type_cross,'ducted') == 1
            disp('Ducted contraction')
            % Input data
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            % Compute height profile using bell-mehta 5th order polynomials
            if height_in > height_out
                A = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
                b = [height_in/2,height_out/2,0,0,0,0];
                x = A\b';
                Contraction_eps = linspace(0,1,NN);
                Contraction_y = (x(1) .* Contraction_eps.^5) + (x(2) .* Contraction_eps.^4) + (x(3) .* Contraction_eps.^3) + (x(4) .* Contraction_eps.^2) + (x(5) .* Contraction_eps.^1) + (x(6));
            elseif height_in == height_out
                Contraction_y = ones(size(linspace(0,1,NN)))*height_in/2;
            else
                Contraction_y = ones(size(linspace(0,1,NN)));
            end
            % Compute width profile using bell-mehta 5th order polynomials
            A1 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b1 = [width_in/2,width_out/2,0,0,0,0];
            x1 = A1\b1';
            Contraction_eps1 = linspace(0,1,NN);
            Contraction_x = (x1(1) .* Contraction_eps1.^5) + (x1(2) .* Contraction_eps1.^4) + (x1(3) .* Contraction_eps1.^3) + (x1(4) .* Contraction_eps1.^2) + (x1(5) .* Contraction_eps1.^1) + (x1(6));
%             figure(4)
%             plot(Contraction_x)
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_y.*2).*(Contraction_x.*2));
%             figure(4)
%             plot((Contraction_y.*2))
%             hold on
%             plot((Contraction_x.*2))
            % Assigning width and height information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_y.*2));
            % Assigning width and height information
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_x.*2));
            % Assigning vertical position of the centreline                               
                      
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
                        
        elseif strcmp(type_cross,'pipe') == 1    
            
            disp('Conical contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j}));
    
            % The pipe section is at the highest position, it will have only
            % its cross sectional height at the contraction as the max-pressure
            datum_in    = str2double(cellstr(Params{8,j}));         
            datum_out   = str2double(cellstr(Params{9,j}));     
    
            % Compute radius profile using bell-mehta 5th order polynomials
            A = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b = [dia_in/2,dia_out/2,0,0,0,0];
            x = A\b';
            Contraction_eps = linspace(0,1,NN);
            Contraction_r = (x(1) .* Contraction_eps.^5) + (x(2) .* Contraction_eps.^4) + (x(3) .* Contraction_eps.^3) + (x(4) .* Contraction_eps.^2) + (x(5) .* Contraction_eps.^1) + (x(6));
    
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2).^2.*(pi/4));
            % Assigning radius information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            % Assigning datum information
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end      
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        
        elseif strcmp(type_cross,'C2R') == 1
            
            disp('Circle to Square Contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            % The pipe section is at the highest position, it will have only
            % its cross sectional height at the contraction as the max-pressure
            datum_in    = str2double(cellstr(Params{8,j}));         
            datum_out   = str2double(cellstr(Params{9,j}));     
    
            % Compute radius profile using bell-mehta 5th order polynomials
            A = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b = [dia_in/2,sqrt(height_out*width_out*4/pi)/2,0,0,0,0];
            x = A\b';
            Contraction_eps = linspace(0,1,NN);
            Contraction_r = (x(1) .* Contraction_eps.^5) + (x(2) .* Contraction_eps.^4) + (x(3) .* Contraction_eps.^3) + (x(4) .* Contraction_eps.^2) + (x(5) .* Contraction_eps.^1) + (x(6));
    
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr(pi/4 * (Contraction_r.*2).^2); % Not simple to calculate, take square section as conservative estimate
            
            % Assigning radius information MOSTLY IRRELEVANT INFORMATION
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            % Assigning datum information
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end      
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        end
            
    
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
    end
    
    if strcmp(type_comp,'ContractionDouble') == 1
        
        
        if strcmp(type_cross,'ducted') == 1
            disp('Ducted contraction')
            % Input data
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            % Contraction ratio parameter
            CR1x         = str2double(cellstr(Params{10,j}));
            CR1y         = str2double(cellstr(Params{11,j}));
            CR2          = str2double(cellstr(Params{12,j}));
            % Length of individual stages
            lengths1         = str2double(cellstr(Params{13,j}));
            lengths2         = str2double(cellstr(Params{14,j}));
            % Modify length
            lengths(j) = lengths1+lengths2;
            
            % Compute height profile using bell-mehta 5th order polynomials
            % First Stage
            if CR1y > 1
                A1 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
                b1 = [height_in/2,(sqrt((height_in)^2/CR1x))/2,0,0,0,0];
                x1 = A1\b1';
                Contraction_eps1 = linspace(0,1,ceil(NN/2));
                Contraction_y1 = (x1(1) .* Contraction_eps1.^5) + (x1(2) .* Contraction_eps1.^4) + (x1(3) .* Contraction_eps1.^3) + (x1(4) .* Contraction_eps1.^2) + (x1(5) .* Contraction_eps1.^1) + (x1(6));
            else
                Contraction_y = ones(size(linspace(0,1,ceil(NN/2))));
            end
            % Second Stage
            A2 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b2 = [(sqrt((height_in)^2/CR1x))/2,height_out/2,0,0,0,0];
            x2 = A2\b2';
            Contraction_eps2 = linspace(0,1,NN-ceil(NN/2));
            Contraction_y2 = (x2(1) .* Contraction_eps2.^5) + (x2(2) .* Contraction_eps2.^4) + (x2(3) .* Contraction_eps2.^3) + (x2(4) .* Contraction_eps2.^2) + (x2(5) .* Contraction_eps2.^1) + (x2(6));
            
            % Compute width profile using bell-mehta 5th order polynomials
            % First Stage
            A3 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b3 = [width_in/2,(sqrt((width_in)^2/CR1y))/2,0,0,0,0];
            x3 = A3\b3';
            Contraction_eps3 = linspace(0,1,ceil(NN/2));
            Contraction_y3 = (x3(1) .* Contraction_eps3.^5) + (x3(2) .* Contraction_eps3.^4) + (x3(3) .* Contraction_eps3.^3) + (x3(4) .* Contraction_eps3.^2) + (x3(5) .* Contraction_eps3.^1) + (x3(6));
            % Second Stage
            A4 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b4 = [(sqrt((width_in)^2/CR1y))/2,width_out/2,0,0,0,0];
            x4 = A4\b4';
            Contraction_eps4 = linspace(0,1,NN-ceil(NN/2));
            Contraction_y4 = (x4(1) .* Contraction_eps4.^5) + (x4(2) .* Contraction_eps4.^4) + (x4(3) .* Contraction_eps4.^3) + (x4(4) .* Contraction_eps4.^2) + (x4(5) .* Contraction_eps4.^1) + (x4(6));
            ContractionWidth_stage2 = Contraction_y4;
            % Assigning length array of component matrix
%             components(1,1+(indexcounter*NN):1+(indexcounter*NN)+(ceil(NN/2)))     = linspace(lengthcounter,lengthcounter + lengths1,ceil(NN/2));
%             components(1,1+(indexcounter*NN)+(ceil(NN/2))+1:NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths2,NN-ceil(NN/2));
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths1 + lengths2,NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr(([Contraction_y2 Contraction_y1].*2)).*fliplr(([Contraction_y4 Contraction_y3].*2));
            % Assigning width and height information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr(([Contraction_y2 Contraction_y1].*2));
            % Assigning width and height information
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr(([Contraction_y4 Contraction_y3].*2));
            % Assigning vertical position of the centreline                               
                      
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
                        
        elseif strcmp(type_cross,'pipe') == 1    
            
            disp('Conical contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j}));
    
            % The pipe section is at the highest position, it will have only
            % its cross sectional height at the contraction as the max-pressure
            datum_in    = str2double(cellstr(Params{8,j}));         
            datum_out   = str2double(cellstr(Params{9,j}));     
    
            CR1         = str2double(cellstr(Params{10,j}));
            CR2         = str2double(cellstr(Params{11,j}));
            
            % Length of individual stages
            length1         = str2double(cellstr(Params{12,j}));
            length2         = str2double(cellstr(Params{13,j}));
            % Modify length
            lengths(j) = lengths1+lengths2;
            
            % Compute radius profile using bell-mehta 5th order polynomials
            % First Stage
            A1 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b1 = [dia_in/2,sqrt((dia_in^2)/CR1)/2,0,0,0,0];
            x1 = A1\b1';
            Contraction_eps1 = linspace(0,1,ceil(NN/2));
            Contraction_r1 = (x1(1) .* Contraction_eps1.^5) + (x1(2) .* Contraction_eps1.^4) + (x1(3) .* Contraction_eps1.^3) + (x1(4) .* Contraction_eps1.^2) + (x1(5) .* Contraction_eps1.^1) + (x1(6));
            % Second Stage
            A2 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b2 = [sqrt((dia_in^2)/CR1)/2,dia_out/2,0,0,0,0];
            x2 = A2\b2';
            Contraction_eps2 = linspace(0,1,NN-ceil(NN/2));
            Contraction_r2 = (x2(1) .* Contraction_eps2.^5) + (x2(2) .* Contraction_eps2.^4) + (x2(3) .* Contraction_eps2.^3) + (x2(4) .* Contraction_eps2.^2) + (x2(5) .* Contraction_eps2.^1) + (x2(6));
    
            
            % Assigning length array of component matrix
%             components(1,1+(indexcounter*NN):1+(indexcounter*NN)+(ceil(NN/2)))     = linspace(lengthcounter,lengthcounter + lengths1,ceil(NN/2));
%             components(1,1+(indexcounter*NN)+(ceil(NN/2))+1:NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths2,NN-ceil(NN/2));
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths1 + lengths2,NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr([(Contraction_r1.*2) (Contraction_r2.*2)].^2.*(pi/4));
            % Assigning radius information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr([(Contraction_r1.*2) (Contraction_r2.*2)]);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr([(Contraction_r1.*2) (Contraction_r2.*2)]);
            % Assigning datum information
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end      
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        
        elseif strcmp(type_cross,'C2R') == 1    
           
            disp('Circle to Square Contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            % The pipe section is at the highest position, it will have only
            % its cross sectional height at the contraction as the max-pressure
            datum_in    = str2double(cellstr(Params{8,j}));         
            datum_out   = str2double(cellstr(Params{9,j}));     
    
            CR1         = str2double(cellstr(Params{10,j}));
            CR2         = str2double(cellstr(Params{11,j}));
            
            % Length of individual stages
            length1         = str2double(cellstr(Params{12,j}));
            length2         = str2double(cellstr(Params{13,j}));
            % Modify length
            lengths(j) = lengths1+lengths2;
            
            % Compute radius profile using bell-mehta 5th order polynomials
            % First stage circle to square
            A1 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b1 = [dia_in/2,(sqrt(height_out^2*CR2)/2),0,0,0,0];
            x1 = A1\b1';
            Contraction_eps1 = linspace(0,1,ceil(NN/2));
            Contraction_r1 = (x1(1) .* Contraction_eps1.^5) + (x1(2) .* Contraction_eps1.^4) + (x1(3) .* Contraction_eps1.^3) + (x1(4) .* Contraction_eps1.^2) + (x1(5) .* Contraction_eps1.^1) + (x1(6));
            % Second stage square cross-sectioned
            A2 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b2 = [(sqrt(height_out^2*CR2)/2),height_out/2,0,0,0,0];
            x2 = A2\b2';
            Contraction_eps2 = linspace(0,1,NN-ceil(NN/2));
            Contraction_r2 = (x2(1) .* Contraction_eps2.^5) + (x2(2) .* Contraction_eps2.^4) + (x2(3) .* Contraction_eps2.^3) + (x2(4) .* Contraction_eps2.^2) + (x2(5) .* Contraction_eps2.^1) + (x2(6));
            
            Contraction_r = [Contraction_r1 Contraction_r2];
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):1+(indexcounter*NN)+(ceil(NN/2)))     = linspace(lengthcounter,lengthcounter + lengths1,ceil(NN/2));
            components(1,1+(indexcounter*NN)+(ceil(NN/2))+1:NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths2,NN-ceil(NN/2));
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2).^2); % Not simple to calculate, take square section as conservative estimate
            % Assigning radius information MOSTLY IRRELEVANT INFORMATION
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            % Assigning datum information
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end      
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ; 
        end
            
    
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
    end
    
    % Bends are assumed straight, but vertical values will still be used
    % for bernoulli calculations. 
    
    % Regarding bubble behaviour, we test
    % terminal velocity related things only in the top section and the
    % rest, we assume that there is no slip between the bubble and the
    % bulk, hence only dissolution is important
    
    
    if strcmp(type_comp,'Bend') == 1
        % Cross-section based information
        
        if strcmp(type_cross,'pipe') == 1
            disp('Pipe bend')
            % Input data
            dia         = str2double(cellstr(Params{4,j}));               % Cross-section based information
            datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
            
            % Modify lengths(j,1)
            lengths(j,1) = pi*dia/4;
            
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4 * (dia)^2;
            % Assigning height and width related information
            components(3:4,1+(indexcounter*NN):NN*(indexcounter+1)) = dia;
            % Assigning vertical position of the centreline
            % Assigning vertical position of the centreline
            
            if datum_in < datum_out
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
                components(5,(1+((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = datum_out;
            elseif datum_in > datum_out
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = datum_in;
                components(5,(1+((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
            else
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = datum_in;
                components(5,(1 + ((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        elseif strcmp(type_cross,'ducted') == 1
            disp('Ducted bend')
            % Input data
            
            
            height          = str2double(cellstr(Params{4,j}));               % Cross-section based information
            width           = str2double(cellstr(Params{6,j}));    
            datum_in        = str2double(cellstr(Params{8,j}));               % Datum related information
            datum_out       = str2double(cellstr(Params{9,j}));               % Datum related information
            
            % Modify lengths(j,1)
            lengths(j,1) = height;
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height*width;
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = height;            
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = width;
            
            
            
            % Assigning vertical position of the centreline
            if datum_in < datum_out
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
                components(5,(1+((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = datum_out;
            elseif datum_in > datum_out
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = datum_in;
                components(5,(1+((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
            else
                components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))       = datum_in;
                components(5,(1 + ((indexcounter*NN) + NN/2)):(NN*(indexcounter+1)))   = linspace(datum_in,datum_out,length(components(5,1+(indexcounter*NN):((indexcounter*NN) + NN/2))));
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        elseif strcmp(type_cross,'arbitrary') == 1
            disp('Arbitrary bend')
            dia         = str2double(cellstr(Params{4,j}));
            alpha0       = str2double(cellstr(Params{5,j}));              % The bending angle in degrees
            datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
            RBend       = str2double(cellstr(Params{6,j}));
            R0= RBend;
            
            % Modify lengths(j,1)
%             lengths(j,1) = (alpha0*R0)

            

            
    
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4 * (dia)^2;
            % Assigning diameter related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            % Assigning datum
            components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
    
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
           
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
    end
    
    if strcmp(type_comp,'Adapter') == 1
        
        % Transition from a circle to a rectangle
        if strcmp(type_cross,'C2R') == 1
            disp('Circle to Square adapter')
            
            
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                           % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                           % Datum related information
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(pi/4*dia_in^2,height_out*width_out,NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(dia_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(dia_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;           
        end
        
        % Transition from a circle to a rectangle
        if strcmp(type_cross,'R2C') == 1
            disp('Circle to Square adapter')
        
            % Input data
            dia_out      = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            height_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            width_in   = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                           % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                           % Datum related information
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(height_in*width_in,pi/4 * dia_out^2,NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(dia_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(dia_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
                       
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
        
    end
    
    
    if strcmp(type_comp,'Mesh') == 1
        
        if strcmp(type_cross,'ducted') == 1
            disp('Mesh in a duct')
            homogenitycheck = 1;
            % Input params
            height      = str2double(cellstr(Params{4,j}));
            width       = str2double(cellstr(Params{6,j}));
            datum_in    = str2double(cellstr(Params{8,j}));           % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));           % Datum related information
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
            homogenitycheck = M/height;
            if homogenitycheck > 0.1
                disp(['Ratio of mesh width to tunnel height is ' num2str(M/height) '. Reduce M if greater than 0.1.'])                
            end   
                        
            
            
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height*width;
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = height;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = width;
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        elseif strcmp(type_cross,'pipe') == 1
            disp('Mesh in a pipe')
             
            homogenitycheck = 1;
            % Input params
            dia         = str2double(cellstr(Params{4,j}));
            d           = str2double(cellstr(Params{10,j}));
            M           = str2double(cellstr(Params{11,j}));
            datum_in    = str2double(cellstr(Params{8,j}));           % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));           % Datum related information
            homogenitycheck = M/dia;
            if homogenitycheck > 0.1
                disp(['Ratio of mesh width to tunnel height is ' num2str(M/height) '. Reduce M if greater than 0.1.'])                
            end   
    
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4 * (dia^2);
            % Assigning diameter related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            % Assigning vertical position of the centreline
            components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
    end
    
    
    if strcmp(type_comp,'Honeycomb') == 1
        
        if strcmp(type_cross,'ducted') == 1
        
        
            height = str2double(cellstr(Params{4,j}));
            width  = str2double(cellstr(Params{6,j}));
                
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
                      
        
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height*width;
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = height;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = width;
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        
        elseif strcmp(type_cross,'pipe')
            % For a honeycomb, we use the data of Lumley's iso-reduction curves
            dia         = str2double(cellstr(Params{4,j}));    
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));      % Datum related information
           
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi* dia^2 / 4;
            % Assigning diameter related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = dia;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = dia;
            % Assigning vertical of the centreline
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;%components(3,:); 
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
    
    end


    
    
    if strcmp(type_comp,'Intake') == 1
        disp('Intake pipe for pipe section')
        % NOTE: For the intake, the length input is, infact, the angular
        % intake from the main pipe or riser
        % Input data
        intercept   = str2double(cellstr(Params{2,j}));
        dia         = str2double(cellstr(Params{4,j}));                        
        alpha0       = str2double(cellstr(Params{10,j}));               % The divergence angle in degrees
        datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
        datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
        % Since we do not input the length, we calculate length as:
        length_intake = abs((datum_out - datum_in)/cos(deg2rad(alpha0)));
        
        
        % Assigning length array of component matrix
        components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + length_intake,NN);
        % Assigning cross-sectional area array of component matrix
        components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4 * (dia)^2;
        % Assigning diameter related information
        components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
        components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
        % Assigning vertical position of the centreline
        components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
                   
        % Volume Calculation
        volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
                
        % f0byfi Calculation
        f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + length_intake;
        
    end
    
    
    if strcmp(type_comp,'Remerge') == 1
        disp('Merging connection for pipe')
        
        % NOTE: For the intake, the length input is, infact, the angular
        % intake from the main pipe or riser
        % Input data
        remerge   = str2double(cellstr(Params{2,j}));
        dia         = str2double(cellstr(Params{4,j}));                        
        alpha0       = str2double(cellstr(Params{10,j}));               % The divergence angle in degrees
        datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
        datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
        % Since we do not input the length, we calculate length as:
        length_intake = abs((datum_out - datum_in)/cos(deg2rad(alpha0)));
        
        
        % Assigning length array of component matrix
        components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + length_intake,NN);
        % Assigning cross-sectional area array of component matrix
        components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4 * (dia)^2;
        % Assigning diameter related information
        components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
        components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
        % Assigning vertical position of the centreline
        components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
                   
        % Volume Calculation
        volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
        
        % f0byfi Calculation
        f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + length_intake;
        
    end
    
    if strcmp(type_comp,'Dynamometer') == 1
        disp('Dynamometer channel')
        % NOTE: We represent the dynamometer with a sphere, for
        % conservative design. We add the channel length it is located in
        % here as well.
        % Input data
        if strcmp(type_cross,'pipe') == 1
            dia         = str2double(cellstr(Params{4,j}));         % Duct height                       
            datum_in    = str2double(cellstr(Params{8,j}));         % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));         % Datum related information


            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = pi/4*(dia^2);
            % Assigning diameter related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = dia;
            % Assigning vertical position of the centreline
            components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);

            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe

            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;

        elseif strcmp(type_cross,'ducted') == 1

            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j,1),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(height_in*width_in, height_out*width_out, NN);
            % Assigning height and width related information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(height_in, height_out, NN);
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(width_in, width_out, NN);
            % Assigning vertical position of the centreline
            
            if datum_in == datum_out
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = datum_in;
            else
                components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
            end
            
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (HTs^2)/components(2,1+(indexcounter*NN)) ;
            
            
            
        end
        
        
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j,1);
        
    end
    
    A_comp(j)    = components(2,((counter)*NN) + 1);
    counter      = counter+1;
    
    
end





end

