function [components, intercept, remerge, volume, f0byfi, A_comp] = geom_tunnel( Params , NN)
% SUMMARY: Code for constructing the geometry matrix for both the main
% tunnel as well as the pipe section, taking information from the cell
% matrix Tunnel_params.

counter=0;
intercept = 0;
remerge   = 0;
n_comp = size(Params,2);

% Assign lengths array
lengths = zeros(size(Params,2));
for k = 1:size(Params,2) 
    lengths(k) = str2double(cellstr(Params{3,k}));
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
    
    
    if strcmp(type_comp,'Straight') == 1 
        
        if strcmp(type_cross,'pipe') == 1
            disp('Straight pipe')
    
            % Input data
            dia_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            dia_out = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                       % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                       % Datum related information
            
                
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
            
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
            
            
                        
        elseif strcmp(type_cross,'ducted') == 1
            
            disp('Straight duct')
    
            
            height_in   = str2double(cellstr(Params{4,j}));                                          % Cross-section based information
            height_out  = str2double(cellstr(Params{5,j}));                                         % Cross-section based information
            width_in    = str2double(cellstr(Params{6,j}));                                            % Cross-section based information
            width_out   = str2double(cellstr(Params{7,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
            
        end
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
    end
    
    if strcmp(type_comp,'Contraction') == 1
        
        
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
            A = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b = [height_in/2,height_out/2,0,0,0,0];
            x = A\b';
            Contraction_eps = linspace(0,1,NN);
            Contraction_y = (x(1) .* Contraction_eps.^5) + (x(2) .* Contraction_eps.^4) + (x(3) .* Contraction_eps.^3) + (x(4) .* Contraction_eps.^2) + (x(5) .* Contraction_eps.^1) + (x(6));
            
            % Compute width profile using bell-mehta 5th order polynomials
            A1 = [[1,1,1,1,1,1];[0,0,0,0,0,1];[0,0,0,0,1,0];[0,0,0,1,0,0];[5,4,3,2,1,0];[20,12,6,2,0,0]];
            b1 = [width_in/2,width_out/2,0,0,0,0];
            x1 = A1\b1';
            Contraction_eps1 = linspace(0,1,NN);
            Contraction_x = (x1(1) .* Contraction_eps1.^5) + (x1(2) .* Contraction_eps1.^4) + (x1(3) .* Contraction_eps1.^3) + (x1(4) .* Contraction_eps1.^2) + (x1(5) .* Contraction_eps1.^1) + (x1(6));
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_y.*2).*(Contraction_x.*2));
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
                        
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
            components(1,1+(indexcounter*NN):NN*(indexcounter+1)) = linspace(lengthcounter,lengthcounter + lengths(j),NN);
            % Assigning cross-sectional area information
            components(2,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2).^2.*(pi/4));
            % Assigning radius information
            components(3,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            components(4,1+(indexcounter*NN):NN*(indexcounter+1)) = fliplr((Contraction_r.*2));
            % Assigning datum information
                      
            % Volume Calculation
            volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
            
            % f0byfi Calculation
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
        end
            
    
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
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
            
            % Modify lengths(j)
            lengths(j) = pi*dia/4;
            
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
        elseif strcmp(type_cross,'ducted') == 1
            disp('Ducted bend')
            % Input data
            
            
            height          = str2double(cellstr(Params{4,j}));               % Cross-section based information
            width           = str2double(cellstr(Params{6,j}));    
            datum_in        = str2double(cellstr(Params{8,j}));               % Datum related information
            datum_out       = str2double(cellstr(Params{9,j}));               % Datum related information
            
            % Modify lengths(j)
            lengths(j) = height;
            
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
        elseif strcmp(type_cross,'arbitrary') == 1
            disp('Arbitrary bend')
            dia         = str2double(cellstr(Params{4,j}));
            alpha       = str2double(cellstr(Params{5,j}));               % The bending angle in degrees
            datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
            
            
            
            a = abs(datum_out - datum_in) * cos(deg2rad(90 - ((180 -alpha)/2))); 
            R0 = (0.5*a)/(sin(deg2rad(alpha)/2));
            % Modify lengths(j)
            lengths(j) = (deg2rad(alpha)*R0);
            
    
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + (deg2rad(alpha)*R0),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
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
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height_out*width_out;
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;           
        end
        
        % Transition from a circle to a rectangle
        if strcmp(type_cross,'R2C') == 1
            disp('Circle to Square adapter')
        
            % Input data
            dia_out      = str2double(cellstr(Params{5,j}));                                           % Cross-section based information
            height_in  = str2double(cellstr(Params{4,j}));                                           % Cross-section based information
            width_in   = str2double(cellstr(Params{6,j}));                                           % Cross-section based information
            
            datum_in    = str2double(cellstr(Params{8,j}));                                           % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));                                           % Datum related information
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
            % Assigning cross-sectional area array of component matrix
            components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height_out*width_out;
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
                       
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
        
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
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
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
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
    end
    
    
    if strcmp(type_comp,'Honeycomb') == 1
        
        if strcmp(type_cross,'ducted') == 1
        
        
            height = str2double(cellstr(Params{4,j}));
            width  = str2double(cellstr(Params{6,j}));
                
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));     % Datum related information
            
                      
        
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
        
        elseif strcmp(type_cross,'pipe')
            % For a honeycomb, we use the data of Lumley's iso-reduction curves
            dia         = str2double(cellstr(Params{4,j}));    
            datum_in    = str2double(cellstr(Params{8,j}));      % Datum related information
            datum_out   = str2double(cellstr(Params{9,j}));      % Datum related information
           
            % Assigning length array of component matrix
            components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
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
            f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
        end
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
    
    end


    
    
    if strcmp(type_comp,'Intake') == 1
        disp('Intake pipe for pipe section')
        % NOTE: For the intake, the length input is, infact, the angular
        % intake from the main pipe or riser
        % Input data
        intercept   = str2double(cellstr(Params{2,j}));
        dia         = str2double(cellstr(Params{4,j}));                        
        alpha       = str2double(cellstr(Params{10,j}));               % The divergence angle in degrees
        datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
        datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
        % Since we do not input the length, we calculate length as:
        length_intake = abs((datum_out - datum_in)/cos(deg2rad(alpha)));
        
        
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
        f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
            
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
        alpha       = str2double(cellstr(Params{10,j}));               % The divergence angle in degrees
        datum_in    = str2double(cellstr(Params{8,j}));               % Datum related information
        datum_out   = str2double(cellstr(Params{9,j}));               % Datum related information
        % Since we do not input the length, we calculate length as:
        length_intake = abs((datum_out - datum_in)/cos(deg2rad(alpha)));
        
        
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
        f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + length_intake;
        
    end
    
    if strcmp(type_comp,'Dynamometer') == 1
        disp('Dynamometer channel')
        % NOTE: We represent the dynamometer with a sphere, for
        % conservative design. We add the channel length it is located in
        % here as well.
        % Input data
        height      = str2double(cellstr(Params{4,j}));         % Duct height                       
        width       = str2double(cellstr(Params{6,j}));         % Duct width
        datum_in    = str2double(cellstr(Params{8,j}));         % Datum related information
        datum_out   = str2double(cellstr(Params{9,j}));         % Datum related information
        
        
        % Assigning length array of component matrix
        components(1,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(lengthcounter,lengthcounter + lengths(j),NN);
        % Assigning cross-sectional area array of component matrix
        components(2,1+(indexcounter*NN):NN*(indexcounter+1))   = height*width;
        % Assigning diameter related information
        components(3,1+(indexcounter*NN):NN*(indexcounter+1))   = height;
        components(4,1+(indexcounter*NN):NN*(indexcounter+1))   = width;
        % Assigning vertical position of the centreline
        components(5,1+(indexcounter*NN):NN*(indexcounter+1))   = linspace(datum_in,datum_out,NN);
                   
        % Volume Calculation
        volume(j) = trapz(components(2,1+(indexcounter*NN):NN*(indexcounter+1))) * (components(1,2+(indexcounter*NN)) - components(1,1+(indexcounter*NN))) ;    % Volume of straight pipe
                
        % f0byfi Calculation
        f0byfi(j) = (0.3^2)/components(2,1+(indexcounter*NN)) ;
        
        indexcounter = indexcounter+1;
        lengthcounter = lengthcounter + lengths(j);
        
    end
    counter   = counter+1;
    A_comp(j)    = components(2,((counter-1)*NN) + 1);
    
end





end

