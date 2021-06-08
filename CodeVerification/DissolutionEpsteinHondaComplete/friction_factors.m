function [ fric_facs, Q_inlet_comp ] = friction_factors( Q, Params, ComponentsMain, NN, nu, Q_pump, Q_pipe, TunnelSolidity, ValveDia )
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

k_friction = 0.01e-3;


for j=1:n_comp
    
    % Input component categorical information
    type_comp = cellstr(Params{1,j});
    type_cross= cellstr(Params{2,j});
    % Find out the inlet flow-rate
    Q_comp    = Q((counter-1)*NN + 1);
    Q_inlet_comp(j) = Q_comp;
    
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
                alphaby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
                n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
                friction = (f/(8*sin(alphaby2))) * (1 - (n_ar1)^-2); 
                
                expansion = 3.2 * (tan(alphaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                
                fric_facs(j) = friction + expansion;
                
            elseif dia_in > dia_out
                D_h = dia_in;
                w = Q_comp/(pi/4 * D_h^2);
                Re = w*D_h / nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                alphaby2 = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
                n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
                friction = (f/(8*sin(alphaby2))) * (1 - (n_ar1)^-2); 
                
                n0 = n_ar1^-1;
                alphap = 0.01745 * 2 * alphaby2;
                contraction = ( (-0.0125 * n0^4) + (0.0224 * n0^3) - (0.00723 * n0^2) + (0.00444 * n0) - 0.00745) * (alphap^3 - (2*pi*alphap^2) - (10*alphap));
                fric_facs(j) = contraction + friction;
            end
                
                
%                 D_h = dia_in;
%                 w = Q_comp/(pi/4 * D_h^2);
%                 Re = w*D_h/nu;
%                 n_ar1 = (max(dia_out,dia_in) / min(dia_out,dia_in))^2;
%                 alpha = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
%                 friction = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2) * (lengths(j)/D_h);
%                 expansion= 3.2 * tan(2*alpha) * tan(2*alpha)^(1/4) * (1 - (1/n_ar1))^2;
%                 fric_facs(j) = (friction + expansion)*1.2; % 1.2 is a non-uniformity correction
%                  
%             end
%                 % A generalised formula found in Idelchik is used here
%                 % (para 38 39 of Idelchik, ch-5, ed-3)
%                 D_h = dia_in;
%                 alpha = 2* rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
%                 x_bar = lengths(j)/D_h;
%                 x_tilde = log(1 + 2*x_bar * tan(deg2rad(alpha/2)))/(2*tan(deg2rad(alpha/2)));
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
%                     Phi = P(alpha,Re_0);
%                 elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 600000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(alpha,Re_0_new);
%                 elseif Re_0 < 0.5e5            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 50000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(alpha,Re_0_new);
%                 end
%                 zeta_un = Phi * (1 - (1/n_ar1))^1.92;
%                  
%                 a = 0.924 / (1 + (1.3 * 1e-5 * alpha^pi));
%                 b = (0.3 + (1.55  * 1.1^(-alpha)))/((1 + (1.03 * 1e-8 * l0_bar^7.5)));
%                 c = 1.05 / (1 + (2.3 * 1e-62 * Re_0^11));
%                 zeta_non = 0.044 * (0.345*alpha)^a * (1 - (0.2*n_ar1 + 0.8)^(-3.82)) * (0.154 * l0_bar)^b * ((2.31 * 1e-6 * Re_0) + (1.1 * 10e-30 * Re_0^5.62))^c;
%                 
%                 lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                 zeta_fr = (lambda / (8 * sin(deg2rad(alpha/ 2) ) ) ) * (1 - (1/n_ar1^2) ); 
%                 zeta_frprime = (1 + (0.5/(1.5 ^ x_tilde)))*zeta_fr;
%                 
%                 % Calculation of friction factor
%                 zeta = zeta_frprime + zeta_non + zeta_un ;
%                 % Assigning to fric_facs matrix
%                 fric_facs(j) = zeta;
%             end
%Method2
% Calculating the diffuser angle (theta/2)
%                 alpha = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j));
%                 n0      = ((min(dia_in, dia_out)/max(dia_in, dia_out)) ^ 2);
%                 D_h = dia_in;
%                 w = Q_comp/(pi/4 * D_h^2);
%                 Re = w*D_h/nu;
% 
%                 if dia_in < dia_out
%                     zeta  = (1 - (n0))^2 * 2.6 * sin(alpha);
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h ;%* 1.5; % 1.5 is a correction factor to account for a non-uniform profile
%                 
%                 elseif dia_in > dia_out
%                     % If the pipe-section is converging (can be expected),
%                     % we apply the formula in diagram 5-23 of Idelchik.
%                     
% %                     alpha   = rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
% %                     alpha_r = 0.01745*alpha;
% %                     zeta_converging    = (((-0.0125*n0^4) + (0.0224*n0^3) + (-0.00723*n0^2) + (0.00444*n0) + (-0.00745)) * ((alpha_r^3) - (10*alpha_r) - (2*pi * alpha_r^2)))  + zeta;
%                     if alpha <= 45
%                         zeta = 0.5*(1 - (n0))^0.75 * 1.6 * sin(alpha);
%                     elseif alpha > 45
%                         zeta = 0.5*(1 - (n0))^0.75 * sqrt( sin(alpha));
%                     end
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h ;%* 1.5; % 1.5 is a correction factor to account for a non-uniform profile
%                 end
%                 
%             end
 


            
                        
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
                D_h = (2 * b*h) / (b + h) * 0.64;            
            
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                        
                % Calculation of friction facgtor from Haaland's explicit
                % formula
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                
                fric_facs(j) = f*lengths(j)/D_h;
                
            elseif (height_in < height_out) && (width_in < width_out)     % Axisymmetric diffuser
                
                h = height_in;
                b = width_in;
                D_h = (2 * b*h) / (b + h) * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                alphaby2 = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
                betaby2 = atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j));
                n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                
                friction = (f/16) * (1 - (n_ar1)^-2) * (sin(betaby2)^(-1) + sin(alphaby2)^(-1)); 
                
                if (2*rad2deg(alphaby2) > 4 && 2*rad2deg(alphaby2) < 12)
                    K = 0.66 + (0.12*2*rad2deg(alphaby2)); 
                    expansionalpha = 3.2 * K * (tan(alphaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(alphaby2) > 12 && 2*rad2deg(alphaby2) < 30)
                    K = 3.3 - (0.03*2*rad2deg(alphaby2));
                    expansionalpha = 3.2 * K * (tan(alphaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionalpha = 0;
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
                
                expansion = (expansionalpha + expansionbeta) / 2;
                fric_facs(j) = friction + expansion;
            
            elseif (height_in == height_out && width_in < width_out) || (height_in < height_out && width_in == width_out)   % Plane diffusers
                
                h = height_in;
                b = width_in;
                D_h = (2 * b*h) / (b + h) * 0.64;            
                
                w = Q_comp/(height_in*width_in);
                Re = w*D_h/nu;
                
                f = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                alphaby2 = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
                betaby2 = atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j));
                if rad2deg(alphaby2) > rad2deg(betaby2)
                    n_ar1 = (max(height_out,height_in) / min(height_out,height_in))^2;
                    friction = (f/(4*sin(alphaby2))) * ( (b/h * (1 - (n_ar1)^-1)) + (0.5*(1 - (n_ar1)^-2)) ); 
                else
                    n_ar1 = (max(width_out,width_in) / min(width_out,width_in))^2;
                    friction = (f/(4*sin(betaby2))) * ((b/h * (1 - (n_ar1)^-1)) + (0.5*(1 - (n_ar1)^-2))); 
                end
                
                
                
                if (2*rad2deg(alphaby2) > 4 && 2*rad2deg(alphaby2) < 12)
                    K = 2 - (0.03*2*rad2deg(alphaby2)); 
                    expansionalpha = 3.2 * K * (tan(alphaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                elseif (2*rad2deg(alphaby2) > 12 && 2*rad2deg(alphaby2) < 20)
                    K = 2 - (0.04*2*rad2deg(alphaby2));
                    expansionalpha = 3.2 * K * (tan(alphaby2))^1.25 * (1 - (n_ar1)^-1)^2;
                else 
                    expansionalpha = 0;
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
                
                
              
                
                expansion = (expansionalpha + expansionbeta) / 2;
                
                fric_facs(j) = friction + expansion;
            
            end
            
        end
% Method 3) Idelchik's own formula
% 
%                 D_h = (2 * b*h) / (b + h);            
%                 w = Q_comp/(height_in*width_in);
%                 Re = w*D_h/nu;
%                 alpha = atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j));
%                 friction = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2) * (lengths(j)/D_h);
%                 expansion= 3.2 * tan(2*alpha) * tan(2*alpha)^(1/4) * (1 - (1/n_ar1))^2;
%                 fric_facs(j) = (friction + expansion)*1.2; % 1.2 is a non-uniformity correction
                
                
% Method 1) A generalised formula found in Idelchik is used here
% (para 40 41 of Idelchik, ch-5, ed-3)
%                 D_h = (2 * height_in * width_in)/(height_in + width_in);
%                 
%                 alpha = 2*rad2deg(atan(((max(height_in,height_out) - min(height_in,height_out))/2) / lengths(j)));
%                 beta = 2*rad2deg(atan(((max(width_in,width_out) - min(width_in,width_out))/2) / lengths(j)));
%                 
%                 a0 = height_in/D_h;
%                 b0 = width_in/D_h;
%                 x_bar = lengths(j)/D_h;
%                 x_tilde = (1/(4*tan(deg2rad(alpha/2))))*log(((4 * x_bar^2 * tan(deg2rad(alpha/2))^2) + (2 * x_bar * (a0 + b0) * tan(deg2rad(alpha/2))) + (a0*b0))/(a0*b0));
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
%                     Phi = P(alpha,Re_0);
%                 elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 600000;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(alpha,Re_0_new);
%                 elseif Re_0 < 0.5e5            % Limit to value of 600000 for Re above this value
%                     Re_0_new = 0.5e5;
%                     data_Phi = csvread('diagram5p2_Idelchik.txt');
%                     P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
%                     Phi = P(alpha,Re_0_new);
%                 end
%                 zeta_un = Phi * (1 - (1/n_ar1))^1.76;
%                  
%                 s = 1.06 / (1 + (2.82 * 1e-3 * alpha^2.24));
%                 t = 0.73 / (1 + (4.31 * 1e-6 * l0_bar^7.31));
%                 u = 1.06 / (1 + (1.1 * 1e-30 * Re_0^5.62));
%                 zeta_non = 0.024 * (0.625*alpha)^s * (1 - (2.81*n_ar1 - 1.81)^(-1.04)) * (0.303 * l0_bar)^t * ((4.8 * 1e-7 * Re_0) + 1.8)^u;
%                 
%                 lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                 zeta_fr = (lambda / 16 ) * (1 - (1 / n_ar1^2) ) * ((1/sin(deg2rad(alpha/2))) + (1/sin(deg2rad(beta/2))) );
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
%                     zeta  = (1 - (n0))^2 * 2.6 * sin(alpha);    
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h;
%                 
%                 elseif height_in > height_out
%                     % If the pipe-section is converging (can be expected),
%                     % we apply the formula in diagram 5-23 of Idelchik.
%                     
% %                     alpha   = rad2deg(atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)));
% %                     alpha_r = 0.01745*alpha;
% %                     zeta_converging    = (((-0.0125*n0^4) + (0.0224*n0^3) + (-0.00723*n0^2) + (0.00444*n0) + (-0.00745)) * ((alpha_r^3) - (10*alpha_r) - (2*pi * alpha_r^2)))  + zeta;
%                     if alpha <= 45
%                         zeta = 0.5*(1 - (n0))^0.75 * 1.6 * sin(alpha);
%                     elseif alpha > 45
%                         zeta = 0.5*(1 - (n0))^0.75 * sqrt( sin(alpha));
%                     end
%                     lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
%                     fric_facs(j) = zeta + lambda*lengths(j)/D_h;
%             end
%             end

%             
%         end
        
    end
    
    if strcmp(type_comp,'Contraction') == 1
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
            
            % Parameters for the calculation
            h = height_in;
            b = width_in;
            D_h = (2 * b*h) / (b + h);            
            w = Q_comp/(height_in*width_in);
            Re = w*D_h/nu;
                    
            Re_sh = (Re/4) * (1 + (b/h));
            Re_long = (Re/4) * ((1 + (b/h))/(b/h));
            lambda_sh = (3.6 * log(Re_sh)   -   2)^(-2);
            lambda_long = (3.6 * log(Re_long)   -   2)^(-2);
            lambda = 4*((b/h)/(1+(b/h)))*(1 + (lambda_sh/lambda_long) * (h/b)) * lambda_long;
            alpha = atan(((height_out - height_in) / 2) / lengths(j));
            % Calculation of friction factor
            f = (lambda / (8 * abs(sin(alpha/2)))) * (1 - 1/((height_in/height_out) * (width_in/width_out)));
            
            fric_facs(j) = f;
            
                        
        elseif strcmp(type_cross,'pipe') == 1    
            
            disp('Conical contraction')
            % Input data
            dia_in      = str2double(cellstr(Params{4,j}));                                          % Contraction area-ratio of 2 is employed
            dia_out     = str2double(cellstr(Params{5,j}));
    
            % A generalised formula found in Idelchik is used here
            % (para 38 39 of Idelchik, ch-5, ed-3)
            D_h = dia_in;
            alpha = atan(((max(dia_in,dia_out) - min(dia_in,dia_out))/2) / lengths(j)); 
            Re_0 = dia_in * (Q_comp/(pi * dia_in^2/ 4))/nu;
               
            % Using look-up table for Phi value
            if Re_0 < 600000
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(alpha),Re_0);
            elseif Re_0 > 600000            % Limit to value of 600000 for Re above this value
                Re_0_new = 600000;
                data_Phi = csvread('diagram5p2_Idelchik.txt');
                P = scatteredInterpolant(data_Phi(:,1),data_Phi(:,2),data_Phi(:,3));
                Phi = P(rad2deg(alpha),Re_0_new);
            end
                
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            zeta_fr = (lambda / (8 * sin(alpha / 2) ) ) * (1 - (1 / ( (max(dia_in, dia_out)/min(dia_in, dia_out)) ^ 2) ) ); 
            
            % Calculation of friction factor
            zeta = zeta_fr ;
            % Assigning to fric_facs matrix
            fric_facs(j) = zeta;
            
        
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
            
            D_h = 2*height*width/(height+width);
            w = Q_comp/(height*width);
            Re = w * D_h/nu;
            k_re = 0.8 + (4.02e4/Re);
            
            if strcmp(corner_mod,'thin_inner_sharp_normal') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.45 * k_re) +  lambda;
                
            elseif strcmp(corner_mod,'thin_inner_sharp_reduced') == 1
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.45 * k_re) +  lambda;
                            
            elseif strcmp(corner_mod,'thin_inner_bevelled') == 1
                % Bevelling is t = 0.25*width_inlet
                
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta = (0.36 * k_re) + (1.28 * lambda);
                
            elseif strcmp(corner_mod,'thin_rounded') == 1
                % Rounding is r = 0.2*width_inlet
                lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
                zeta_fr = (1 + (1.57*0.2)) * lambda;
                zeta = (k_re * 0.14) + zeta_fr;         
            end

            fric_facs(j) = zeta;
            
            
        elseif strcmp(type_cross,'arbitrary') == 1
            disp('Arbitrary bend')
            dia         = str2double(cellstr(Params{4,j}));
            alpha       = str2double(cellstr(Params{5,j}));                                       % The bending angle in degrees
            R0          = str2double(cellstr(Params{10,j}));
            
            % Modify lengths(j)
            lengths(j) = (deg2rad(alpha)*R0);
            
            %Diagram 6-1 of Idelchik
            %A1
            if alpha <= 70
                A1 = 0.9*sin(deg2rad(alpha));
            elseif alpha == 90
                A1 = 1;
            elseif alpha >= 100
                A1 = 0.7 + (0.35*alpha/90);
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
            lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);

            zeta_loc = A1*B1*C1;
            zeta_fr  = 0.0175 * alpha * lambda * R0/dia;
            
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
            
            D_h = (2 * b*h) / (b + h);            
            w = Q_comp/(height_out*width_out);
            
            % Calculation of friction factor
            Re_0 = D_h*w/nu;               
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            % Calculation of friction factor
            zeta = lambda * lengths(j)/D_h;
            % Assigning to fric_facs matrix
            fric_facs(j) = zeta;

                    
                      
            
                       
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
            
            % Calculation of friction factor
            Re_0 = D_h*w/nu;               
            lambda = (-1.8 * log10(6.9/Re_0 + ((k_friction/D_h)/3.7)^1.11))^(-2);
            % Calculation of friction factor
            zeta = lambda * lengths(j)/D_h;
            % Assigning to fric_facs matrix
            fric_facs(j) = zeta;        
                       
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
                disp(['Mesh Reynolds number of ' num2str(Re_D) ' is lower than 1000. Increase.'])
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


    
    
    if strcmp(type_comp,'Intake') == 1
        disp('Intake pipe for pipe section')
        % NOTE: For the intake, the length input is, infact, the angular
        % intake from the main pipe or riser
        % Input data
        intercept   = str2double(cellstr(Params{2,j}));
        dia         = str2double(cellstr(Params{4,j}));                        
        alpha       = str2double(cellstr(Params{10,j}));               % The divergence angle in degrees
        
        
        
        % Friction calculations
        A_pipe      =   pi/4 * dia^2 ;
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
        
        zeta_cs     = A_prime * (1 + ((Q_pipe / (Q_pump - Q_pipe))*(A_pipe / A_tunnel))^2 - (2 * ((Q_pipe / (Q_pump - Q_pipe))*(A_pipe / A_tunnel)) * cos(deg2rad(alpha))));
        zeta_cst    = tau_st * (Q_pipe/(Q_pump - Q_pipe))^2;
        
        zeta_s  = zeta_cs / ((Q_pipe/Q_pump)*(A_tunnel/A_pipe))^2;
        zeta_st = zeta_cst / ((1 - (Q_pipe/Q_pump))^2  * (A_tunnel/A_pipe)); 
        
        zeta = ((Q_pipe/Q_pump) * zeta_s) + (((Q_pump - Q_pipe)/Q_pump) * zeta_st);
                  
        % Add length of pipe as friction
        D_h = dia;
        w = Q_pipe/(pi/4 * D_h^2);
        Re = w*D_h/nu;
        lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
               
        % Calculation of friction factor
            
        f = lambda / D_h;
        fric_facs(j) = f*lengths(j)+zeta;
        
        
    end
    
    
    if strcmp(type_comp,'Remerge') == 1
        disp('Merging connection for pipe')
        
        % NOTE: For the intake, the length input is, infact, the angular
        % intake from the main pipe or riser
        % Input data
        remerge     = str2double(cellstr(Params{2,j}));
        dia         = str2double(cellstr(Params{4,j}));                        
        
        % Friction calculations
        A_pipe      =   pi/4 * dia^2; 
        [x, index]  =   unique(ComponentsMain(1,:));
        A_tunnel    =   interp1(x,ComponentsMain(2,index),35);    % Finding the cross-section of the main tunnel at the location of the intake
        
      
        
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
        
        zeta_cs     = A * (1 + ((Q_pipe / (Q_pump - Q_pipe))*(A_pipe / A_tunnel))^2 - (2 * (1 - (Q_pipe / (Q_pump - Q_pipe))))    - (1.7 * ((Q_pipe / (Q_pump - Q_pipe)) * (A_pipe / A_tunnel)^2)) );
        zeta_cst    = -(1.7 * ((Q_pipe / (Q_pump - Q_pipe)) * (A_pipe / A_tunnel)^2));
        
        zeta_s  = zeta_cs / ((Q_pipe/Q_pump)*(A_tunnel/A_pipe))^2;
        zeta_st = zeta_cst / (1 - (Q_pipe/Q_pump)^2); 
        
        zeta = ((Q_pipe/Q_pump) * zeta_s) + (((Q_pump - Q_pipe)/Q_pump) * zeta_st);
                  
        % Add length of pipe as friction
        D_h = dia;
        w = Q_pipe/(pi/4 * D_h^2);
        Re = w*D_h/nu;
        lambda = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
               
        % Calculation of friction factor
            
        f = lambda / D_h;
        fric_facs(j) = f*lengths(j)+zeta;
                
        
        
    end
    
    if strcmp(type_comp,'Dynamometer') == 1
        disp('Dynamometer channel')
        % NOTE: We represent the dynamometer with a sphere, for
        % conservative design. We add the channel length it is located in
        % here as well.
        % Input data
        sphere_dia  = str2double(cellstr(Params{2,j}));         % Sphere diameter in the tunnel
        height      = str2double(cellstr(Params{4,j}));         % Duct height                       
        width       = str2double(cellstr(Params{6,j}));         % Duct width
        
        %Friction calculations
        
        % Friction by the channel
        w   = Q_comp/(height*width);
        D_h = 2 * height * width/ (height+width);
        Re  = w*D_h/nu;
        
        lambda      = (-1.8 * log10(6.9/Re + ((k_friction/D_h)/3.7)^1.11))^(-2);
        f_channel   = lambda*lengths(j)/D_h;  
        
        % Friction by the sphere
        Re_sph   = sphere_dia * w/nu;
        
        % From the correlation of Morrison:
        Cd       = (24/Re_sph) + ((2.6*Re_sph/5)/(1 + (Re_sph/5)^1.52)) + ((0.411 * (Re_sph/ 2.63e5)^(-7.94) )/(1 + ((Re_sph/ 2.63e5))^-8)) + ((0.25*Re_sph/1e6)/(1 + (Re_sph/1e6)));
        
        % Adding individual contributions
        zeta     = Cd + f_channel;
        fric_facs(j) = zeta;
            
        
    end    
    
    counter = counter+1;
end



if TunnelSolidity > 0
    theta    = asin(TunnelSolidity * 0.3 / ValveDia); % Opening angle of the valve
    % Display warning for lightly opened throttling valve
    if rad2deg(theta) < 30
        disp(['The valve opening angle ' num2str(theta) ' is less than 30 degrees. Friction factor will not be accurately calculated.'])
    end
    Re      = ((Q_pump - Q_pipe)/(0.3^2)) * 0.3 / nu;
    Dd_bar  = ValveDia/0.3; 
    A       = abs(60 * ((1 + (0.5 * Dd_bar * (1 + sin(theta)))) / (1 - (Dd_bar^2 * sin(theta))^2 )));
    zeta_qu = (((1.56) / (1 - (Dd_bar*sin(theta)))) - 1)^2;
    
    zeta_throttling = abs((A/Re) + ((1- (50/Re)) * zeta_qu));
    
    fric_facs = [fric_facs zeta_throttling]; 
        
end


end

