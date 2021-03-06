function alpha = SW_Diffusivity(T,uT,S,uS)
    % SW_Diffusivity    Thermal diffusivity of seawater
    %=========================================================================
    % USAGE:  alpha = SW_Diffusivity(T,uT,S,uS)
    %
    % DESCRIPTION:
    %   Thermal diffusivity of seawater at atmospheric pressure (0.1 MPa) using
    %   thermal conductivty, specific heat, and density correlations given in [1].
    %   Values at temperature higher than the normal boiling temperature
    %   are calculated at the saturation pressure.
    %
    % INPUT:
    %   T  = temperature
    %   uT = temperature unit
    %        'C'  : [degree Celsius] (ITS-90)
    %        'K'  : [Kelvin]
    %        'F'  : [degree Fahrenheit]
    %        'R'  : [Rankine]
    %   S  = salinity
    %   uS = salinity unit
    %        'ppt': [g/kg]  (reference-composition salinity)
    %        'ppm': [mg/kg] (in parts per million)
    %        'w'  : [kg/kg] (mass fraction)
    %        '%'  : [kg/kg] (in parts per hundred)
    %
    %   Note: T and S must have the same dimensions
    %
    % OUTPUT:
    %   alpha = thermal diffusivity [m^2/s]
    %
    %   Note: alpha will have the same dimensions as T and S
    %
    % VALIDITY: 0 < T < 180 C; 0 < S < 150 g/kg;
    %
    % ACCURACY: 3% (estimated at average value within the range)
    %
    % REVISION HISTORY:
    %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
    %               - Initial version
    %   2012-06-06: Karan H. Mistry (mistry@alum.mit.edu), MIT
    %               - Allow T,S input in various units
    %               - Allow T,S to be matrices of any size
    %   2015-07-01: Kishor G. Nayar (kgnayar@mit.edu)
    %               - Updated density and specific heat capacity functions
    %                 to include the input variable of pressure
    %   2016-04-10: Karan H. Mistry (mistry@alum.mit.edu), MIT
    %               - Allow T,S to be matrices of any size
    %
    % DISCLAIMER:
    %   This software is provided "as is" without warranty of any kind.
    %   See the file sw_copy.m for conditions of use and licence.
    %
    % REFERENCES:
    %   [1] M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair, Desalination
    %       and Water Treatment, 16, 354-380, 2010. (http://web.mit.edu/seawater/)
    %=========================================================================

    % RANGE CHECKING NOT REQUIRED SINCE OTHER PROPERTIES ARE CALLED

    %% BEGIN


    P0 = SW_Psat(T,'C',S,'ppt')/1E6;
    P0(find(T<100)) = 0.101325;

    uP = 'MPa';
    cp    = SW_SpcHeat(T,uT,S,uS,P0,uP);
    rho   = SW_Density(T,uT,S,uS,P0,uP);
    k     = SW_Conductivity(T,uT,S,uS);
    alpha = k./(rho.*cp);

end
