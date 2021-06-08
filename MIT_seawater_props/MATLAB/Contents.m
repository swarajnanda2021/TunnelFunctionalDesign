%========== SEAWATER THERMOPHYSICAL PROPERTIES LIBRARY =================
% Version 3.1.2 2016-07-08
% Version 3.1.1 2016-04-15
% Version 3.1   2016-03-03
% Version 3.0   2015-07-22
% Version 2.0   2012-06-06
% Version 1.0   2009-12-18
%
%
%           ***********************************************
%           *               SEAWATER PROPERTIES           *
%           *                                             *
%           *        John H. Lienhard V, Ph.D., P.E.      *
%           * Abdul Latif Jameel Professor of Water & Food*
%           *     Professor of Mechanical Engineering     *
%           *    Director, Center for Clean Water and     *
%           *       Clean Energy                          *
%           *    MASSACHUSETTS INSTITUTE OF TECHNOLOGY    *
%           *  77 Massachusetts Ave., Cambridge MA 02139  *
%           *              lienhard@mit.edu               *
%           *                                             *
%           *      Mostafa H. Sharqawy, Ph.D., P.Eng.     *
%           * Assistant Professor - School of Engineering *
%           *             University of Guelph            *
%           *             melsharq@uoguelph.ca            *
%           *                                             *
%           *           Kishor G. Nayar, S.M.             *
%           *            Ph.D. Candidate                  *
%           *     Department of Mechanical Engineering    *
%           *    Massachusetts Institute of  Technology   *
%           *                kgnayar@mit.edu              *   
%           *                                             *
%           ***********************************************
%
% LIST OF ROUTINES:
%
%     SW_COPY               Copyright and Licence file
%     SW_BPE                Boiling point elevation of seawater
% 	  SW_ChemPot_s          Chemical potential of salts in seawater
% 	  SW_ChemPot_w          Chemical potential of water in seawater
% 	  SW_SChemPot_s         Product of salinity and chemical potential of salt
%     SW_Conductivity       Thermal conductivity of seawater
%     SW_ConductivityP      Pressure-dependent thermal conductivity of seawater
%     SW_Density            Density of seawater
%     SW_Diffusivity        Thermal diffusivity of seawater
%     SW_Enthalpy           Specific enthalpy of seawater
%     SW_Entropy            Specific entropy of seawater
%     SW_FlowExergy         Specific flow exergy of seawater
%     SW_Gibbs              Specific Gibbs energy of seawater
%     SW_IntEnergy          Specific internal energy of seawater
%     SW_IsobExp            Isobaric expansivity of seawater
%     SW_IsothComp          Isothermal compressibility of seawater
%     SW_Kviscosity         Kinematic viscosity of seawater
%     SW_LatentHeat         Latent heat of vaporization of seawater
%     SW_OsmCoeff           Osmotic coefficient of seawater
%     SW_OsmPress           Osmotic pressure of seawater
%     SW_Prandtl            Prandtl of seawater
%     SW_Psat               Vapor pressure of seawater
%     SW_SpcHeat            Specific heat at constant pressure of seawater
%     SW_SurfaceTension     Surface tension of seawater
%     SW_Viscosity          Dynamic viscosity of seawater
%     SW_Volume             Specific volume of seawater
%
%
% REQUEST FROM AUTHORS:
% While our seawater properties library is freely distributed, we request 
% that users of this library please cite the following two journal articles 
% that comprehensively present how correlations in the library were developed:
%
%     a. K.G. Nayar, M.H. Sharqawy, L.D. Banchik and J.H. Lienhard V, 
%        Thermophysical properties of seawater: A review and new correlations
%        that include pressure dependence, Desalination, 390 (2016) 1-24. 
%        doi:10.1016/j.desal.2016.02.024.
%    
%     b. M.H. Sharqawy, J.H. Lienhard, S.M. Zubair, Thermophysical properties 
%        of seawater: a review of existing correlations and data, Desalination and
%        Water Treatment, 16 (2010) 354–380. doi:10.5004/dwt.2010.1079.
%
%
% Users using the following select functions developed by the authors are 
% encouraged to also cite the original sources:
%
%     a. Surface tension of seawater:
%        K.G. Nayar, D. Panchanathan, G.H. McKinley and J.H. Lienhard V, Surface 
%        Tension of Seawater, J. Phys. Chem. Ref. Data. 43 (2014) 043103. 
%        doi:10.1063/1.4899037.
%   
%    b. Thermal conductivity of seawater (atmospheric pressures):
%        M.H. Sharqawy, New correlations for seawater and pure water thermal 
%        conductivity at different temperatures and salinities, Desalination,
%        313 (2013) 97–104. doi:10.1016/j.desal.2012.12.010.
%     
%     c. Osmotic pressure of seawater (extrapolating to low salinities):
%        L.D. Banchik, M.H. Sharqawy, J.H. Lienhard V, Effectiveness-mass 
%        transfer units (ε-MTU) model of a reverse osmosis membrane mass exchanger,
%        J. Memb. Sci. 458 (2014) 189–198. doi:10.1016/j.memsci.2014.01.039.
%
%
% VERSION HISTORY:  
%     Version 3.1.2   2016-07-08
%       - Extended range of flow exergy function to S = 0 g/kg and updated description
%       - New function for product of salinity and chemical potential of salt
%
%     Version 3.1.1   2016-04-15
%       - Updated documentation of functions,range of chemical potential of
%         water function and surface tension, osmotic coefficient and references
%         to Nayar et al. (2016).
%
%     Version 3.1   2016-03-03                                       
%       - First public release version of revised package 
%       - New function for osmotic pressure (cacluated from osmotic 
%         coefficient) and pressure dependent thermal conductivity
%
%     Version 3.0   2015-07-22                                       
%       - Beta version
%       - Incorporated pressure dependence for density, specific heat 
%         capacity, enthalpy, entropy, Gibbs, chemical potential of 
%         salt and water in seawater, internal energy and specific volume
%       - New correlations for isothermal compressibility, isobaric
%         expansivity, Gibbs energy, chemical potential of salt and water
%         in seawater, osmotic coefficient and surface tension
%
%     Version 2.0   2012-06-06
%       - Allowed T and S arrays to be manipulated and unit conversion
%     Version 1.0   2009-12-18
%       - Original library dependent on temperature and salinity only
%
% ACKNOWLEDGEMENT:
% The development of this package of correlations and codes was primarily 
% sponsored through the Center for Clean Water and Clean Energy at MIT and 
% KFUPM during the period 2009 through 2016.  
%
% The authors wish to acknowledge Prof. Syed M. Zubair, Dr. Karan H. Mistry,
% Leonardo D. Banchik, Adam M. Wiener and Joanna Zhu for their contributions
% to developing various aspects of this library. Contributions to individual 
% functions have been acknowledged in the 'Revision history' of the function. 
%       
%=======================================================================


