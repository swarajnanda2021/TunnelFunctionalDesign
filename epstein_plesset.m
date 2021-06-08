function DrDt = epstein_plesset(t, r, f, time_tunnel, D_AB, tau, P_infty,rho_infty,C_S,sigma,M_a)


rho_infty_t     = interp1(time_tunnel,rho_infty,t);
P_infty_t       = interp1(time_tunnel,P_infty,t);
C_s_t           = interp1(time_tunnel,C_S,t);

rho_R           = (M_a/(8.314*298)) * (P_infty_t  + 2*sigma/r);
d               = C_s_t/rho_R;
% C_i_t           = f*d*rho_infty_t;
% C_i_t           = f*C_s_t;
% d = 0.02;

% DrDt = ((D_AB*(C_i_t - C_s_t))/(rho_infty_t + ( (2*tau)/(3*r) ))) * ( (1/r) + ( 1 / sqrt(pi*D_AB*t) ) );

DrDt = -D_AB*d* ( (1 - f + (tau/(r*rho_infty_t))) / (1 + (2*tau / (3*r*rho_infty_t))) ) * ( (1/r) + ( 1 / sqrt(pi*D_AB*t) ) );



end

