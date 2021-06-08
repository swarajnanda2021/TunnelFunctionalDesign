function DrDt = epstein_plesset(r, t, f, d_tunnel, time_tunnel, D_AB, tau, rho_f)


d   = interp1(time_tunnel,d_tunnel,t);


DrDt = - D_AB * d * ((1 - f + (tau/(r * rho_f) ) )/(1 + (2*tau/(3 * r * rho_f) ))) * ((1/r) + (1/sqrt(pi*D_AB*t)));

end

