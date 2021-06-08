function DepsDxprime = epstein_plesset_nondim_2(xprime, epsilon, xprime_dt, f , b, delta)



b_t = interp1(xprime_dt,b,xprime);         % Interpolate the data set (b,xprime_dt) at non-dimensional time xprime

%extra
delta_t = interp1(xprime_dt,delta,xprime);


DepsDxprime = -((1 - f + (delta_t/epsilon))/(1 + (2*delta_t/(3*epsilon)))) * ((xprime/epsilon) + (2*b_t));

end

