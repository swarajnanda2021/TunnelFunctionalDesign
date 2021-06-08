function dDdt = Levich(d,t,Dc,D_AB,M_a,sigma,PInfty)

U = ((-4.17e-1) + (1.12e-2 * (d*1e6)) + (1.42e-6 * (d*1e6)^2)) * 0.01;
dDdt = ( -3.02 * (8.314 * 298 * Dc * D_AB * (U^(1/3)) * (d^(1/3)) ) / (M_a * ((1.5 * d * PInfty) - (2*sigma))));

end

