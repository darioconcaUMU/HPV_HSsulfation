% Parameters for simulation
% 
% 
% 
R = sqrt(1.818655106633252e+03/pi)/1000;
mu = 1.1*1e3*pi .* R.^2;
koff = 1./[2.50E+02	3.59E+05	6.89E+08	1.49E+12	3.42E+15]; % s^-1 for 1..5 bonds
lifetimes = simulate_poisson_multibond_lifetimes(mu, koff, 'rng', round(rand*100));
[t1, S1] = survival_from_lifetimes(lifetimes);
R = sqrt(1.818655106633252e+03/pi)/1000;
mu = 1.3*1e3*pi .* R.^2;
koff = 1./[6.67E+01	4.05E+04	3.27E+07	2.98E+10	2.89E+13]; % s^-1 for 1..5 bonds
lifetimes = simulate_poisson_multibond_lifetimes(mu, koff, 'rng', round(rand*100));
[t2, S2] = survival_from_lifetimes(lifetimes);
R = sqrt(1.818655106633252e+03/pi)/1000;
mu = 1.5*1e3*pi .* R.^2;
koff = 1./[8.33E+00  1.89E+02  5.58E+03  1.86E+05  6.58E+06]; % s^-1 for 1..5 bonds
lifetimes = simulate_poisson_multibond_lifetimes(mu, koff, 'rng', round(rand*100));
[t3, S3] = survival_from_lifetimes(lifetimes);
