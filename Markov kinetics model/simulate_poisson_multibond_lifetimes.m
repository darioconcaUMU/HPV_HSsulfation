function lifetimes_s = simulate_poisson_multibond_lifetimes(mu, koff_per_bond, varargin)
%SIMULATE_POISSON_MULTIBOND_LIFETIMES Simulate detachment lifetimes for particles
%   - 100 particles
%   - #bonds ~ Poisson(mu), conditioned on >=1 (no 0-bond particles)
%   - Cap bonds at 5 (values >5 are set to 5)
%   - koff depends on bond count (1..5): koff_per_bond(n) in s^-1
%   - Particle is "sampled" every dt=15 s; lifetime is first sample time at/after detachment.
%
% Usage:
%   koff = [0.20 0.10 0.06 0.04 0.03]; % s^-1 for 1..5 bonds
%   lifetimes = simulate_poisson_multibond_lifetimes(2.5, koff);
%
% Optional name-value:
%   'N'   : number of particles (default 100)
%   'dt'  : sampling interval in seconds (default 15)
%   'rng' : RNG seed (default [])
%
% Output:
%   lifetimes_s : Nx1 vector of observed lifetimes (seconds), quantized to dt.

% ----------------- parse inputs -----------------
p = inputParser;
p.addRequired('mu', @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addRequired('koff_per_bond', @(x) isnumeric(x) && isvector(x) && numel(x) == 5 && all(x > 0));
p.addParameter('N', 1000, @(x) isnumeric(x) && isscalar(x) && x == round(x) && x > 0);
p.addParameter('dt', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('rng', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.parse(mu, koff_per_bond, varargin{:});

N  = p.Results.N;
dt = p.Results.dt;
seed = p.Results.rng;

if ~isempty(seed)
    rng(seed);
end

koff_per_bond = koff_per_bond(:); % ensure column

% ----------------- draw bond counts (Poisson conditioned on >=1) -----------------
bonds = zeros(N,1);
for i = 1:N
    k = 0;
    % rejection sampling until k>=1
    while k == 0
        k = poissrnd(mu);
    end
    % cap at 5
    if k > 5
        k = 5;
    end
    bonds(i) = k;
end

% ----------------- simulate continuous detachment times -----------------
% For a memoryless detachment with rate koff, lifetime T ~ Exp(koff):
%   T = -log(U)/koff
U = rand(N,1);
koff = koff_per_bond(bonds);          % per-particle koff (s^-1)
t_detach_cont = -log(U) ./ koff;      % continuous detachment time (s)

% ----------------- sample every dt: observed lifetime is first sample >= detachment -----------------
% If we "check" at times dt, 2dt, 3dt, ... then observed time is ceil(T/dt)*dt
lifetimes_s = ceil(t_detach_cont ./ dt) .* dt;

% ----------------- (optional) quick display -----------------
% Uncomment if you want a quick summary printed:
% fprintf('Mean bonds: %.2f\n', mean(bonds));
% fprintf('Mean lifetime (observed): %.2f s (%.2f min)\n', mean(lifetimes_s), mean(lifetimes_s)/60);

end
