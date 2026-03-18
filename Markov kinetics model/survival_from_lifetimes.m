function [t, S] = survival_from_lifetimes(lifetimes_s)
%SURVIVAL_FROM_LIFETIMES Survival probability vs time (Kaplan–Meier style)
%   lifetimes_s : Nx1 lifetimes in seconds
%   t           : time vector (s)
%   S           : survival probability S(t)

Tmax = 0.5*3600;              % 3 hours in seconds

% censor everything beyond Tmax
lifetimes_s = lifetimes_s(:);
is_event = lifetimes_s < Tmax;      % true detachment
lifetimes_s(~is_event) = Tmax;      % censor at Tmax

% unique event times (including Tmax)
t = unique(lifetimes_s);
t = sort(t);

N = numel(lifetimes_s);
S = ones(size(t));

n_at_risk = N;
S_curr = 1;

for i = 1:numel(t)
    ti = t(i);

    % number of detachments at this time
    d = sum(is_event & lifetimes_s == ti);

    % number censored at this time
    c = sum(~is_event & lifetimes_s == ti);

    % Kaplan–Meier update
    if n_at_risk > 0
        S_curr = S_curr * (1 - d/n_at_risk);
    end
    S(i) = S_curr;

    % update risk set
    n_at_risk = n_at_risk - d - c;
end

end
