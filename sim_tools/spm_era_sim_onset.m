function ons = spm_era_sim_onset(P)

ons     = 0;
switch P.distribution
    case 'normal'
        while any(ons<P.bounds(1)) || any(ons>P.bounds(2))
            ons = normrnd(P.mu,P.std,P.ntrials,1);
        end
end
% snap to nearest sample
ons = round(ons / P.dt) * P.dt;