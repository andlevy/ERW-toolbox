function u = spm_era_u(type,pst,ons,dur)

switch type
    case 'normal'
        u  = exp(-((pst - ons).^2) / (2 * dur^2));
    case 'gamma'
        a     = 2;
        u     = gampdf(pst-ons+dur,a,dur);
        u     = u./max(u);
    otherwise
        error('requested distribution not recognised')
end
