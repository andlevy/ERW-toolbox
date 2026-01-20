function ERA = spm_era_load(era)


if ischar(era) && size(era,1) == 1
    tmp    = load(era,'ERA');
    ERA    = tmp.ERA;
elseif iscell(era) && ischar(era{1})
    for i = 1:numel(era)
        tmp    = load(era{i},'ERA');
        ERA{i} = tmp.ERA;
    end
elseif iscell(era) && isstruct(era{1})
    ERA = era;
elseif isstruct(era)
    ERA = era;
end

