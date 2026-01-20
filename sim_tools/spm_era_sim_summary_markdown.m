
function spm_era_sim_summary_markdown(sim_num)

T = readtable(fullfile('tables', [sim_num '_results.xlsx']));
sim_levels = unique(T.SimLevel);

txtfile = fullfile('tables', [sim_num '_writeup.md']);
fid = fopen(txtfile, 'w');

fprintf(fid, '## Results Summary: Simulation %s\n\n', sim_num);

for i = 1:numel(sim_levels)

    this_sim = sim_levels{i};
    rows     = strcmp(T.SimLevel, this_sim);
    subt     = T(rows, :);

    fprintf(fid, '### %s\n\n', this_sim);
    fprintf(fid, 'For the %s condition, we compared four averaging methods across three modulation levels.\n', this_sim);

    % Identify best-performing row (contains '*')
    is_best = contains(subt.Improvement, '*');
    best_row = subt(is_best, :);

    % Describe best method
    fprintf(fid, 'The best-performing method was **%s** with %s modulation, achieving a mean %s of %s (SD = %s), representing a %s improvement relative to the unmodulated CA baseline.\n', ...
        best_row.Kappa{1}, ...
        lower(best_row.Alpha{1}), ...
        lower(strrep(best_row.Properties.VariableNames{4}, '_', ' ')), ...
        best_row.Mean_CI{1}, ...
        best_row.SD{1}, ...
        erase(best_row.Improvement{1}, '*'));

    % Describe others
    for j = 1:height(subt)
        if is_best(j), continue; end

        fprintf(fid, 'The %s method with %s modulation yielded %s (SD = %s), with a %s improvement.\n', ...
            subt.Kappa{j}, ...
            lower(subt.Alpha{j}), ...
            subt.Mean_CI{j}, ...
            subt.SD{j}, ...
            subt.Improvement{j});
    end

    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Writeup saved to %s\n', txtfile)
end

