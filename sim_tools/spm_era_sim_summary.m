function SUM = spm_era_sim_summary(sim_num,metric,savedir)

load(fullfile(savedir,['SIM' sim_num '.mat']),'SIM')

nvar    = size(SIM,1);
nsim    = size(SIM,2);
nkappa  = size(SIM{1,1},2);
nalpha  = size(SIM{1,1},1);

E    = [];
Yu   = [];
Yw   = [];
Yg   = [];
GT   = [];
G    = [];

g = 1;
for Iv = 1:nvar
    for Ia = 1:nalpha
        for Ik = 1:nkappa
            for Is = 1:nsim

                if isfield(SIM{Iv,Is}(Ia,Ik).E.(metric),'wm')
                    E(Is,Ia,Ik,Iv)     = SIM{Iv,Is}(Ia,Ik).E.(metric).wm;
                    Yw(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).Aw.wm;
                else
                    E(Is,Ia,Ik,Iv)     = SIM{Iv,Is}(Ia,Ik).E.(metric).mn;
                    Yw(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).Aw.mn;
                end

                if Ia == 1 && Ik == 1

                    Yg(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).Ag.mn;
                    Ye(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).Ae.mn;
                    Yu(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).Au.mn;
                    GT(:,Is,Ia,Ik,Iv)  = SIM{Iv,Is}(Ia,Ik).GT2;

                end

                G(Is,Ia,Ik,Iv)     = g;

            end
            g = g+1;
        end
    end
end

% Calculate improvement
%--------------------------------------------------------------------------
m             = 10000;       % Number of bootstrap replicates
alpha         = 0.05;        % Confidence level
IV            = nan(nalpha,nkappa,nvar);
CI5_IV        = nan(nalpha,nkappa,nvar);
CI95_IV       = nan(nalpha,nkappa,nvar);

for Iv = 1:nvar
    base = squeeze(E(:,1,1,Iv));  % Reference: CA

    for Ia = 1:nalpha
        for Ik = 1:nkappa
            if Ia == 1 && Ik == 1
                IV(Ia,Ik,Iv)      = 0;
                CI5_IV(Ia,Ik,Iv)  = 0;
                CI95_IV(Ia,Ik,Iv) = 0;
                continue
            end

            comp = squeeze(E(:,Ia,Ik,Iv));
            imp  = 100 * (base - comp) ./ base;

            % Mean improvement
            IV(Ia,Ik,Iv) = mean(imp);

            % Bootstrap CIs
            bootstat           = bootstrp(m, @mean, imp);
            CI5_IV(Ia,Ik,Iv)   = prctile(bootstat, 100 * alpha/2);        % Lower bound
            CI95_IV(Ia,Ik,Iv)  = prctile(bootstat, 100 * (1 - alpha/2));  % Upper bound
        end
    end
end


[~,i]     = min(sum(abs(E - 0),[2 3 4]));
m         = 10000;                      % bootstrap replicates
bootstat  = bootstrp(m,@mean,E);        % sampling distribution
alpha     = 0.05;
ci_low    = prctile(bootstat,100*alpha/2);      % 0.5 %
ci_high   = prctile(bootstat,100*(1-alpha/2));  % 99.5 %
mean_x    = spm_vec(mean(E));
lower_err = reshape(mean_x - ci_low',nalpha,nkappa,nvar);           % asymmetric distances
upper_err = reshape(ci_high' - mean_x,nalpha,nkappa,nvar);

switch metric
    case 'nrmse'
        metric_label = 'nRMSE';
end

% Package
%--------------------------------------------------------------------------
SUM.E       = E;
SUM.Yw      = Yw;
SUM.Yg      = Yg;
SUM.Ye      = Ye;
SUM.Yu      = Yu;
SUM.GT      = GT;
SUM.G       = G;
SUM.ME      = squeeze(mean(E,1));
SUM.MD      = squeeze(median(E,1));
SUM.SD      = squeeze(std(E,[],1));
SUM.SE      = (squeeze(std(E,[],1))./sqrt(nsim));
SUM.CI5     = lower_err;
SUM.CI95    = upper_err;
SUM.IV      = IV;
SUM.sumnum  = sim_num;
SUM.metric  = metric_label;
SUM.Isim    = i;
SUM.pst     = SIM{1, 1}(1, 1).pst;

% Tabulate
%--------------------------------------------------------------------------
kapu =  char(0x03BA);
sigu = char(0x03C3);

switch sim_num
    case 'noise_nrmse'
        level_names = {[sigu ': 25 ms'],[sigu ': 75 ms'],[sigu ': 125 ms']};
    case 'noise_less'
        level_names = {[sigu ': 25 ms'],[sigu ': 75 ms'],[sigu ': 125 ms']};
    case 'one'
        level_names = {[sigu ': 25 ms'],[sigu ': 75 ms'],[sigu ': 125 ms']};
    case 'two'
        level_names = {'Skew(0)','Skew(-1)','Skew(1)'};
    case 'three'
        level_names = {'linear','quadratic','cubic'};
    case 'four'
        level_names = {'linear','quadratic','cubic'};
end

avg_names   = {'CA',['WA(' kapu ':1)'],['WA(' kapu ':2)'],['WA(' kapu ':3)']};
% mod_names   = {'Unmodulated','Partially Modulated', 'Fully Modulated'};
row_idx     = 1;
TAB         = struct();
for n = 1:size(SUM.ME,3)

    best      = 0;
    best_idx  = 1;

    for l = 1%:size(SUM.ME,1)
        for m = 1:size(SUM.ME,2)

            TAB.SimLevel{row_idx,1} = level_names{n};
            TAB.Kappa{row_idx,1}    = avg_names{m};
            % TAB.Alpha{row_idx,1}    = mod_names{l};

            % Summary statistics
            mu      = SUM.ME(l,m,n);
            lo      = SUM.ME(l,m,n)-SUM.CI5(l,m,n);
            hi      = SUM.ME(l,m,n)+SUM.CI95(l,m,n);
            sigma   = SUM.SD(l,m,n);
            improv  = SUM.IV(l,m,n);

            mu_imp  = IV(l,m,n);
            lo_imp  = CI5_IV(l,m,n);
            hi_imp  = CI95_IV(l,m,n);


            if improv > best
                best     = improv;
                best_idx = row_idx;
            end

            if improv == 0
                improv_str = '--';
            elseif improv < 0
                improv_str = sprintf('%.2f%% [%.2f%%, %.2f%%] Higher', mu_imp,lo_imp,hi_imp);
            else
                improv_str = sprintf('%.2f%% [%.2f%%, %.2f%%] Lower', mu_imp,lo_imp,hi_imp);
            end

            TAB.Mean_CI{row_idx,1}          = sprintf('%.3f [%.3f, %.3f]', mu, lo, hi);
            TAB.SD{row_idx,1}               = sprintf('%.3f', sigma);
            TAB.Improvement_CI{row_idx,1}   = improv_str;

            % % Store for ranking
            % ME_for_rank(row_idx) = mu;
            % simlevel_id(row_idx) = n;

            row_idx = row_idx +1;
        end
    end

    TAB.Improvement_CI{best_idx,1} = [TAB.Improvement_CI{best_idx,1} '*'];
end

TAB = struct2table(TAB);
filePath = fullfile(savedir,[sim_num '_results.xlsx']);

% Delete existing table file if it exists
if exist(filePath, 'file')
    delete(filePath);
end

% Write new table from scratch
writetable(TAB, filePath, 'Sheet', 1);


pause(3)
%% Format table
%==========================================================================
% Requires Excel on Windows
xlsfile = fullfile(savedir, [sim_num '_results.xlsx']);

% Open Excel
ex = actxserver('Excel.Application');
ex.DisplayAlerts = false;
wb = ex.Workbooks.Open(xlsfile);
ws = wb.Worksheets.Item(1);

try
    % Insert a new first row
    ws.Rows.Item(1).Insert;

    % Set base font for entire used range
    ur           = ws.UsedRange;
    ur.WrapText  = true;
    ur.Font.Name = 'Palatino Linotype';
    ur.Font.Size = 12;
    ur.HorizontalAlignment   = -4108;           % xlCenter
    ur.VerticalAlignment     = -4108;           % xlCenter
    ur.Interior.Color = 16777215;  % RGB(255,255,255)


    % --- headers ---
    ws.Range('A1:A2').Merge;
    ws.Range('B1:B2').Merge;
    ws.Range('C1:D1').Merge;
    ws.Range('E1:E2').Merge;

    ws.Range('A1:A2').Value = '';

    ws.Range('B1').Value                 = sprintf('Averaging%cMethod',10);
    ws.Range('B1').Font.Name             = 'Palatino Linotype';
    ws.Range('B1').Font.Size             = 12;
    ws.Range('B1').Font.Bold             = true;
    ws.Range('B1').HorizontalAlignment   = -4108;           % xlCenter
    ws.Range('B1').VerticalAlignment     = -4108;
    ws.Range('B1').WrapText              = true;

    % Write, then format the merged C1:D1 via MergeArea
    r                       = ws.Range('C1:D1');                % top-left of the merged block
    r.Merge;
    r.Value                 = SUM.metric;      % e.g. 'nRMSE'
    r.Font.Name             = 'Palatino Linotype';
    r.Font.Size             = 12;
    r.Font.Bold             = true;
    r.HorizontalAlignment   = -4108;           % xlCenter
    r.VerticalAlignment     = -4108;           % xlCenter

    ws.Range('C2').Value                 = 'Mean [95% CI]';
    ws.Range('C2').Font.Name             = 'Palatino Linotype';
    ws.Range('C2').Font.Size             = 12;
    ws.Range('C2').Font.Bold             = true;
    ws.Range('C2').HorizontalAlignment   = -4108;           % xlCenter
    ws.Range('C2').VerticalAlignment     = -4108;


    ws.Range('D2').Value     = 'SD';
    ws.Range('D2').Font.Name             = 'Palatino Linotype';
    ws.Range('D2').Font.Size             = 12;
    ws.Range('D2').Font.Bold             = true;
    ws.Range('D2').HorizontalAlignment   = -4108;           % xlCenter
    ws.Range('D2').VerticalAlignment     = -4108;

    ws.Range('E1').Value     = ['Pairwise improvement in nRMSE ' ...
        'compared with CA [95% CI]'];
    ws.Range('E1').Font.Name             = 'Palatino Linotype';
    ws.Range('E1').Font.Size             = 12;
    ws.Range('E1').Font.Bold             = true;
    ws.Range('E1').HorizontalAlignment   = -4108;           % xlCenter
    ws.Range('E1').VerticalAlignment     = -4108;
    ws.Range('E1:E2').WrapText = true;

    % Merge SimLevel blocks in column A and make bold
    nlevels  = numel(level_names);
    nmethods = numel(avg_names);
    startRow = 3;
    for L = 1:nlevels
        r1 = startRow + (L-1)*nmethods;
        r2  = r1 + nmethods - 1;
        rngA = sprintf('A%d:A%d', r1, r2);
        ws.Range(rngA).Merge;
        ws.Range(sprintf('A%d', r1)).Font.Bold = true;
        ws.Range(sprintf('A%d', r1)).HorizontalAlignment   = -4108;           % xlCenter
        ws.Range(sprintf('A%d', r1)).VerticalAlignment     = -4108;


    end

    % Convert D3:D14 from text to numeric and format to 3 d.p.
    for r = 3:14
        c = ws.Range(sprintf('D%d', r));
        v = c.Value;
        if ischar(v) || isstring(v)
            % strip anything non-numeric just in case (e.g. spaces)
            v = regexprep(char(v), '[^0-9eE\.\+\-]', '');
            num = str2double(v);
            if ~isnan(num)
                c.Value2 = num;                 % write as numeric
            end
        end
    end
    ws.Range('D3:D14').NumberFormat = '0.000';

    % Column widths
    ws.Columns.Item('A:A').ColumnWidth = 13;
    ws.Columns.Item('B:B').ColumnWidth = 13;
    ws.Columns.Item('C:C').ColumnWidth = 19;
    ws.Columns.Item('D:D').ColumnWidth = 13;
    ws.Columns.Item('E:E').ColumnWidth = 34;

    % Row height for all used rows
    lastRow = ws.UsedRange.Rows.Count + ws.UsedRange.Row - 1;
    ws.Rows.Item(sprintf('1:%d', lastRow)).RowHeight = 19;

    % Thick horizontal rules
    % xlEdgeBottom = 9, xlContinuous = 1, xlThick = 4
    thickRows = {'A2:E2','A6:E6','A10:E10','A14:E14'};
    for k = 1:numel(thickRows)
        r = ws.Range(thickRows{k});
        b = r.Borders.Item(9);
        b.LineStyle = 1;
        b.Weight    = 4;
    end

    % Make starred cells in column E red and bold
    for r = 3:lastRow
        cellE = ws.Range(sprintf('E%d', r));
        val = cellE.Value;
        if ischar(val) || isstring(val)
            if contains(string(val), "*")
                cellE.Font.Bold = true;
                cellE.Font.Color = 255;  % RGB(255,0,0)
            end
        end
    end


    % Save and close
    wb.Save;
catch ME
    wb.Close(false);
    ex.Quit; delete(ex);
    rethrow(ME);
end

wb.Close(false);
ex.Quit;
delete(ex);
