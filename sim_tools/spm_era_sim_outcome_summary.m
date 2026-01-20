function [SUM,ESIM] = spm_era_sim_outcome_summary(SIM,export_table)

if nargin < 2
    export_table = [];
end


av_labels  = {'Mean','wMean'};%,'Robust','wRobust','Median','wMedian'};
av_fields  = {'mn','wm'};%,'ra','rb','md','wd'};
nfield       = numel(av_fields);
met_req    = {'nrmse'};
met_labels = {'nRMSE'};

nsim       = size(SIM,1);
nrep       = size(SIM,2);

nwgt       = size(SIM{1,1}.Aw.wm,2);
mnam       = fieldnames(SIM{1,1}(1).E);
im         = arrayfun(@(x) find(strcmp(mnam, met_req{x})), 1:numel(met_req));
nmet       = length(im);

% kappa  = extractfield(SIM{1,1},'Kapppa');
% beta   = extractfield(SIM{1,1},'Beta');
wgt_label = SIM{1,1}.Aw.wlabel;%arrayfun(@(k,b) sprintf('\\kappa: %.2f x \\beta: %.2f', k, b), kappa, beta, 'UniformOutput', false);
label     = ['mn';wgt_label(:)];
% wgt_label  = 1:50;
% Collate outcome measures
%--------------------------------------------------------------------------
E = {}; % E{outcome measure}{sim variant}[nav x nsim x nkappa]
for i = 1:nsim
    for j = 1:nrep
        for k = 1:nmet
            for m = 1:nfield

                e            = SIM{i,j}.E.(mnam{im(k)}).(av_fields{m});
                E{i}{m}(j,:) = e;

            end
        end
    end
end

% Compute and tabulate summary stats
%--------------------------------------------------------------------------
nav     = numel(E{1});  % number of alignment variants
row_idx = 1;
SUM     = table;
PUB     = table;

for i = 1:nsim


    data     = spm_cat(E{i});

    Imn      = find(strcmpi(label,'mn'));
    g1       = squeeze(data(:,Imn));
    Iav      = 1:size(data,2);
    Iav(Imn) = [];
    Ikp      = 1:nwgt;

    P=[];G=[];D=[];C=[];DM=[];RM=[];
    g=1;

    for ikp = 1:length(Ikp)

        g2          = squeeze(data(:,Iav(ikp)));
        D(:,ikp)    = g1-g2;
        R(:,ikp)    = ((g1-g2) ./ g1) *100;
        G(:,ikp)    = ones(size(D(:,ikp)))*g;
        P(ikp)      = sum(g1 > g2)./length(g1);
        DM(ikp)     = mean(D(:,ikp));
        RM(ikp)     = mean(R(:,ikp));
        n               = length(g1);
        se              = std(D(:,ikp)) / sqrt(n);            % standard error
        tval            = tinv(0.975, n - 1);                     % 97.5th percentile t-value for 95% CI
        CI_low          = DM(ikp) - tval * se;                 % lower bound
        CI_high         = DM(ikp) + tval * se;                 % upper bound
        C(:,ikp)        = [CI_low, CI_high];                    % store CI as 3rd dimension

        g = g +1;

    end

    figure
    subplot(3,1,1)
    bar(P)
    title('Proportion reduction in nRMSE')
    subplot(3,1,2)
    bar(DM)
    title('Difference nRMSE')
    subplot(3,1,3)
    bar(P.*DM)
    xtickangle(45)
    % set(gca, 'TickLabelInterpreter', 'latex')
    title('Proprotion weighted Difference nRMSE')


    % CI_low     = squeeze(C(1,:,:));  % lower bounds
    % CI_high    = squeeze(C(2,:,:));  % upper bounds
    % ngroups    = size(DM,1);
    % nbars      = size(DM,2);
    % gw         = min(0.8, nbars/(nbars + 1.5));
    %
    % for ib = 1:nbars
    %     x = (1:ngroups) - gw/2 + (2*ib-1) * gw / (2*nbars);
    %     errorbar(x, DM(:,ib), DM(:,ib) - CI_low(:,ib), CI_high(:,ib) - DM(:,ib), ...
    %         'k', 'linestyle', 'none', 'linewidth', 1);
    % end


    R = squeeze(R);
    G = squeeze(G);

    % figure
    % hold on
    % yline(0, 'k--', 'LineWidth', 1); % horizontal zero line: black dashed
    % h = daviolinplot(D(:), ...
    %     'groups',G(:), ...
    %     'violin','full', ...
    %     'violinwidth',.4, ...
    %     'box',3,...
    %     'boxwidth', 5, ...
    %     'boxalpha', 1,...
    %     'outliers',0, ...
    %     'scatter',0);%
    % ...
    % 'colors',cols);%,...
    % Add significance markers
    % unique_groups = unique(G(:));
    % y_max = max(D(:))*1.5;
    % for i = 1:length(unique_groups)
    %     this_group = unique_groups(i);
    %     if H(this_group)
    %         text(this_group, y_max + 0.02 * range(ylim), '*', ...
    %             'HorizontalAlignment', 'center', ...
    %             'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
    %     end
    % end

    % for v = 1:nav
    %     for l = 1:nwgt
    %
    %         d           = squeeze(data(:,l));
    %         q           = quantile(d,[0.25 0.75]);
    %
    %         SUM.SimLevel{row_idx,1}    = i;
    %         SUM.Kappa{row_idx,1}       = wgt_label(l);
    %         SUM.AvVariant{row_idx,1}   = av_labels{v};
    %         SUM.Metric{row_idx,1}      = met_labels{k};
    %         SUM.Mean(row_idx,1)        = mean(d);
    %         SUM.SD(row_idx,1)          = std(d);
    %         SUM.Median(row_idx,1)      = median(d);
    %         SUM.Q1(row_idx,1)          = q(1);
    %         SUM.Q3(row_idx,1)          = q(2);
    %         SUM.IQR(row_idx,1)         = q(2) - q(1);
    %
    %         PUB.SimLevel{row_idx,1}    = i;
    %         PUB.Kappa{row_idx,1}       = wgt_label(l);
    %         PUB.AvVariant{row_idx,1}   = av_labels{v};
    %         PUB.Mean{row_idx,1}        = sprintf('%.03f (%.03f)',SUM.Mean(row_idx,1),SUM.SD(row_idx,1));
    %         PUB.Median{row_idx,1}      = sprintf('%.03f (%.03f)',SUM.Median(row_idx,1),SUM.IQR(row_idx,1));
    %
    %         if strcmpi(SUM.AvVariant{row_idx,1},'mean')
    %
    %             Base     = SUM.Mean(row_idx,1);
    %             PUB.Change{row_idx,1} = '--';
    %
    %         else
    %
    %             Change   = ((SUM.Mean(row_idx,1) - Base) ./ Base) *100;
    %
    %             if sign(Change) < 0
    %                 dir = 'lower';
    %             else
    %                 dir = 'higher';
    %             end
    %
    %             PUB.Change{row_idx,1} = sprintf('%.02f%s',abs(Change),['% ' dir]);
    %
    %         end
    %         row_idx = row_idx + 1;
    %     end
    % end
end


disp(PUB)

% vars_to_round = {'Mean','SD','Median','Q1','Q3','IQR',};
% for v = vars_to_round
%     SUM.(v{1}) = round(SUM.(v{1}), 6);  % or whatever precision you like
% end
% 
% if ~isempty(export_table)
%     writetable(SUM, [ export_table '_summary.xlsx'], 'Sheet', 1)
%     writetable(PUB, [ export_table '_results.xlsx'], 'Sheet', 1)
% end
% 
% % Get index for examplars
% %--------------------------------------------------------------------------
% 
% % index for nRMSE
% %-----------------
% idxStat = find(ismember(SUM.Metric,'nRMSE'));
% idxSim  = cell2mat(SUM.SimLevel(idxStat));
% 
% % Identify which sim level performed worst
% %-----------------------------------------
% SimMed = []; SimQ3 = [];
% for i = unique(idxSim)'
% 
%     SimMed(:,i) = (SUM.Median(idxStat(idxSim==i),1));
%     SimQ3(:,i)  = (SUM.Q3(idxStat(idxSim==i),1));
% 
% end

ESIM = struct();
% [~,iSim] = max(max(SimMed));
%
% % Find examplar SIM
% %-----------------------------------------
% [~,iMed] = min(sum(abs( E{1}{iSim} - SimMed(:,iSim) ),1));
% [~,iQ3]  = min(sum(abs( E{1}{iSim} - SimQ3(:,iSim)  ),1));
%
%
% ESIM.Med.SIM      = SIM{iSim,iMed};
% ESIM.Med.descript = sprintf('SimLevel:%i',iSim);
% ESIM.Q3.SIM       = SIM{iSim,iQ3};
% ESIM.Q3.descript  = sprintf('SimLevel:%i',iSim);