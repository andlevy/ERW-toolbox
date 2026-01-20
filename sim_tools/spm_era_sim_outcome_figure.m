function fig = spm_era_sim_outcome_figure(type,SIM,varargin)


switch lower(type)

    case 'simulated responses'

        legendN = varargin{1}; % legend names

        col = [13 152 186; 255 174  66; 199  21 133]./255;

        fig = figure('Units','normalized','Position',[.1 .1 .5 .5]);
        hold on;
        legendH = [];
        for i = 1:size(SIM,1)
            pst              = SIM{i,1}.pst;
            h                = plot(pst,SIM{i,1}.Y,'Color',[col(i,:) .5],'LineWidth',3);
            legendH(end+1)   = h(1);

        end
        h1             = plot(pst,SIM{1,1}.GT1,'k:','LineWidth', 4);
        legendH(end+1) = h1;
        h2             = plot(pst, SIM{1,1}.GT2, 'k', 'LineWidth', 4);
        legendH(end+1) = h2;
        legend(legendH, legendN,'location','best');
        xlabel('Peristimulus time (s)');
        ylabel('Amplitude');
        axis tight square

    case {'evaluation metrics','evaluation metrics -kappa'}


        % met_req = {'r2','nrmse'};
        % met_lab = {'R^2','nRMSE'};
        % met_sup    = {[0 1],[0 inf]};
        % met_lim    = {[.85 1],[0 .5]};

        met_req = {'nrmse'};
        met_lab = {'nRMSE'};
        met_sup = {[0 inf]};


        nsim       = size(SIM,1);
        nrep       = size(SIM,2);
        nwgt       = numel(SIM{1,1});
        mnam       = fieldnames(SIM{1,1}(1).E);
        im         = arrayfun(@(x) find(strcmp(mnam, met_req{x})), 1:numel(met_req));
        nmet       = length(im);

        E = {};
        if contains(type,'-kappa')

            for i = 1:nsim
                for j = 1:nrep
                    for k = 1:nmet

                        E{k}{i}(j,1) = SIM{i,j}(1).E.(mnam{im(k)}).mn;

                        for l = 1:nwgt

                            E{k}{i}(j,l+1) = SIM{i,j}(l).E.(mnam{im(k)}).wm;

                        end
                    end
                end
            end
            xlabels    = {'CA','wCA:1','wCA:2','wCA:3','wCA:4'};
            G          = [];
            G(:,:,1)   = repmat([1:length(xlabels)]',1,nrep)';
            ncol       = 1;
            cold       = [142 69 163;141 182 20]./255;
            colb(1,:)  = brighten(cold(2,:),.3);
            colb(2,:)  = brighten(cold(2,:),.6);
            colb(3,:)  = brighten(cold(2,:),.9);
            cols       = [cold(1,:);cold(2,:);colb];


        else

            for i = 1:nsim
                for j = 1:nrep
                    for k = 1:nmet
                        for l = 1:nwgt

                            E{k,l}{i}(:,j) = [
                                SIM{i,j}(l).E.(mnam{im(k)}).mn,...
                                SIM{i,j}(l).E.(mnam{im(k)}).wm]';

                        end
                    end
                end
            end
            xlabels   = {'CA','wCA'};%,'RA','wRA','MA','wMA'};
            G        = [];
            G(:,:,1) = repmat([1:length(xlabels)]',1,nrep);


            [a,b]   = size(E);
            if a*b ~= max([a b])
                error('Can only support multi measure or multi weight not both')
            else
                ncol   = nwgt * nmet;
            end

            cold       = [142 69 163; 141 182 20; 242 102 20]./255;
            colb       = brighten(cold,.75);
            cols       = [cold(1,:);colb(1,:);cold(2,:);colb(2,:);cold(3,:);colb(3,:)];

        end

        for iE = 1:numel(E)
            for iF = 1:numel(E{iE})
                met_lim(iE,iF) = max(spm_vec(E{iE}{iF}));
            end
        end

        if size(met_lim,1) > 1
            met_lim = ceil(max(met_lim)*10)./10;
        end

        fig      = figure('Units','normalized','Position',[.51,.04,.2,.89]);
        set(fig, 'Color', 'w', ...
            'DefaultAxesFontName','Roboto', ...
            'DefaultAxesFontSize',14, ...
            'DefaultTextFontName','Roboto', ...
            'DefaultTextFontSize',14, ...
            'DefaultLineLineWidth',1.5, ...
            'DefaultLineMarkerSize',6);

        if contains(type,'-kappa')

            for i = 1:nsim
                for j = 1:ncol
                    s = subplot(nsim, ncol, (i - 1) * ncol + j);
                    h = daviolinplot(E{j}{i}(:), ...
                        'groups',G(:), ...
                        'violin','full', ...
                        'violinwidth',.4, ...
                        'box',3,...
                        'boxwidth', 5, ...
                        'boxalpha', 1,...
                        'outliers',0, ...
                        'scatter',0,...
                        'colors',cols);%,...
                    % 'support', met_sup);
                    % ylabel(met_lab{j})
                    box(s,'on')
                    axis  square
                    xlim([.5 5.5])
                    ylim([0 met_lim(i)])

                    if i == nsim
                        xticklabels(xlabels)
                    else
                        xticklabels(xlabels)
                        % xticks([])
                    end
                end
            end


        else

            for i = 1:nsim
                for j = 1:ncol
                    s = subplot(nsim, ncol, (i - 1) * ncol + j);
                    h = daviolinplot(E{j}{i}(:), ...
                        'groups',G(:), ...
                        'violin','full', ...
                        'violinwidth',.4, ...
                        'box',3,...
                        'boxwidth', 5, ...
                        'boxalpha', 1,...
                        'outliers',0, ...
                        'scatter',0,...
                        'colors',cols);%,...
                    % 'support', met_sup);
                    % ylabel(met_lab{j})
                    box(s,'on')
                    axis  square
                    xlim([.5 2.5])
                    ylim([0 met_lim(i)])

                    if i == nsim
                        xticklabels(xlabels)
                    else
                        xticklabels(xlabels)
                        % xticks([])
                    end
                end
            end

        end

    case 'fitted responses'

        % Establihsh if entire SIM space is provided or examplar SIMs
        %---------------------------------------------------------------
        if isstruct(SIM) % Exemplars

            fnames = fieldnames(SIM);

            FIM = {};
            for f = 1:length(fnames)
                FIM{f} = SIM.(fnames{f}).SIM;
            end

        else

            k    = randi(size(SIM,2),1,1);
            if ~isempty(varargin)
                iSIM = varargin{1}; % which SIMs to plot
            end

        end

        plabels  = {'ground truth','mean','w.mean','rob.av','w.rob.av','median','w.median'};
        cold     = [142 69 163; 141 182 20; 242 102 20]./255;
        colb     = brighten(cold,.75);
        cols     = [cold(1,:);colb(1,:);cold(2,:);colb(2,:);cold(3,:);colb(3,:)];

        fig   = figure('Units','normalized','Position',[0.01,.05,.94,.69]);
        set(fig, 'Color', 'w', ...
            'DefaultAxesFontName','Roboto', ...
            'DefaultAxesFontSize',14, ...
            'DefaultTextFontName','Roboto', ...
            'DefaultTextFontSize',14, ...
            'DefaultLineLineWidth',1.5, ...
            'DefaultLineMarkerSize',6);


        nFim  = numel(FIM);

        for i = 1:nFim
            % ──────────────────────────────────────────────────────────
            pst = FIM{i}.pst;
            GT  = FIM{i}.GT2;
            X   = [
                FIM{i}.A.mn, ...
                FIM{i}.A.b,  ...
                FIM{i}.A.ra, ...
                FIM{i}.A.rb, ...
                FIM{i}.A.md, ...
                FIM{i}.A.wd ];

            V                   = max(abs(X - GT), [], 2);
            dt                  = pst(2) - pst(1);
            k                   = fix(0.075 / dt);
            W                   = spm_convmtx(ones(k,1), length(V), 'square');
            W(:, sum(W,1) < k)  = [];
            [~, m]              = max(V' * W);
            w                   = logical(W(:, m));

            % ──────────────────────────────────────────────────────────
            main_ax = subplot(1, nFim, i);
            set(main_ax, 'Units','normalized', ...
                'LineWidth',1, ...
                'Box','off', ...
                'TickDir','out', ...
                'TickLength',[0.02 0.02]);
            axis(main_ax,'tight');

            xlabel(main_ax,'pst (s)');
            ylabel(main_ax,'amplitude (au)')

            hold(main_ax, 'on');

            xlims_main = [ min(pst), max(pst) ];
            ylims_main = [ min(GT)*1.1, max(GT)*1.1 ];

            plot(main_ax, pst, GT,         ':', 'Color',[1 0 0], 'LineWidth', 7);
            plot(main_ax, pst, FIM{i}.A.mn, 'Color', cols(1,:), 'LineWidth', 3);
            plot(main_ax, pst, FIM{i}.A.b,  'Color', cols(2,:), 'LineWidth', 3);
            plot(main_ax, pst, FIM{i}.A.ra, 'Color', cols(3,:), 'LineWidth', 3);
            plot(main_ax, pst, FIM{i}.A.rb, 'Color', cols(4,:), 'LineWidth', 3);
            plot(main_ax, pst, FIM{i}.A.md, 'Color', cols(5,:), 'LineWidth', 3);
            plot(main_ax, pst, FIM{i}.A.wd, 'Color', cols(6,:), 'LineWidth', 3);

            xlabel('PST (s)')
            ylabel('Amplitude (au)')
            xlim(main_ax, xlims_main);
            ylim(main_ax, ylims_main);

            if i == nFim

                l = legend(plabels);
                set(l,'Position',[.91,.73,.08,.19]);
                set(l,'Box','off','FontSize',12);

            end

            % ──────────────────────────────────────────────────────────
            main_pos = get(main_ax, 'Position');
            x_marg   = 0;
            y_marg   = 0.07;
            inset_wf = 0.30;

            inset_pos = [ ...
                main_pos(1) + (1 - inset_wf - x_marg) * main_pos(3), ...  % left
                main_pos(2) + (    y_marg           ) * main_pos(4), ...  % bottom
                .15 .39];

            inset_ax = axes('Units','normalized', 'Position', inset_pos);
            hold(inset_ax, 'on');
            box(inset_ax, 'on');

            xlims_ins = [ min(pst(w)), max(pst(w)) ];
            ylims_ins = [ min(GT(w)) - 2, max(GT(w)) + 2 ];

            plot(inset_ax, pst(w), GT(w),         ':', 'Color',[1 0 0], 'LineWidth', 7);
            plot(inset_ax, pst(w), FIM{i}.A.mn(w),     'Color', cols(1,:), 'LineWidth', 3);
            plot(inset_ax, pst(w), FIM{i}.A.b(w),      'Color', cols(2,:), 'LineWidth', 3);
            plot(inset_ax, pst(w), FIM{i}.A.ra(w),     'Color', cols(3,:), 'LineWidth', 3);
            plot(inset_ax, pst(w), FIM{i}.A.rb(w),     'Color', cols(4,:), 'LineWidth', 3);
            plot(inset_ax, pst(w), FIM{i}.A.md(w),     'Color', cols(5,:), 'LineWidth', 3);
            plot(inset_ax, pst(w), FIM{i}.A.wd(w),     'Color', cols(6,:), 'LineWidth', 3);
            xlim(inset_ax, xlims_ins);
            ylim(inset_ax, ylims_ins);

        end

        % print(fig,'‑dpng','‑r300','myFigure.png');

end