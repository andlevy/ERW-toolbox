function fig = spm_era_sim_results_figure(SUM,fig_type,save_fig,save_dir)


cold    = [
    142  69 163;
    141 182  20;
    242 102  20;
    20 149 170] / 255;
cold = brighten(cold,-.33);
colb1   = brighten(cold,.66);
colb2   = brighten(colb1,.66);
shades  = {cold, colb1, colb2};
% err     = brighten(cold,-.65);
shades = {colb1};

fig = figure('Units','centimeters','Position',[0 10 35 10]);
set(fig, 'Color', 'w', ...
    'DefaultAxesFontName','Roboto', ...
    'DefaultAxesFontSize',14, ...
    'DefaultTextFontName','Roboto', ...
    'DefaultTextFontSize',14, ...
    'DefaultLineLineWidth',1.5, ...
    'DefaultLineMarkerSize',6);

switch fig_type

    case 'violin'
        s   = 1;

        up = max(squeeze(SUM.E(:,1,:,:)),[],'all');
        up = max(up,.7);
        for l = 1:3
            ax = subplot(1,3,s);
            hold on
                yline(0.2956, 'r-','LineWidth',4);
            D  = squeeze(SUM.E(:,1,:,l));   % data (nG x nS)
            G  = squeeze(SUM.G(:,1,:,1));
            h  = daviolinplot(D(:), ...
                'groups',G(:), ...
                'violin','full', ...
                'violinwidth',.4, ...
                'box',3,...
                'boxwidth', 5, ...
                'boxalpha', 1,...
                'outliers',0, ...
                'scatter',0,...
                'colors',shades{1});%,...
            axis  square
        
            xlim([.5 4.5])
            ylim([0.1 up])
            xticks([])
            box off

            if l ==1
                ylabel([SUM.metric ' (w.r.t. ground truth)'])
            end

            if l == 3

                algNames = { ...
                    'Conventional Avg', ...
                    'Weighted Avg (κ: 1)', ...
                    'Weighted Avg (κ: 2)', ...
                    'Weighted Avg (κ: 3)'};
                hAlg = gobjects(4,1);
                for k = 1:4
                    hAlg(k) = bar(NaN, NaN, 'FaceColor', shades{1}(k,:));
                    hold on
                end

                lg          = legend(hAlg, algNames);
                lg.Position =  [0.13,0.79,.15,.19];
                lg.FontSize = 10;
                lg.Box      = 'off';
                hold off

            end

            s = s + 1;
        end

        fig_title = ['nRMSE_Violin_Sim' upper(SUM.sumnum)];

        set(fig,'name',fig_title)

    case 'bar'

        s   = 1;

        for l = 1:3
            subplot(1,3,s)

            D           = SUM.ME(1,:,l)';   % data (nG x nS)
            E           = SUM.SD(1,:,l)';   % standard errors (nG x nS) same orientation
            G           = SUM.G(1,:,1)';
            CI5         = SUM.CI5(1,:,l)';
            CI95        = SUM.CI95(1,:,l)';
            [nG,nS]     = size(D);
            groupWidth  = 0.9;
            barWidth    = groupWidth / nS;
            xBase       = 1:nG;

            hold on
            for j = 1:nS

                x    = xBase - groupWidth/2 + (j-0.5)*barWidth;
                h(j) = bar(x, D(:,j), barWidth, ...
                    'FaceColor','flat', 'CData', shades{j});%, ...
                % 'EdgeColor','none');
                for i = 1:length(x)
                    errorbar(x(i), D(i,j), CI5(i,j), CI95(i,j), ...
                        'color', [.2 .2 .2],...%err(i,:), ...                 % colour of error bar lines
                        'LineStyle','none', ...
                        'LineWidth',5, ...
                        'CapSize',0);           % adjust / set to 0 for no caps
                end

            end


            % ylim([0 0.4])
            xlim([0.5 nG+0.5])
            xticks([])
            box off

            if l ==1
                ylabel(SUM.metric)
            end

            if l == 3

                algNames = { ...
                    'Conventional Avg', ...
                    'Weighted Avg (κ: 0.5)', ...
                    'Weighted Avg (κ: 1)', ...
                    'Weighted Avg (κ: 2)'};
                hAlg = gobjects(4,1);
                for k = 1:4
                    hAlg(k) = bar(NaN, NaN, 'FaceColor', cold(k,:));
                    hold on
                end

                lg          = legend(hAlg, algNames);
                lg.Position =  [.79,.79,.18,.19];
                lg.FontSize = 10;
                lg.Box      = 'off';

                hold off

            end

            s = s + 1;
        end

        fig_title = ['nRMSE_Bar_Sim' upper(SUM.sumnum)];

        set(fig,'name',fig_title)

    case 'time'


        for v = 1:3

            Y      = struct();
            Y.ge   = squeeze(SUM.Yg(:,:,1,1,v));
            Y.g    = squeeze(SUM.GT(:,:,1,1,v));
            Y.wcm  = squeeze(SUM.Yw(:,:,1,1,v));
            Y.wwm  = squeeze(SUM.Yw(:,:,1,2,v));
            % Y.wjm  = squeeze(SUM.Yw(:,:,3,3,v));
            Y.u    = squeeze(SUM.Yu(:,:,1,1,v));
            Y.e    = squeeze(SUM.Ye(:,:,1,1,v));

            G    = Y.ge;
            Wcm  = Y.wcm;
            Wwm  = Y.wwm;
            % Wjm  = Y.wjm;
            U    = Y.u;
            N    = Y.e;

            % residuals
            Ewcm = Wcm - G;
            Ewwm = Wwm - G;
            % Ewjm = Wjm - G;
            Eu   = U - G;

            % preallocate
            RMSE_wcm   = zeros(size(G,1),1);
            % RMSE_wjm   = zeros(size(G,1),1);
            RMSE_wwm   = zeros(size(G,1),1);
            RMSE_u     = zeros(size(G,1),1);
            RMSE_noise = zeros(size(G,1),1);

            for ti = 1:size(G,1)
                % RMSE
                RMSE_noise(ti) = sqrt(mean(N(ti,:).^2));
                if RMSE_noise(ti) == 0
                    RMSE_noise(ti) = 2;
                end
                RMSE_wcm(ti)   = sqrt(mean(Ewcm(ti,:).^2))./ RMSE_noise(ti);
                RMSE_wwm(ti)   = sqrt(mean(Ewwm(ti,:).^2))./ RMSE_noise(ti);
                % RMSE_wjm(ti)   = sqrt(mean(Ewjm(ti,:).^2))./ RMSE_noise(ti);
                RMSE_u(ti)     = sqrt(mean(Eu(ti,:).^2))./ RMSE_noise(ti);

            end
            T = size(G,1);
            Fs = 500;
            t = ((0:T-1)/Fs - .1);  % time in seconds
            % Plot

            ax = subplot(1,3,v);
            hold on
            plot(t, RMSE_u, 'color',[0 0 0], 'LineWidth',4);
            plot(t, RMSE_wcm, 'color',colb1(1,:), 'LineWidth',4);
            % plot(t,RMSE_wjm,'color',colb1(2,:), 'LineWidth',4);
            plot(t,RMSE_wwm, 'color',colb1(2,:), 'LineWidth',4);

            yline(1, '--','LineWidth',4);
            xlabel('PST (s)');
            ylim([0 5])
            xr          = xlim(ax);
            rawStep     = (xr(2)-xr(1))/6;
            prettyStep  = round(rawStep,1);
            xt          = xr(1):prettyStep:xr(2);
            set(ax,'XTick',    xt)
            if v ==1
                ylabel('nRMSE (w.r.t. noise)');
            end
            if v == 3
                lg = legend('Unwarped','Warped (conventional)','Warped (weighted κ: 1)');
                lg.FontSize = 10;
                lg.Position = [.12,.79,.18,.19];
                lg.Box      = 'off';
            end
        end

        fig_title = ['nRMSE_Time_Sim' upper(SUM.sumnum)];
        set(fig,'name',fig_title)

    case 'response'

        Is = 1;
        pst = SUM.pst;
        for r = 1:3

            ax = subplot(1,3,Is);
            hold on

            for Ia = 1
                for Ik = 2
                    plot(pst,squeeze(SUM.Yw(:,SUM.Isim,Ia,Ik,r)),'color',colb1(2,:),'LineWidth',4)
                    axis tight
                end
            end

            plot(pst,squeeze(SUM.Yu(:,SUM.Isim,1,1,r)) ,'--','color',[0 0 0],'LineWidth',4)
            plot(pst,squeeze(SUM.Yg(:,SUM.Isim,1,1,r)),'--','color',[1 0 0 1],'LineWidth',4)
            xlabel('PST (s)');
            if Is == 1
            ylabel('Amplitude (a.u.)')
            end
            xr          = xlim(ax);
            rawStep     = (xr(2)-xr(1))/6;
            prettyStep  = round(rawStep,1);
            xt          = xr(1):prettyStep:xr(2);
            set(ax,'XTick',    xt)
            if Is == 3
                lg = legend('Warped (weighted κ: 1)','Unwarped',' Ground Truth');
                lg.FontSize = 10;
                lg.FontWeight = 'bold';
                lg.Position = [.83,.23,.15,.19];
                lg.Box      = 'off';
            end

            Is = Is +1;

        end

        fig_title = ['RecoveredResponse_Sim' upper(SUM.sumnum)];
        set(fig,'name',fig_title)

end

if save_fig

    set(fig,'Renderer','painters');
    print(fig,fullfile(save_dir,['Result_' fig_title]),'-dsvg','-r600');

end
