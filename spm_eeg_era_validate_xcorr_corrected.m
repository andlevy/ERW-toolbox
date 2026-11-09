function [D,ERA,S] = spm_eeg_era_validate_xcorr_corrected(S)
% Inter-channel lag preservation: ERW-aligned vs reference averages
%
% Corrected version of spm_eeg_era_validate_xcorr. Changes:
%   1. Channel-pair mask now retains strongly correlated pairs
%      (|r| >= threshold in both conditions), excluding the diagonal.
%      The original mask retained weakly correlated (|r| < threshold),
%      near-zero-lag pairs, i.e. the pairs it was meant to prune.
%   2. The lag error (LE) and sign-agreement proportion (LC) are computed
%      inside the threshold loop, over the masked pairs only, so that each
%      threshold uses its own mask.
%   3. If fewer than two pairs survive the mask, LE is returned as NaN.


%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename);
spm_figure('FindWin','Interactive'); spm_clf('Interactive');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    S.D = D;
end
D = spm_eeg_load(D);

% Get ERA structure
%--------------------------------------------------------------------------
ERA = spm_era_load(S.ERA);


% Calculate channel cross correlation
%--------------------------------------------------------------------------
V   = struct();
Ic  = D.indchantype('EEG','GOOD');
LN  = nan(length(D.indchantype('EEG')),length(D.indchantype('EEG')),2);
XCOV.LL          = {};
XCOV.LX          = {};

for is = 1:numel(S.conds)

    % Get data (averaged ERW-aligned and reference conditions)
    %----------------------------------------------------------------------
    Ie          = D.indtrial({S.conds(is).warp, S.conds(is).ref});
    It          = D.indsample(S.timewin(1)):D.indsample(S.timewin(2));
    Y           = D(Ic,It,Ie);

    [Nc,Nt,Ne]  = size(Y);
    maxlag      = S.xcorr_maxlag /(1/D.fsample);
    lags        = (-maxlag:maxlag) *(1/D.fsample);
    Ycov        = zeros(Nc,Nc,length(lags),2);

    for i = 1:Nc
        for j = 1:Nc
            for k = 1:2
                Ycov(i,j,:,k) = xcov(squeeze(Y(i,:,k)), squeeze(Y(j,:,k)), maxlag, 'coeff');
            end
        end
    end

    % Identify max xcorr and lag
    %----------------------------------------------------------------------
    [M,I]            = max(abs(Ycov),[],3);
    M                = squeeze(M);          % [Nc x Nc x 2]
    I                = squeeze(I);
    L                = lags(I);             % [Nc x Nc x 2]

    xc = S.xcorr_thresh;
    LE = nan(length(xc),1);
    LC = nan(length(xc),1);
    for ix = 1:length(xc)

        % retain strongly correlated pairs in both conditions (off-diagonal)
        %------------------------------------------------------------------
        mask               = all(abs(M) >= xc(ix), 3) & ~eye(Nc);

        % optional: also exclude near-zero (volume-conducted) lags in the
        % reference condition
        % mask             = mask & abs(L(:,:,2)) >= 1/30;

        LX                 = L;
        LX(repmat(~mask,[1 1 2])) = 0;

        % lag-lead directions (agreement with reference), masked pairs only
        %------------------------------------------------------------------
        LL                 = sign(LX);
        LD                 = LL(:,:,1) == LL(:,:,2);

        % lag-lead magnitudes (differences to reference) and reference lags
        %------------------------------------------------------------------
        LM                 = LX(:,:,2) - LX(:,:,1);
        LR                 = LX(:,:,2);

        % normalised RMS lag error over masked pairs
        %------------------------------------------------------------------
        if nnz(mask) > 1
            LE(ix)         = sqrt(mean(LM(mask).^2)) ./ std(LR(mask));
            LC(ix)         = mean(LD(mask));
        end

        % store
        %------------------------------------------------------------------
        XCOV.LL{ix}          = LN;
        XCOV.LX{ix}          = LN;
        XCOV.LL{ix}(Ic,Ic,:) = LL;
        XCOV.LX{ix}(Ic,Ic,:) = LX;
        XCOV.mask{ix}        = mask;
        XCOV.npairs(ix)      = nnz(mask)/2;   % unique pairs (matrix is symmetric)

    end

    XCOV.L(Ic,Ic,:)  = L;

    % Package
    %----------------------------------------------------------------------

    % High level variables
    %----------------------
    V(is).type                = 'xcov';
    V(is).cond_warp           = S.conds(is).warp;
    V(is).cond_ref            = S.conds(is).ref;
    V(is).metric_type         = 'sRMSE';
    V(is).metric              = LE;

    % Low level variables
    %--------------------
    V(is).aux = struct( ...
        'proportion', LC,...
        'lags',lags,...
        'corr_thresh',xc,...
        'XCOV',XCOV );

end

if isfield(ERA,'V')
    ERA.V   = spm_cat_struct(ERA.V,V);
else
    ERA.V    = V;
end

% Save ERA
%--------------------------------------------------------------------------
save(fullfile(ERA.swd,ERA.name),'ERA')

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG cross correlation error: done'); spm('Pointer','Arrow');