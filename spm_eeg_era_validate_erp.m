function [D,ERA,S] = spm_eeg_era_validate_erp(S)


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

% Calculate RMSE wrt reference ERP
%--------------------------------------------------------------------------
V   = struct();
Ic  = D.indchantype('EEG','GOOD');

for is = 1:numel(S.conds)

    It          = D.indsample(S.timewin(1)):D.indsample(S.timewin(2));
    Ie          = D.indtrial({S.conds(is).warp, S.conds(is).ref});

    X           = D(Ic,It,Ie(1));
    R           = D(Ic,It,Ie(2));
  
    Y           = cat(3,X,R);
    diff_Y      = diff(Y, 1, 3);  % Difference across trials (should be only 2 trials)

    nRMSE(Ic) = squeeze((sqrt( mean(diff_Y.^2, 2) ) ./ std(Y(:,:,2), 0, 2)));  % Normalised RMSE

    V(is).type              = 'erp';
    V(is).cond_warp         = S.conds(is).warp;
    V(is).cond_ref          = S.conds(is).ref;
    V(is).metric_type       = 'nRMSE';
    V(is).metric            = nRMSE;
    V(is).aux               = struct();

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
