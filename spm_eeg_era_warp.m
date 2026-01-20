function [Da,ERA,S] = spm_eeg_era_warp(S)
% Applying warp trajectories to spmeeg timeseries
% FORMAT [ERA,S] = spm_eeg_era_warp(S)
%
% S - input structure
%   fields S:
%       S.D    - MEEG object or filename of M/EEG mat-file with 
%                epoch data
%       S.ERA  - ERA structure of path name, can be cell array for group 
%                template construction
%__________________________________________________________________________
% Andrew Levy

try S.prefix;           catch, S.prefix   = 'A';   end
try S.channels;         catch, S.channels = 'all'; end

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
if ~isequal(D.type, 'single')
    error('Design definition requires epoched dataset as input.');
end

% Get ERA structure
%--------------------------------------------------------------------------
ERA = spm_era_load(S.ERA);

% Apply warp
%--------------------------------------------------------------------------
switch ERA.method
    case {'parabolic','fltw'}

        It          = extractfield(ERA.U,'itrial');
        if strcmpi(S.channels,'all')
            Ic      = 1:D.nchannels;
        else
            Ic      = D.indchantype(S.channels);
        end
        Y           = D(Ic,:,It);
        Yw          = spm_era_warp_fltw(Y,ERA.Tau);
        Ya          = D(:,:,:);
        Ya(Ic,:,It) = Yw;

end

% Write warped spmeeg object
%--------------------------------------------------------------------------
Da        = clone(D,fullfile(D.path,spm_file(D.fname,'prefix',S.prefix)),size(D));
Da(:,:,:) = Ya;
save(Da)

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG temporal aligning: done'); spm('Pointer','Arrow');


