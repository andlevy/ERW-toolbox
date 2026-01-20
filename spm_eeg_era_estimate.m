function [ERA,S] = spm_eeg_era_estimate(S)
% Estimation of warping trjacetories between design and template
% FORMAT [ERA,S] = spm_eeg_era_estimate(S)
%
% S - input structure
%   fields S:
%       S.ERA     - ERA structure of path name
%       S.method  - string of method to esimtate trajectory [default: fltw]


%-Startup
%--------------------------------------------------------------------------
% spm('sFnBanner', mfilename);
% spm_figure('FindWin','Interactive'); spm_clf('Interactive');

try  S.save;   catch,  S.save   = false;  end
try  S.review; catch,  S.review = true;   end
try  S.method; catch,  S.method = 'linear'; end

% Load ERA structure
%--------------------------------------------------------------------------
ERA = spm_era_load(S.ERA);

X               = ERA.X;
F               = ERA.F;
XT              = double(ERA.XT);
T                = size(XT,1);
Y               = XT(:,:,end);
Tau             = []; A = [];

if S.review
    fig             = spm_figure;
    clf
end

for i = 1:size(XT,3)

    switch S.method
        case 'linear'
            [y, X, tau,a] = spm_era_estimate_fltw(X,XT(:,:,i),F);
        case 'hermite'
            [y, X, tau,a] = spm_era_estimate_hermite(X,XT(:,:,i),F);
        case 'parabolic'
            [y, X, tau,a] = spm_era_estimate_parabolic(X,XT(:,:,i),F);
        case 'sinc'
            [y, X, tau,a] = spm_era_estimate_ttw(X,XT(:,:,i),F);
        otherwise
            error('method not recognised')
    end

    A(:,:,i)    = a;
    Tau(:,:,i)  = tau';

    % Review
    %------------------------------------------------
    if S.review
        D = []; c =[];
        % for j = 1:size(X,2)
        %     [cost,d, p] = spm_era_dtw(X(:,j),Y(:,j),200);
        %     c(j)        = cost;
        %     D(:,:,j)    = d;
        % end

        % [~, j] = max(c);
        figure(fig)
        clf
        subplot(2,1,1)
        hold on
        % imagesc(D(:,:,j))
        plot((tau-[1:T])','LineWidth',2,'Color',[0 1 0])
        axis square tight
        colormap('hot')

        subplot(2,1,2)
        hold on
        plot(X,'Color',[0 0 0 .2])
        plot(XT(:,:,end),'LineWidth',1.2,'Color',[1 0 0 .7])
        axis square tight

    
    end

end

% tau = spm_era_fx_tau(mean(A,3),F);
% Store warps
%------------------------------------------------------------------
ERA.Tau    = Tau;
ERA.A      = A;
ERA.MSE    = mean(abs(X-XT(:,:,end)),1);
ERA.method = S.method;


if S.save
    % Save ERA
    %----------------------------------------------------------------------
    save(fullfile(ERA.swd,ERA.name),'ERA')
end

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG epoching: done'); spm('Pointer','Arrow');
