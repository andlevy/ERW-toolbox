function [ERA,S] = spm_eeg_era_template(S)
% Creation of temporal warping templates
% FORMAT [ERA,S] = spm_eeg_era_template(S)
%
% S - input structure
%   fields S:
%       S.ERA                         - ERA structure of path name, can be cell 
%                                       array for group template construction
%       S.excludebad                  - exclude bad trial from template
%                                       construction [default: true]
%       S.grouptemplate               - create group average template
%                                       [default: true if supplied with cell array]
%       S.review                      - graphical review of template [default: true]
%       S.averagemethod               - method for average over trials [default: 'median']
%       S.alignevents                 - events between templates to align
%       S.templatesteps               - number of intermediate warping
%                                       templates.
%       S.templatedef.conditionlabel  - condition label for template, can be more 
%                                       than one condition but must contain the same event types
%__________________________________________________________________________
% Andrew Levy

%-Startup
%--------------------------------------------------------------------------
% spm('sFnBanner', mfilename);
% spm_figure('FindWin','Interactive'); spm_clf('Interactive');

try S.excludebad;       catch, S.excludebad     = true;     end
try S.grouptemplate;    catch, S.grouptemplate  = true;     end
try S.review;           catch, S.review         = true;     end
% try S.averagemethod;    catch, S.averagemethod  = 'median'; end
try S.alignevents;      catch, S.alignevents     = {};      end
try S.templatesteps;    catch, S.templatesteps  = 8;        end
try S.save;             catch, S.save           = true;     end

% Load ERA structure
%--------------------------------------------------------------------------
ERA = spm_era_load(S.ERA);
if ~iscell(ERA); ERA = {ERA}; end
nE  = numel(ERA);

% Check design has been created
%--------------------------------------------------------------------------
for e = 1:nE
    if ~isfield(ERA{e},'U')
        error('%s - design has not been specified',fullfile(ERA{e}.swd,ERA{e}.name))
    end
end

% Define template based on subject onset and duration averages
%--------------------------------------------------------------------------
for e = 1:nE

    U     = ERA{e}.U;
    UT    = [];
    Ucond = extractfield(U,'condition');

    for t = 1:numel(S.templatedef)

        u = find(ismember(Ucond,S.templatedef(t).conditionlabel));
        if ~isequal(U(u).events) || ~isequal(U(u).dist)
            error('trial specification is not identical')
        end

        nevents = numel(U(u(1)).events);

        if strcmpi(S.templatedef(t).onset,'mean')
            UT(t).onset     = mean(reshape([U(u).onset],nevents,[])',1);
        elseif strcmpi(S.templatedef(t).onset,'median')
            UT(t).onset     = median(reshape([U(u).onset],nevents,[])',1);
        elseif isnumeric(S.templatedef(t).onset)
            UT(t).onset     = S .templatedef(t).onset;
        end

        if strcmpi(S.templatedef(t).duration,'mean')
            UT(t).duration     = mean(reshape([U(u).duration],nevents,[])',1);
        elseif strcmpi(S.templatedef(t).duration,'median')
            UT(t).duration     = median(reshape([U(u).duration],nevents,[])',1);
        elseif isnumeric(S.templatedef(t).duration)
            UT(t).duration     = S.templatedef(t).duration;
        end

        UT(t).condition = S.templatedef(t).conditionlabel;
        UT(t).events    = U(u(1)).events;
        UT(t).dist      = U(u(1)).dist;
        UT(t).itrial    = u;
        UT(t).type      = 'subject';

    end

    if ~isempty(S.alignevents)

        if ~iscell(S.alignevents); S.alignevents = {S.alignevents}; end
        events = extractfield(UT,'events');

        for a = 1:numel(S.alignevents)

            % Find events
            %-------------------------------------------------
            ons = 0; dur = 0; idx = [];
            for b = 1:numel(events)
                ab = find(strcmp(events{b},S.alignevents{a}));
                if any(ab)
                   ons = ons + UT(b).onset(ab);
                   dur = dur + UT(b).duration(ab);
                   idx = [idx b];
                end
            end

            % Average onsets and duration and replace
            %-----------------------------------------------
            ons = ons/length(idx);
            dur = dur/length(idx);

            for b = idx
                ab                 = find(strcmp(events{b},S.alignevents{a}));
                UT(b).onset(ab)    = ons;
                UT(b).duration(ab) = dur;
            end
        end
    end

    ERA{e}.UT = UT;

end

% Update template spefication if group average requested
%--------------------------------------------------------------------------
if S.grouptemplate

    for t  = 1:numel(S.templatedef)

        % Calculate group average
        %------------------------
        ons = 0;
        dur = 0;
        for e = 1:nE
            ons = ons + ERA{e}.UT(t).onset;
            dur = dur + ERA{e}.UT(t).duration;
        end
        ons = ons./nE;
        dur = dur./nE;

        for e = 1:nE

            ERA{e}.UT(t).onset    = ons;
            ERA{e}.UT(t).duration = dur;
            ERA{e}.UT(t).type     = 'group';

        end
    end
end

% Create transitioned templates
%--------------------------------------------------------------------------
nStep = S.templatesteps;
for e = 1:nE

    % Initialise template
    %---------------------------------
    XT  = zeros(size(ERA{e}.X));
    pst = ERA{e}.pst;
    UT  = ERA{e}.UT;
    U   = ERA{e}.U;

    for i = 1:nStep



        for j = 1:numel(UT)
            for k = UT(j).itrial

                ons = U(k).onset    + (i*(UT(j).onset    - U(k).onset)/nStep);
                dur = U(k).duration + (i*(UT(j).duration - U(k).duration)/nStep);

                % dur = dur_step;

                x = 0;
                for l = 1:length(ons)
                    x = x + spm_era_u(UT(j).dist{l},pst,ons(l),dur(l));
                end

                XT(:,k,i) = x(:);

            end
        end


    end

    ERA{e}.XT = single(XT);

end


if S.save
    for e = 1:nE
        % Save ERA
        %------------------------------------------------------------------
        saveERA(ERA{e})
    end
end


%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG epoching: done'); spm('Pointer','Arrow');

if S.review

    spm_figure('GetWin','ERA Template-Design');
    clf
    UT = ERA{1}.UT;
    nU = numel(UT);
    for t = 1:nU
        subplot(nU,2,1 + (t-1)*2)
        hold on
        plot(pst,ERA{1}.X(:,UT(t).itrial),'Color',[0 0 0 .2])
        plot(pst,ERA{1}.XT(:,UT(t).itrial,end),'LineWidth',2)
        xlabel('Peristimulus')
        axis square
        title(UT(t).condition)

        subplot(nU,2,2+ (t-1)*2)
        X1     = spm_vec(ERA{1}.X(:,UT(t).itrial));
        X2     = spm_vec(ERA{1}.XT(:,UT(t).itrial,end));
        X0     = round(X1,1) ==0 & round(X2,1) ==0;
        X1(X0) = [];
        X2(X0) = [];
        hist3([X1 X2],[5 5], 'CdataMode', 'auto');
        view(2); % View as 2D
        colormap('gray')
        axis square tight
    end

end

if nE == 1; ERA = ERA{1}; end

function saveERA(ERA)
save(fullfile(ERA.swd,ERA.name),'ERA')

