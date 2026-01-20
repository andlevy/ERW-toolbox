function [ERA,S] = spm_eeg_era_design(S)
% Defintion of design based on events within a trial
% FORMAT [ERA,S] = spm_eeg_era_design(S)
%
% S - input structure
%   fields S:
%       S.D           - MEEG object or filename of M/EEG mat-file with 
%                       epoch data
%       S.excludebad  - exclude bad trials from the design matrix   [default: true]
%       S.review      - graphic review of design matrix [default: true]
%       S.designdef   - structure array for design definition with fields
%           S.designdef.conditionlabel      - string label for condition
%           S.designdef.onset(i).eventtype  - string for event within trial
%           S.designdef.onset(i).eventvalue - string or numeric
%           S.designdef.onset(i).eventdist  - distribution type (see spm_era_u)
%           S.designdef.onset(i).eventdur   - distribution dispersion (see spm_era_u)
%               First onset is treated as the trial locked event i.e onset t=0
%
% Output:
% ERA                     - Event-Related Alignment structure (ERA)
%   fields:
%       ERA.D      - path to spmeeg object with data
%       ERA.swd    - save working directory
%       ERA.name   - save name for ERA structure 
%       ERA.U      - design specification
%       ERA.X      - design matrix
%__________________________________________________________________________
% Andrew Levy

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
    if ~sts, trl = []; conditionlabels = {}; S = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);
if ~isequal(D.type, 'single')
    error('Design definition requires epoched dataset as input.');
end

%-Get input parameters
%--------------------------------------------------------------------------
try S.excludebad; catch S.excludebad = true; end
try S.review;     catch S.review     = true; end

event    = events(D,1);
if isempty(event)
    error('No event information was found in the input.');
end

conds = extractfield(S.designdef,'conditionlabel');
if ~all(ismember(conds,D.condlist))
    m = find(~ismember(conds,D.condlist));
    error('Conditionlabel(%s) definition not found in spmeeg \n',conds{m})
end


% Create design structure
%--------------------------------------------------------------------------
U = []; u =1;

for i = 1:D.ntrials

    % skip bad trials
    %----------------------------------------------------------------------
    if ismember(i,D.badtrials)
       continue
    end

    % skip trials not specified
    %----------------------------------------------------------------------
    if ~ismember(D.conditions(i),conds)
        continue
    end

    c        = find(strcmp(conds,D.conditions(i)));
    onset    = S.designdef(c).onset;
    Otype    = extractfield(onset,'eventtype');
    Ovalue   = extractfield(onset,'eventvalue');

    event    = events(D,i); 
    event    = event{1};
    Etype    = extractfield(event,'type');
    Evalue   = extractfield(event,'value');
    Etime    = extractfield(event,'time');

    trlons   = D.trialonset(i) - D.time(1);

    ons  = []; dur = []; dist = {};
    for o = 1:numel(onset)

        e = find(strcmp(Etype,Otype{o}) ...
            & strcmp(string(Evalue),string(Ovalue(o))));

        if length(e)>1

            [~,m]  = min(abs(Etime(e) - trlons));
            e      = e(m);

            warning(['Trial %04d, contains more that one event that matches ' ...
                'onset defintion - the event closest to trial onset has ' ...
                'been selected'],i)
        end

        ons  = [ons Etime(e)];
        dur  = [dur onset(o).eventdur];
        dist = [dist onset(o).eventdist];

    end

    U(u).condition = D.conditions{i};
    U(u).onset     = ons-ons(1);
    U(u).events    = Otype;
    U(u).duration  = dur;
    U(u).dist      = dist;
    U(u).itrial    = i;

    u = u +1;
end

ev      = unique([U.events]);
onsets  = cell(1,numel(ev));

for u = 1:numel(U)
    ue = find(ismember(U(u).events,ev));
    for  i  = 1:length(ue)
        onsets{ue(i)}   = [onsets{ue(i)} U(u).onset(i) ];
    end
end

for o = 1:numel(onsets)

    om(o) = mean(onsets{o});
    ou(o) = om(o)+2*std(onsets{o});
    ol(o) = om(o)-2*std(onsets{o});  

end

bad = [];
for u = 1:numel(U)

    ue  = find(ismember(U(u).events,ev));
    ons = U(u).onset;

    if any(ons > ou(ue)) || any(ons < ol(ue))
        bad = [bad u];
    end
end

itrial = extractfield(U,'itrial');
badt   = [badtrials(D) itrial(bad)];
D      = badtrials(D,badt,1);
D      = save(D);
U(bad) = [];

% Create design matrix
%--------------------------------------------------------------------------
X   = [];
pst = D.time;
for i = 1:numel(U)

    % create event-related basis set
    %----------------------------
    x = 0;
    for j = 1:numel(U(i).events)
        x = x + spm_era_u(U(i).dist{j},pst,U(i).onset(j),U(i).duration(j));
    end
    X = [X x(:)];

end

% Add basis set
%--------------------------------------------------------------------------
if isfield(S,'basisset')

    switch S.basisset.type
        case 'sin'
            ERA.F = spm_era_dstmtx(size(X,1),S.basisset.order,S.basisset.tukey);
    end
end

if S.review
    spm_figure('GetWin','ERA Trial-Design');
    subplot(2,1,1)
    plot(pst,X)
    xlabel('Peristimulus')
    axis square
    subplot(2,1,2)
    imagesc([],pst,X)
    ylabel('Peristimulus')
    xlabel('Trial')
    axis square
end

ERA.D    = D.fullfile;
ERA.swd  = D.path; 
ERA.name = 'ERA.mat';
ERA.U    = U;
ERA.X    = X;
ERA.pst  = pst;

% Save ERA
%--------------------------------------------------------------------------
save(fullfile(ERA.swd,ERA.name),'ERA')

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG epoching: done'); spm('Pointer','Arrow');

