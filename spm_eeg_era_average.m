function D = spm_eeg_era_average(S)

% Average each channel over trials or trial types using robust, trimmed, or barycentric averaging
% FORMAT D = spm_eeg_era_average(S)
%
% S         - input struct with fields:
%   S.D            - MEEG object or filename of M/EEG mat-file with epoched data
%   S.ERA          - ERA structure returned by spm_era_load (required for barycentric weighting)
%   S.robust       - (optional) structure to enable robust averaging
%       .savew         - save the weights in an additional dataset (default: 0)
%       .bycondition   - compute weights separately for each condition (default: 0)
%       .ks            - offset of the weighting function (default: 3)
%       .removebad     - ignore bad samples during averaging (default: 0)
%   S.barycentre   - (optional) apply barycentric weighting from ERA.Tau (default: 0)
%   S.trim         - (optional) trim mean by a percentile (e.g. 10% trim: S.trim=10) (default: 0)
%   S.prefix       - prefix for the output file (default: 'm')
%
% Output:
%   D         - MEEG object containing the averaged data (also saved to disk)
%__________________________________________________________________________
%
% Updated for full support of robust and barycentric averaging combinations
% Stefan Kiebel and Andrew Levy
% Copyright (C) 2008-2025 Wellcome Centre for Human Neuroimaging

%-Startup
%-----------------------------------------------------------------------
spm('FnBanner', mfilename);
spm('FigName','M/EEG ERA averaging'); spm('Pointer','Watch');

%-Get input parameters
%-----------------------------------------------------------------------
if ~isfield(S, 'prefix'),       S.prefix     = 'm';           end
if ~isfield(S, 'robust'),       S.robust     = 0;             end
if ~isfield(S, 'weight'),   S.weight = 0;             end
if ~isfield(S, 'trim'),         S.trim       = 0;             end
if ~isfield(S, 'median'),       S.median     = 0;             end

%-Configure robust averaging
%-----------------------------------------------------------------------
if isstruct(S.robust)

    if ~isfield(S.robust, 'savew'),        S.robust.savew       = 0;   end
    if ~isfield(S.robust, 'bycondition'),  S.robust.bycondition = 0;   end
    if ~isfield(S.robust, 'ks'),           S.robust.ks          = 3;   end
    if ~isfield(S.robust, 'removebad'),    S.robust.removebad   = 0;   end

    robust      = 1;
    savew       = S.robust.savew;
    bycondition = S.robust.bycondition;
    ks          = S.robust.ks;
    removebad   = S.robust.removebad;

else
    robust = 0;
end

if isstruct(S.weight)
    weight = 1;
    kappa = S.weight.kappa;
    beta  = S.weight.beta;
else
    weight = 0;
end

% barycentre = S.weight;
med        = S.median;

%-Load ERA structure
%-----------------------------------------------------------------------
ERA = spm_era_load(S.ERA);

%-Get MEEG object
%-----------------------------------------------------------------------
D           = spm_eeg_load(S.D);
montage_ind = D.montage('getindex');
D           = D.montage('switch', 0);

%-Check that there is any good data available
%-----------------------------------------------------------------------
if ntrials(D)==0 || isempty(indtrial(D, D.condlist, 'GOOD'))
    warning('No good trials for averaging were found. Nothing to do.');
    return;
end

%-Generate new MEEG object with new files
%-----------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nsamples D.nconditions]);
Dnew = type(Dnew, 'evoked');

if weight

    %-Configure tau weighting
    %----------------------------------------------------------------------
    Tau         = ERA.Tau;
    nsamp       = size(Tau,1);
    b           = sum(abs(sum( Tau-[1:nsamp]', 3)),1);
    % b           = b - mean(b);

    % Gaussian Weighting
    %-------------------
    % kappa   = 1;
    % beta    = 2;
    % b       = b - mean(b);
    % sig     = std(b);  % adapt width to data
    b       = exp(-(abs(b)./(sqrt(2)*kappa.*std(b))).^beta);

    if all(b==0) % if no onset variability
        b = ones(size(b));
    end
    % Normalise
    %-----------------
    b = b./sum(b);
else
    b = ones(1,D.ntrials);

end




% b           = 1 ./ (b + eps);                  % inverse deviation
% lambda      = 0.1 * max(b);                % regularisation scale
% b           = b ./ (b + lambda);           % squash high values
% b           = b / sum(b);

utrial      = extractfield(ERA.U,'itrial');
dcl         = D.condlist;
uc          = extractfield(ERA.U,'condition');
ucl         = unique(uc);
B           = cell(1,D.nconditions);
Bi          = cell(1,D.nconditions);

for i = 1:D.nconditions

    % Parition weights
    %-------------------------------------------------------------------
    if ismember(dcl(i),ucl)

        j       = ismember(uc,dcl(i));

        % Check indexes align
        %---------------------------------------------------------------
        assert(isempty(setdiff(utrial(j),D.indtrial(dcl(i),'GOOD'))), ...
            'Mismatch between trial indices')

        Bi{i}   = D.indtrial(dcl(i),'GOOD');
        B{i}    = b(j)./sum(b(j));

    else

        Bi{i}   = D.indtrial(dcl(i),'GOOD');
        B{i}    = ones(1,length(Bi{i}))/length(Bi{i});

    end
end

ni = zeros(1, D.nconditions);
for i = 1:D.nconditions

    w       = indtrial(D, deblank(dcl{i}), 'GOOD');
    ni(i)   = length(w);

    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end

end

goodtrials  = indtrial(D, dcl, 'GOOD');
chanind     = D.indchantype('All');%D.indchantype({'MEEG', 'LFP'});

if prod(size(D))*8 < spm('memory')

    %-Average everything at once if there is enough memory
    %===================================================================
    spm_progress_bar('Init', D.nconditions, 'Averages done');
    if D.nconditions > 100,
        Ibar = floor(linspace(1, D.nconditions, 100));
    else
        Ibar = 1:D.nconditions;
    end

    if robust && ~bycondition

        W1                  = ones(D.nchannels, D.nsamples, length(goodtrials));
        % B1                  = ones(D.nchannels, D.nsamples, length(goodtrials));
        Y                   = D(chanind, :, goodtrials);

        if removebad
            bad    = badsamples(D, chanind, ':', goodtrials);
            Y(bad) = NaN;
        end

        if weight
            B1                  = spm_vec(B);
            B1                  = repmat(B1,1,length(chanind),D.nsamples);
            B1                  = permute(B1,[2 3 1]);
            [Y, W2]             = spm_robust_barycentric(Y, 3, ks,B1);
        else
            [Y, W2]             = spm_robust_average(Y, 3, ks);
        end

        W1(chanind, :, :)   = W2;
        W                   = zeros(size(D));
        W(:, :, goodtrials) = W1;

        if savew
            Dw(:, :, goodtrials) = W1;
        end

    end

    for i = 1:D.nconditions

        w = indtrial(D, deblank(dcl{i}), 'GOOD');

        if isempty(w)
            continue;
        end

        if ~robust

            if weight
                Y               = D(chanind, :, w);
                B1              = spm_vec(B{i});
                if med
                    Dnew(:, :, i) = medianw( Y, B1, 3 );
                else
                    
                  
                    B1              = repmat(B1,1,length(chanind),D.nsamples);
                    B1              = permute(B1,[2 3 1]);
                    Dnew(:, :, i)   = sum(B1.*Y,3)./sum(B1,3);
                end

            else
                if med
                    Dnew(:, :, i) = mean(D(:, :, w), 3);
                else
                    Dnew(:, :, i) = median(D(:, :, w), 3);
                end

            end

        else

            if bycondition

                W       = ones(D.nchannels, D.nsamples, length(w));
                Y       = D(chanind, :, w);

                if removebad

                    bad     = badsamples(D, chanind, ':', w);
                    Y(bad)  = NaN;

                end

                if weight

                    B1                  = spm_vec(B{i});
                    B1                  = repmat(B1,1,length(chanind),D.nsamples);
                    B1                  = permute(B1,[2 3 1]);
                    [Y, W1]             = spm_robust_barycentric(Y, 3, ks,B1);

                else

                    [Y, W1]             = spm_robust_average(Y, 3, ks);

                end

                W(chanind, :, :)    = W1;
                Dnew(chanind, :, i) = Y;

                if length(chanind)<D.nchannels
                    Dnew(setdiff(1:D.nchannels, chanind), :, i) = mean(D(setdiff(1:D.nchannels, chanind), :, w), 3);
                end

                if savew
                    Dw(:, :, w)   = W;
                end

            else

                B1              = repmat(spm_vec(B{i}),1,length(chanind),D.nsamples);
                B1              = permute(B1,[2 3 1]);
                X               = D(:, :, w);
                X(isnan(X))     = 0;
                Dnew(:, :, i)   = sum(W(:, :, w).*B1.*X, 3)./sum(B1.*W(:, :, w), 3);

            end

        end

        if ismember(i, Ibar),
            spm_progress_bar('Set', i);
        end

    end

else

    %-Averaging by channel
    %======================================================================
    spm_progress_bar('Init', D.nchannels, 'Channels completed');
    if D.nchannels > 100,
        Ibar = floor(linspace(1, D.nchannels, 100));
    else
        Ibar = [1:D.nchannels];
    end

    for j = 1:D.nchannels

        if robust && ~bycondition
            if ismember(j, chanind)

                Y       = D(j, :, goodtrials);

                if removebad
                    bad     = badsamples(D, j, ':', goodtrials);
                    Y(bad)  = NaN;
                end

                if weight

                    B1      = spm_vec(B);
                    B1      = permute(repmat(B1, 1, 1, D.nsamples), [1 3 2]); % 1 x nsamples x ntrials
                    [Y, W1] = spm_robust_barycentric(Y, 3, ks, B1);

                else

                    [Y, W1] = spm_robust_average(Y, 3, ks);

                end

                W                   = zeros([1 D.nsamples D.ntrials]);
                W(1, :, goodtrials) = W1;

            else

                W1 = ones(1, D.nsamples, length(goodtrials));

            end

            if savew
                Dw(j, :, goodtrials) = W1;
            end

        end

        for i = 1:D.nconditions

            w = indtrial(D, deblank(cl{i}), 'GOOD');

            if isempty(w)
                continue;
            end

            if ~robust || ~ismember(j, chanind)

                Y = D(j, :, w);

                if weight

                    B1 = spm_vec(B{i});
                    B1 = permute(repmat(B1, 1, 1, D.nsamples), [1 3 2]); % 1 x nsamples x ntrials
                    Y  = sum(B1 .* Y, 3) ./ sum(B1, 3);

                else

                    if S.trim > 0
                        Y = trimmean(Y, S.trim, 3);
                    else
                        Y = mean(Y, 3);
                    end

                end

                Dnew(j, :, i) = Y;

            else

                Y = D(j, :, w);

                if removebad
                    bad = badsamples(D, j, ':', w);
                    Y(bad) = NaN;
                end

                if bycondition
                    if weight
                        B1 = spm_vec(B{i});
                        B1 = permute(repmat(B1, 1, 1, D.nsamples), [1 3 2]);
                        [Y, W] = spm_robust_barycentric(Y, 3, ks, B1);
                    else
                        [Y, W] = spm_robust_average(Y, 3, ks);
                    end
                    Dnew(j, :, i) = Y;

                    if savew
                        Dw(j, :, w) = W;
                    end

                else
                    X = D(j, :, w);
                    X(isnan(X)) = 0;
                    Dnew(j, :, i) = sum(W(1, :, w) .* X, 3) ./ sum(W(1, :, w), 3);
                end
            end
        end

        if ismember(j, Ibar),
            spm_progress_bar('Set', j);
        end

    end
end

spm_progress_bar('Clear');

%-Update some header information
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, ':', dcl);
Dnew = repl(Dnew, ':', ni);
Dnew = montage(Dnew, 'switch', montage_ind);

%-Display averaging statistics
%--------------------------------------------------------------------------
fprintf('%s: Number of replications per contrast:\n', Dnew.fname);      %-#
for i = 1:D.nconditions
    fprintf('  average %s: %d trials\n', dcl{i}, ni(i));                 %-#
end

if robust
    disp('Robust averaging might have introduced high frequencies in the data. It is advised to re-apply low-pass filter.');
end

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

if robust && savew
    Dw = Dw.history(mfilename, S);
    save(Dw);
end

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG averaging: done'); spm('Pointer','Arrow');


function m = medianw( x, w, dim )
% Fast weighted median.
%
% Computes the weighted median of a set of samples.
%  http://en.wikipedia.org/wiki/Weighted_median
% "A weighted median of a sample is the 50% weighted percentile."
% For matrices computes median along each column (or dimension dim).
% If all weights are equal to 1 gives identical results to median.
%
% USAGE
%  m = medianw( x, w, [dim] )
%
% INPUTS
%  x      - vector or array of samples
%  w      - vector or array of weights
%  dim    - dimension along which to compute median
%
% OUTPUTS
%  m      - weighted median value of x
%
% EXAMPLE - simple toy example
%  x=[1 2 3]; w=[1 1 5]; medianw(x,w)
%
% EXAMPLE - comparison to median
%  n=randi(100); m=randi(100);
%  x=rand(n,m); w=ones(n,m);
%  m1=median(x); m2=medianw(x,w);
%  assert(isequal(m1,m2))
%
% See also median
%
% Piotr's Computer Vision Matlab Toolbox      Version 3.24
% Copyright 2014 Piotr Dollar.  [pdollar-at-gmail.com]
% Licensed under the Simplified BSD License [see external/bsd.txt]

if(nargin<3), 
    dim = find(size(x)~=1,1); 
end
d   = dim; 
nd  = ndims(x); 
n   = numel(x);

if( n==1 || size(x,d)==1 )

    m = x;

elseif (length(x)==n)

    [x,o]   = sort(x); 
    w       = w(o); 
    w       = cumsum(w);
    w       = w/w(end); 
    [~,j]   = min(w<=.5);

    if(j==1 || w(j-1)~=.5)
        m = x(j);
    else 
        m = (x(j-1)+x(j))/2; 
    end

else

    if (d>1) 
        p = [d 1:d-1 d+1:nd]; 
        x = permute(x,p); 
        w = permute(w,p); 
    end

    [x,o]       = sort(x);
    w           = w(o);
    w           = cumsum(w);
    is          = {':'};
    is          = is(ones(1,nd-1));
    w           = bsxfun(@rdivide,w,w(end,is{:}));
    [~,j]       = min(w<=.5);
    s           = size(x);
    s           = reshape(((1:n/s(1))-1)*s(1),size(j));
    j0          = max(1,j-1);
    j0          = j0+s; 
    j           = j+s;
    same        = w(j0)~=.5;
    j0(same)    = j(same);
    m           = (x(j0)+x(j))/2;
    if(d>1)
        p       = [2:d 1 d+1:nd];
        m       = permute(m,p);
    end
end

