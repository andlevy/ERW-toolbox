


function [Y, pst] = spm_era_gen_transient(trans)
% spm_era_gen_transient - Generate synthetic transient signals composed of
% time-locked components (e.g., sine bursts, exponentials, sigmoids)
%
% FORMAT: [Y, pst] = spm_era_gen_transient(trans)
%
% INPUT:
%   trans - structure array defining transient waveforms.
%           Each element describes a simulated trial and must contain:
%             .pst        - [start, end] peri-stimulus time window (s)
%             .dt         - sampling interval (s)
%             .elem       - array of transient elements with fields:
%                 .type      - 'sin' | 'exp' | 'sig'
%                 .set       - structure with parameters:
%                     .onset   - onset time (s)
%                     .twin    - [start, end] duration rel. to onset (s)
%                     .amp     - amplitude
%                     .hz      - frequency (for 'sin')
%                     .tau     - time constant (for 'exp')
%                     .k, .t0  - slope & inflection point (for 'sig')
%
% OUTPUT:
%   Y   - ns × ntr matrix of synthetic signals (one column per trial)
%   pst - 1 × ns peri-stimulus time vector (s)
% Assume all trials share the same pst and dt
dt  = trans(1).dt;
pst = trans(1).pst(1) : dt : trans(1).pst(2);
ns  = numel(pst);
ntr = numel(trans);

Y = zeros(ns, ntr);


for i = 1:ntr
    y = zeros(ns, 1);

    for e = 1:numel(trans(i).elem)

        S     = trans(i).elem(e).set;
        type  = trans(i).elem(e).type;
        twin  = S.twin;
        if ~any(round(pst,3) == round(round(S.onset/dt)*dt,3))
            error('onset time not aligned with sampling rate')
        end
        onset = round(S.onset/dt)*dt;

        % 1) Build mask on global pst
        stim = (pst >= (onset + twin(1))) ...
            & (pst <= (onset + twin(2)));
        idx  = find(stim);           % e.g. [235:286] for a 52-sample window
        seg_t = pst(idx) - onset;  % perfectly multiples of dt

        % 1) template at zero
        % seg_t = round(twin(1)/dt)*dt:dt:round(twin(2)/dt)*dt;
        seg   = generate_waveform(seg_t,type, S);  % your sin/exp/sig switch, as column
        % L     = numel(seg);

        % 4) Apply continuity / carry‑in logic
        switch type
            case 'sin'
                seg = S.amp*seg + y(idx(1));
            case 'sig'
                baseline = y(idx(1)-1);
                % seg is currently in [0,1]; scale it by S.amp and add baseline
                seg = baseline + S.amp * seg;

                %   x1  = y(idx(1)-1);
                % x2  = y(idx(end)) + S.amp;
                % seg = x1 + (x2 - x1) * seg;
            case 'exp'
                % nothing extra
        end
        % plot(seg_t,seg)
        % 5) Insert into y
        x = zeros(ns,1);
        x(idx) = seg;
        y(idx) = 0;    % clear any previous energy there
        y      = y + x;

        % % 2) find onset and offset index
        % idx_ons = round(((twin(1) +onset) - pst(1))/dt);
        % idx_off = round(((twin(2) +onset) - pst(1))/dt);
        % if idx_off > ns
        %     trim    = idx_off - ns;
        %     L       = L - trim;
        %     seg     = seg(1:L);
        %     seg_t   = seg_t(1:L);
        %     idx_off = ns;
        % end
        %
        % % plot(seg_t,seg)
        % if length(idx_ons:idx_off) ~= L
        %     error('err')
        % end
        %
        % % 3) make stim mask in global pst
        % stim                  = false(ns,1);
        % stim(idx_ons:idx_off) = true;   % (you may wish to guard bounds)
        % % plot(pst,stim)
        %
        % % 4) build x
        %
        % switch type
        %     case 'sig'
        %         x1   = y(idx_ons-1);
        %         x2   = y(idx_off) + S.amp;
        %         seg  = x1 + (x2 - x1) * seg;
        %     case 'exp'
        %         seg      =  S.amp*seg;
        %     case 'sin'
        %         baseline = y(idx_ons);
        %         seg      =  S.amp*seg + baseline;
        % end
        %
        % x        = zeros(ns,1);
        % x(stim)  = seg;
        % y(stim) = 0;
        % y = y +x;

    end
    Y(:, i) = y;
end


function seg = generate_waveform(t,type, S)
switch type
    case 'sin'
        seg =  sin(2*pi*S.hz * t);
    case 'exp'
        seg =  (exp(-t/S.tau) - 1);
    case 'sig'
        seg =  (1./(1+exp(-S.k*(t-S.t0))));
end
seg = seg(:);

% %{
%
% Y = [];
% for i = 1:numel(trans)
%
%     % ns   = 1+(trans(i).pst(end)-trans(i).pst(1)) ./trans(1).dt;
%     % pst  =  (0:(ns-1)) *trans(1).dt + trans(i).pst(1);
%     pst  = trans(i).pst(1):trans(i).dt:trans(i).pst(2);
%     ns   = length(pst);
%     y    = zeros(ns,1);
%
%     for e = 1:numel(trans(i).elem)
%
%         ons  = trans(i).elem(e).set.onset;
%         t    = pst-ons;
%         stim = spm_vec(t>=trans(i).elem(e).set.twin(1)&t<=trans(i).elem(e).set.twin(2));
%
%         switch trans(i).elem(e).type
%             case 'sin'
%
%                 x  = spm_vec(sin(2*pi*trans(i).elem(e).set.hz*t));
%                 x  = x.*trans(i).elem(e).set.amp;
%                 x1 = y(find(stim,1,'first'));
%                 x  = x + x1;
%
%             case 'exp'
%
%                 x = spm_vec(exp(-t/trans(i).elem(e).set.tau))-1;
%                 x = x.*trans(i).elem(e).set.amp;
%
%             case 'sig'
%
%                 x  = spm_vec(1./(1+exp(-trans(i).elem(e).set.k*(t-trans(i).elem(e).set.t0))));
%                 x1 = y(find(stim,1,'first')-1);
%                 x2 = y(find(stim,1,'last')) + trans(i).elem(e).set.amp;
%                 x  = x1 + (x2 - x1) * x;
%
%
%         end
%         x = x.*stim;
%         y = y.*~stim;
%         y = y + x;
%     end
%     Y(:,i) = y;
% end
% %}
