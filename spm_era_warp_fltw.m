
function Yw = spm_era_warp_fltw(Y,Tau)%,alpha,A,F)

% if nargin <3
%     alpha = 0; % no gain correction
% end


nChans               = size(Y,1);
[P, nTrials, nSteps] = size(Tau);
% TauFinal             = zeros(P, nTrials);
% G                    = zeros(P, nSteps, nTrials, 'like', Tau);   % gain for inspection
% Yw                   = zeros(size(Y),            'like', Y);      % warped result
% grid                 = (1:P)';   % column vector
%
%
%
% % Combination of warp to a single trajectory
% %------------------------------------------------------------------------
% combine = true;
% if combine
%     [Tau,~,~,J]  = spm_era_fx_tau(sum(A,3),F);
%     Tau          = Tau';
%     J            = J';
%     V            = var(J,[],2);
%     V            = V./max(V);
%     J            = (J-1).*V +1;
%
%
%     nSteps       = 1;
%     G            = ones(P, nSteps, nTrials, 'like', Tau);   % gain for inspection
% end
%
%
% % Calculate gain
% %--------------------------------------------------------------------------
%
% % Smoothing kernel
% %------------------
% s1       = 12;
% w1       = ceil(6*s1)+1;
% x1       = -(w1-1)/2 : (w1-1)/2;
% g1       = exp(-0.5*(x1/s1).^2);
% g1       = g1 / sum(g1);
%
% for t = 1:nTrials
%
%     % dτ/dx
%     % tauMat = squeeze(Tau(:,t,:));
%     % JMat = (gradient(tauMat'))' ;
%     JMat  = squeeze(J(:,t,:));
%
%     % Gaussian smoothing & extrema management
%     JgMat = bsxfun(@minus, JMat, JMat(1,:)) + 1;
%     JgMat = imfilter(JgMat, g1(:), 'symmetric', 'conv', 'same');
%     d     = .125;
%     JgMat = 1+ d*tanh((JgMat-1)./d);
%     JgMat = bsxfun(@rdivide, JgMat, mean(JgMat));
%
%     % Calculate gain
%     gainMat = (JgMat.^-alpha);
%     gainMat = bsxfun(@rdivide, gainMat, mean(gainMat));     % normalise
%
%     % Store
%     G(:,:,t) = gainMat;
%
% end
%
% % Apply warp to timeseries with optional gain correction
% %--------------------------------------------------------------------------
% for t = 1:nTrials
%     for c = 1:nChans
%         XS           = squeeze(Y(c,:,t));
%         for s = 1:nSteps
%             tau         = Tau(:,t,s);
%             XS          = G(:,s,t).*interp1(grid, XS, tau, 'spline', 'extrap');
%             % XS          = sign(XS) .* (abs(XS) .* G(:,s,t)); % gain only on magnitude
%         end
%
%         Yw(c,:,t)   = XS';
%
%     end
% end

% apply warp to timeseries
%----------------------------------------------------------------------
Yw = Y;
for s = 1:nSteps
    for t = 1:nTrials

        tau  = Tau(:,t,s);
        % dtau = gradient(tau);
        % dtau = bsxfun(@minus, dtau, dtau(1,:)) + 1;
        % dtau = 1+ d*tanh((dtau-1)./d);
        % dtau = bsxfun(@rdivide, dtau, mean(dtau));
        % gain = dtau.^alpha;
        % gain = gain./mean(gain);
        gain= 1;
        for c = 1:nChans

            X           = squeeze(Yw(c,:,t))';
            T           = size(X,1);
            % tauF        = floor(tau);
            % tauC        = ceil(tau);
            % XF_index    = min(max(tauF, 1), T);
            % XC_index    = min(max(tauC, 1), T);
            % XF          = X(XF_index);
            % XC          = X(XC_index);
            % YS          = (tau-tauF).*XC + (tauC-tau).*XF;

            % YS = interp1(1:T, X, tau, 'pchip');
            YS          = gain.*interp1(1:T, X, tau, 'makima');
            Yw(c,:,t)   = YS';

        end
    end
end



