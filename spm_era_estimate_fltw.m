function [y, XS, tau,a1] = spm_era_estimate_fltw(X,XT,F,set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast learnable time warping algorithm for DTW averaging
% Inputs
%     X: input time-series; it is an N-by-T matrix where N is the
%        number of time-series and T is the length of them. Note that
%        this algorithm assumes all time-series have the same length.
%        If your time-series have different lengths, use an interpolation
%        technique (e.g., spline function) to make them equi-length.
%     order: order of the DST used to model the warping functions. With
%            higher order, we consider more details of the input signals
%            and we learn more complex warping functions. However, a
%            high order value will cause a more complex optimization
%            landscape, which increases the likelihood of converging to
%            a weak local optimum point.
%    max_iter_num: maximum number of iterations in our gradient-based
%                  optimization algorithm
%
% Outputs
%    y: final average signal. It is a 1-by-T vector.
%    XS: synchronized time-series. y is equal to mean(XS, 1).
%    tau: learned warping functions. It is an N-by-T matrix. tau(n, :)
%         is the warping function used to obtain XS(n, :) from X(n, :).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Soheil Khorram (khorram.soheil@gmail.com)
% Please feel free to contact me if you have any question.
%
% @article{khorram2019trainable,
%   title={Trainable Time Warping: Aligning Time-Series in the Continuous-Time Domain},
%   author={Khorram, Soheil and McInnis, Melvin G and Provost, Emily Mower},
%   journal={arXiv preprint arXiv:1903.09245},
%   year={2019}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try set.max_iter; catch, set.max_iter = 500; end

% Setup Optimisation
%--------------------------------------------------------------------------
[T,N]   = size(X);
K       = size(F,1);
a0      = zeros(N, K);
% opt     = optimoptions('fminunc', ...
%     'SpecifyObjectiveGradient',true, ...
%     'HessianApproximation',"bfgs", ...
%     'FiniteDifferenceType','central',...
%     'StepTolerance',5e-4, ...
%     'PlotFcn','optimplotfirstorderopt',...
%     'ObjectiveLimit',5e-5,...
    % 'Display','iter');
opt     = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianApproximation',"bfgs", ...
    'FiniteDifferenceType','central',...
    'StepTolerance',5e-4, ...
    'PlotFcn',[],...
    'ObjectiveLimit',5e-5,...
    'Display','final');

a1             = fminunc(@(a)calc_mse_grad(a,X,XT,F,T,N),a0,opt);
[~,~,y,XS,tau] = spm_era_fx_linear(a1,X,XT,F,T,N);

end

% function [mse, da, y, XS, tau] = forward_backward(a,X,XT,F,T,N)
% 
% [tau,B,G]   = spm_era_fx_tau(A,F);
% tauF        = floor(tau);
% tauC        = ceil(tau);
% 
% XF_index    = min(max(tauF, 1), T)';
% XC_index    = min(max(tauC, 1), T)';
% XF          = X;
% XC          = X;
% for n = 1:N
%     XF(:,n) = X(XF_index(:,n),n);
%     XC(:,n) = X(XC_index(:,n),n);
% end
% XS          = (tau'-tauF').*XC + (tauC'-tau').*XF;
% % 
% % XS = zeros(T,N);
% % for n = 1:N
% %     XS(:,n) = interp1(1:T, X(:,n), tau(n,:), 'makima');
% % end
% 
% y           = mean(XS, 2);
% e           = XT - XS;
% mse         = mean2(e.^2);
% 
% % Calculate gradients
% %--------------------------------------------------------------------------
% dtau  = (-2/N)*(e.*(XC-XF))';
% dtau  = dtau(:,2:end);
% dp    = cumsum(dtau,2, 'reverse');
% gsum  = sum(g, 2);
% dg    = ((T-1)/T)*((dp.*repmat(gsum,1,T-1))-repmat(sum((dp.*g),2),1,T-1))./repmat(gsum.^2,1,T-1);
% db    = dg .* sign(b);
% da    = db * F';
% 
% % disp(['mse: ' num2str(mse)])
% 
% end

function [mse, grad] = calc_mse_grad(a,X,XT,F,T,N)
[mse,grad,~,~,~] = spm_era_fx_linear(a,X,XT,F,T,N);
end


