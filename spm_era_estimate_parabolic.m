function [y, XS, tau,a1,mse] = spm_era_estimate_parabolic(X,XT,F,set)

try set.max_iter; catch, set.max_iter = 500; end

% Setup Optimisation
%--------------------------------------------------------------------------
[T,N]   = size(X);
K       = size(F,1);
a0      = zeros(N, K);
opt     = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianApproximation','bfgs', ...
    'FiniteDifferenceType','central',...
    'StepTolerance',5e-4, ...
    'ObjectiveLimit',5e-5, ...
    'Display','final');

a1               = fminunc(@(a)calc_mse_grad_parabolic(a,X,XT,F,T,N),a0,opt);
[mse,~,y,XS,tau] = spm_era_fx_parabolic(a1,X,XT,F,T,N);

end


function [mse, grad] = calc_mse_grad_parabolic(a,X,XT,f,T,N)
[mse, grad, ~, ~, ~] = spm_era_fx_parabolic(a,X,XT,f,T,N);
end

