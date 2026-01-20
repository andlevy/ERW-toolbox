function [y, XS, tau,a1] = spm_era_estimate_hermite(X, XT,F, set)



try set.max_iter; catch, set.max_iter = 500; end
try set.order;    catch, set.order    = 8;   end
% 
% % Construct temporal basis set
% K           = set.order;
% [T, N]      = size(X);
% f           = (1:(K))' * (1:(T-1)) * pi / T;
% f           = sin(f);
% w           = tukeywin(T-1, 1);
% f           = (f'.*w)';

% Initialise
[T,N]   = size(X);
K       = size(F,1);
a0 = zeros(N, K);
opt = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', "bfgs", ...
    'FiniteDifferenceType', 'central', ...
    'StepTolerance', 5e-4, ...
    'ObjectiveLimit', 5e-5, ...
    'Display', 'final');
% opt     = optimoptions('fminunc', ...
%     'SpecifyObjectiveGradient',true, ...
%     'HessianApproximation',"bfgs", ...
%     'FiniteDifferenceType','central',...
%     'StepTolerance',5e-4, ...
%     'PlotFcn','optimplotfirstorderopt',...
%     'ObjectiveLimit',5e-5,...
%     'Display','iter');
a1 = fminunc(@(a) calc_mse_grad_hermite(a, X, XT, F, T, N), a0, opt);

% opt = optimoptions('fmincon', ...
%     'Algorithm', 'interior-point', ...
%     'SpecifyObjectiveGradient', true, ...
%     'HessianApproximation', 'lbfgs', ...
%     'Display', 'final', ...
%     'StepTolerance', 5e-4, ...
%     'OptimalityTolerance', 5e-5, ...
%     'MaxIterations', set.max_iter);

% Unconstrained bounds
% a1 = fmincon(@(a) calc_mse_grad_hermite(a, X, XT, f, T, N), a0, [], [], [], [], [], [], [], opt);


[~, ~, y, XS, tau] = forward_backward_hermite(a1, X, XT, F, T, N);

end
function [mse, da, y, XS, tau] = forward_backward_hermite(a, X, XT, f, T, N)

b = a * f + 1;
g = abs(b);
p = (T - 1) * g ./ sum(g, 2);
p = [ones(N, 1), p];

tau = cumsum(p, 2);
tau = correct_T(tau);

% Precompute Hermite interpolation indices and coefficients
k = floor(tau);
k = max(2, min(T - 2, k));
tt = tau - k;

h00  = (1 + 2 .* tt) .* (1 - tt).^2;
h10  = tt .* (1 - tt).^2;
h01  = tt.^2 .* (3 - 2 .* tt);
h11  = tt.^2 .* (tt - 1);

dh00 = 6 .* tt .* (tt - 1);
dh10 = 3 .* tt.^2 - 4 .* tt + 1;
dh01 = -dh00;
dh11 = 3 .* tt.^2 - 2 .* tt;

% Vectorised Hermite interpolation
rows = (1:T)';
cols = 1:N;
XS = zeros(T, N);
dXS_dtau = zeros(T, N);

for n = 1:N
    k_n   = k(n, :)';
    tt_n  = tt(n, :)';
    h00_n = h00(n, :)'; h10_n = h10(n, :)'; h01_n = h01(n, :)'; h11_n = h11(n, :)';
    dh00_n = dh00(n, :)'; dh10_n = dh10(n, :)'; dh01_n = dh01(n, :)'; dh11_n = dh11(n, :)';

    X0 = X(k_n, n);
    X1 = X(k_n+1, n);
    D0 = (X(k_n+1, n) - X(k_n-1, n)) / 2;
    D1 = (X(k_n+2, n) - X(k_n, n)) / 2;

    XS(:, n) = h00_n .* X0 + h10_n .* D0 + h01_n .* X1 + h11_n .* D1;
    dXS_dtau(:, n) = dh00_n .* X0 + dh10_n .* D0 + dh01_n .* X1 + dh11_n .* D1;
end

y   = mean(XS, 2);
e   = XS - XT;
mse = mean2(e.^2);

dMSE_dtau = (2 / (T * N)) * (e .* dXS_dtau);  % T x N

% d tau / d p = lower triangular matrix
L = tril(ones(T, T-1));
dp = (dMSE_dtau') * L;  % N x (T-1)

% d p / d g
dg = zeros(N, T-1);
for n = 1:N
    g_n  = g(n, :);
    gsum = sum(g_n);
    dg(n,:) = ((T - 1)/gsum^2) * (gsum * eye(T - 1) - g_n' * ones(1, T - 1)) * dp(n,:)';
end

% d g / d b
signb = sign(b);
db = dg .* signb;

% d b / d a
da = db * f';

end


% 
% function [mse, da, y, XS, tau] = forward_backward_hermite(a, X, XT, f, T, N)
% 
% b = a * f + 1;
% g = abs(b);
% p = (T - 1) * g ./ sum(g, 2);
% p = [ones(N, 1), p];
% 
% tau = cumsum(p, 2);
% tau = correct_T(tau);
% 
% XS = zeros(T, N);
% dXS_dtau = zeros(T, N);
% e = zeros(T, N);
% 
% for n = 1:N
%     for t = 1:T
%         tau_t = tau(n, t);
%         k = max(2, min(T - 2, floor(tau_t)));
%         h = 1;
%         tt = tau_t - k;
% 
%         y0 = X(k, n);
%         y1 = X(k+1, n);
%         d0 = (X(k+1, n) - X(k-1, n)) / 2;
%         d1 = (X(k+2, n) - X(k, n)) / 2;
% 
%         h00 = (1 + 2*tt)*(1 - tt)^2;
%         h10 = tt*(1 - tt)^2;
%         h01 = tt^2*(3 - 2*tt);
%         h11 = tt^2*(tt - 1);
% 
%         XS(t, n) = h00*y0 + h10*d0 + h01*y1 + h11*d1;
% 
%         dh00 = 6*tt*(tt - 1);
%         dh10 = 3*tt^2 - 4*tt + 1;
%         dh01 = -dh00;
%         dh11 = 3*tt^2 - 2*tt;
% 
%         dXS_dtau(t, n) = dh00*y0 + dh10*d0 + dh01*y1 + dh11*d1;
%     end
% end
% 
% y = mean(XS, 2);
% e = XS - XT;
% mse = mean2(e.^2);
% 
% dMSE_dtau = (2 / (T * N)) * (e .* dXS_dtau);  % T x N
% 
% % d tau / d p = lower triangular matrix
% L = tril(ones(T, T-1));
% dp = zeros(N, T - 1);
% for n = 1:N
%     dp(n, :) = dMSE_dtau(:, n)' * L;
% end
% 
% % d p / d g
% dg = zeros(N, T-1);
% for n = 1:N
%     g_n = g(n, :);
%     gsum = sum(g_n);
%     dg(n, :) = ((T - 1)/gsum^2) * (gsum * eye(T - 1) - g_n' * ones(1, T - 1)) * dp(n, :)';
% end
% 
% % d g / d b
% signb = sign(b);
% db = dg .* signb;
% 
% % d b / d a
% da = db * f';
% 
% end

function [mse, grad] = calc_mse_grad_hermite(a, X, XT, f, T, N)
[mse, grad, ~, ~, ~] = forward_backward_hermite(a, X, XT, f, T, N);
end

function T = correct_T(T)
TF = floor(T);
TC = ceil(T);
T = T + (TF == TC) .* (1e-5);
end

function w = tukeywin(N, alpha)
    if nargin < 2
        alpha = 0.5;
    end
    w = zeros(N, 1);
    for n = 1:N
        t = (n - 1) / (N - 1);
        if t < alpha / 2
            w(n) = 0.5 * (1 + cos(pi * (2*t/alpha - 1)));
        elseif t <= 1 - alpha / 2
            w(n) = 1;
        else
            w(n) = 0.5 * (1 + cos(pi * (2*t/alpha - 2/alpha + 1)));
        end
    end
end
%}