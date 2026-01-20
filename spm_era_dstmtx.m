function F = spm_era_dstmtx(T,order,tukey)

if nargin < 3
    tukey = 0;
end

% Construct temporal basis set
%--------------------------------------------------------------------------
K       = order;
% [T,N]   = size(X);
f       = (1:(K))' * (1:(T-1)) * pi / (T);
f       = sin(f);

% Apply Tukey window
%--------------------------------------------------------------------------
w       = tukeywin(T-1,tukey);
 % w     = gausswin(T-1, 4);
F       = (f'.*w)';

function w = tukeywin(N, alpha)
if nargin < 2
    alpha = 0.5;
end
w = zeros(N,1);
for n = 1:N
    t = (n - 1) / (N - 1);
    if t < alpha/2
        w(n) = 0.5*(1 + cos(pi*(2*t/alpha - 1)));
    elseif t <= 1 - alpha/2
        w(n) = 1;
    else
        w(n) = 0.5*(1 + cos(pi*(2*t/alpha - 2/alpha + 1)));
    end
end
