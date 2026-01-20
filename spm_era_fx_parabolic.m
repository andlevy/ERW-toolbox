function [mse, da, y, XS, tau] = spm_era_fx_parabolic(A,X,XT,F,T,N)


[tau,B,G]   = spm_era_fx_tau(A,F);
kF          = floor(tau);
kF          = max(2, min(T-1, kF));
tt          = tau - kF;
XS          = zeros(T,N);
dXS_dtau    = zeros(T,N);
for n = 1:N
    for t = 1:T
        k             = kF(n,t);
        y_1           = X(k-1,n);
        y0            = X(k,n);
        y1            = X(k+1,n);

        a_            = 0.5*(y1 - 2*y0 + y_1);
        b_            = 0.5*(y1 - y_1);
        c_            = y0;

        t_            = tt(n,t);
        XS(t,n)       = a_*t_^2 + b_*t_ + c_;
        dXS_dtau(t,n) = 2*a_*t_ + b_;
    end
end

y   = mean(XS,2);
e   = XT - XS;
mse = mean2(e.^2);

% Gradient
%--------------------------------------------------------------------------
dtau = (-2/N) * (e .* dXS_dtau)';
dtau = dtau(:,2:end);
dp   = cumsum(dtau,2,'reverse');
gsum = sum(G,2);
dg   = ((T-1)/T)*((dp.*gsum) - sum(dp.*G,2))./gsum.^2;
db   = dg .* sign(B);
da   = db * F';

function T = correct_T(T)
TF = floor(T);
TC = ceil(T);
T  = T + (TF == TC) .* (1e-5);


