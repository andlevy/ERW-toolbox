function [mse, da, y, XS, tau] = spm_era_fx_linear(A,X,XT,F,T,N)


[tau,B,G]   = spm_era_fx_tau(A,F);
tauF        = floor(tau);
tauC        = ceil(tau);

XF_index    = min(max(tauF, 1), T)';
XC_index    = min(max(tauC, 1), T)';
XF          = X;
XC          = X;
for n = 1:N
    XF(:,n) = X(XF_index(:,n),n);
    XC(:,n) = X(XC_index(:,n),n);
end
% XS          = (tau'-tauF').*XC + (tauC'-tau').*XF;
XS = [];
for n = 1:N
    XS(:,n) = interp1([1:T]',X(:,n),tau(n,:),'makima');
end

y           = mean(XS, 2);
e           = XT - XS;
mse         = mean2(e.^2);

% Calculate gradients
%--------------------------------------------------------------------------
dtau  = (-2/N)*(e.*(XC-XF))';
dtau  = dtau(:,2:end);
dp    = cumsum(dtau,2, 'reverse');
gsum  = sum(G, 2);
dg    = ((T-1)/T)*((dp.*repmat(gsum,1,T-1))-repmat(sum((dp.*G),2),1,T-1))./repmat(gsum.^2,1,T-1);
db    = dg .* sign(B);
da    = db * F';

% disp(['mse: ' num2str(mse)])
