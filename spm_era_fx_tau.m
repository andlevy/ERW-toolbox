function [tau,B,G,P] = spm_era_fx_tau(A,F)

B           = A * F + 1;
G           = abs(B);
P           = size(F,2) * G ./ sum(G,2);
P           = [ones(size(A,1),1), P];
tau         = cumsum(P,2);
tau         = correct_T(tau);

function T = correct_T(T)
TF = floor(T);
TC = ceil(T);
T  = T + (TF == TC) .* (1e-5);
