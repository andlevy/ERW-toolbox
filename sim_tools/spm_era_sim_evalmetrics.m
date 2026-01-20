
function E = spm_era_sim_evalmetrics(GT, A,W)

% fs       = 500;
% lpf_cutoff = 50;
% hpf_cutoff = .1;
% N        = size(GT,1);
% f       = (0:N-1)' * (fs / (2*N));
% mask    = f <= lpf_cutoff & f >= hpf_cutoff;
% D       = dctmtx(N);
% 
% GT      = D * GT;             % Forward DCT
% GT      = GT .* mask;
% GT      = D' * GT;
% GT      = GT - mean(GT(1:10,:),1);

% Extract the field names from the structure A, representing different models
fn                      = fieldnames(A);
fn(strcmp(fn,'wlabel')) = [];
small_value = 0.01 * mean(abs(GT(:)));  % 1% of the average of GT

if isempty(W)
    W = ones(length(GT),1);
else
    W = double(W(:));
end

for f = 1:numel(fn)
    pred  = W.*A.(fn{f});
    truth = W.*GT(:);
    
    % Calculate Root Mean Squared Error
    E.rms.(fn{f})   = sqrt(mean((truth - pred).^2));
    
    % Calculate Coefficient of Variation of RMSE
    E.cvrms.(fn{f}) = E.rms.(fn{f}) / mean(truth);
    
    % Calculate Maximum Absolute Deviation
    E.mad.(fn{f})   = max(abs(truth - pred));
    
    % Calculate R^2 (Coefficient of Determination)
    ss_res          = sum((truth - pred).^2);
    ss_tot          = sum((truth - mean(truth)).^2);
    E.r2.(fn{f})    = 1 - (ss_res / ss_tot);
    
    % Calculate MAPE (Mean Absolute Percentage Error)
    E.mape.(fn{f}) = 100 * mean(abs(truth - pred) ./ (abs(truth) + small_value));
    
    % Calculate MAE (Mean Absolute Error)
    E.mae.(fn{f}) = mean(abs(truth - pred));
    
    % Calculate Pearson correlation coefficient
    r_numerator   = sum((truth - mean(truth)) .* (pred - mean(pred)));
    r_denominator = sqrt(sum((truth - mean(truth)).^2) * sum((pred - mean(pred)).^2));
    E.pearsonr.(fn{f}) = r_numerator / r_denominator;

     % Calculate L2 norm difference
    E.l2.(fn{f}) = norm(truth - pred, 2);

    % Calculate normalised L2 norm
    E.l2norm.(fn{f}) = norm(truth - pred, 2) / (norm(truth, 2) + eps);  % add eps to avoid division by zero

    % Normalised RMSE by L2 norm of ground truth
    E.nrmse.(fn{f}) = E.rms.(fn{f}) / std(truth);
end

end

