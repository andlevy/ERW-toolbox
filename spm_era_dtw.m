function [D, dtwMatrix, path] = spm_era_dtw(X, Y, w)
    % X and Y are the two time series vectors.
    % w is the optional window size for local restriction of DTW calculation.

    % Get the lengths of the time series
    n = length(X);
    m = length(Y);

    % Check if the window size is given, otherwise set it large enough to cover the whole matrix
    if nargin < 3 || isempty(w)
        w = max(n, m); % default: no window restriction
    end

    % Initialize the DTW matrix with infinite values
    dtwMatrix = nan(n+1, m+1);
    dtwMatrix(1, 1) = 0;

    % Populate the DTW matrix within the window constraints
    for i = 2:n+1
        for j = max(2, i-w):min(m+1, i+w)
            cost = abs(X(i-1) - Y(j-1));
            localMin = min([dtwMatrix(i-1, j), dtwMatrix(i, j-1), dtwMatrix(i-1, j-1)]);
            dtwMatrix(i, j) = cost + localMin;
        end
    end

    % The DTW distance is the bottom right corner of the matrix
    D = dtwMatrix(n+1, m+1);

    % Traceback to find path
    i = n + 1;
    j = m + 1;
    path = [i, j];

    while i > 1 && j > 1
        [~, idx] = min([dtwMatrix(i-1, j), dtwMatrix(i, j-1), dtwMatrix(i-1, j-1)]);
        switch idx
            case 1
                i = i - 1;
            case 2
                j = j - 1;
            case 3
                i = i - 1;
                j = j - 1;
        end
        path = [i j; path];
    end

    % Adjust path to proper index (remove the initial +1 offset)
    path = path - 1;
end
