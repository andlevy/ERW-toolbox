function E = spm_era_gen_noise(noise,Y)


% Create filter
%---------------
fs           = 500;
lpf_cutoff   = 50;
[b_lp, a_lp] = butter(4, lpf_cutoff / (fs/2));
signal_std   = std(Y,0,1);
[Nt,Nc]      = size(Y);

N = {}; k  = numel(noise);
for n = 1:k

    if noise(n).snr == 0

        N{n} = zeros(size(Y));

    else
        
        % Generate noise
        %------------------------------------------------------------------
        switch lower(noise(n).type)

            case 'pink'

                % b       = [0.049922, -0.095993, 0.050612, -0.004408];
                % a       = [1, -2.494956, 2.017265, -0.522189];
                % e       = randn(size(Y'));
                % N{n}    = filter(b, a, e, [], 2)';

                % Build the 1/sqrt(f) gain vector for FFT-based pink noise
                f                = (0:Nt-1)' * (fs/Nt);
                m                = floor(Nt/2);
                H                = zeros(Nt,1);
                H(2:m)           = 1 ./ sqrt(f(2:m));        % positive freqs up to Nyquist
                H(m+1)           = 1 ./ sqrt(f(m+1));        % Nyquist or nearest bin
                H(Nt-m+2 : Nt)   = flipud(H(2:m));            % mirror for negative freqs
                H(1)             = 0;                        % DC component zero

                % FFT-shaped pink noise
                W = fft(randn(Nt, Nc));          % white in freq domain
                W = bsxfun(@times, W, H);          % apply 1/sqrt(f) shaping
                N{n} = real(ifft(W));               % back to time domain


            case 'ar1'

                phi     = 0.98;
                e       = randn(size(Y')) * 1;
                N{n}    = filter(1, [1, -phi], e, [], 2)';

        end

  
        % % Apply Filter
        % %-------------
        % N{n} = filter(b_lp, a_lp, N{n});

        % Scale noise
        %------------------------------------------------------------------
        noise_std   = std(N{n},0,1);
        scale       = (signal_std ./ noise_std) / (noise(n).snr * sqrt(k));
        N{n}        = N{n} .* scale;  % scale per channel

    end
  

end

% Add noise
%------------------------------------------------------------------
E = zeros(size(Y));
for n = 1:numel(N)
    E = E + N{n};
end

E = filter(b_lp, a_lp, E);

% Verify SNR
%--------------------------------------------------------------------------
SNR = std(Y,0,1) ./ std(E,0,1);
disp('Amplitude SNR per channel:');
disp(SNR);


