# Event-Related Warping (ERW) Toolbox

MATLAB toolbox for SPM that does temporal alignment and jitter correction of M/EEG data from sequential experimental paradigms.

In sequential tasks the timing of events within a trial (cues, targets, responses) varies from trial to trial, so ordinary epoch averaging smears the responses that are locked to later events. ERW estimates a temporal warp for each trial by aligning *template functions built from the event timings*, not the neural signals. The same warp is then applied to every channel, which keeps the timing relationships between channels intact.

## Preprint

> Levy, A. D., Zeidman, P., & Friston, K. (2026). **Event-Related Warping: A Toolbox for Temporal Alignment and Jitter Correction in Sequential Experimental Paradigms.** *bioRxiv*. https://doi.org/10.64898/2026.03.24.713943

Full text: https://www.biorxiv.org/content/10.64898/2026.03.24.713943v1

If you use this toolbox, please cite the preprint above.

## Requirements

- MATLAB
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/) on the MATLAB path (ERW works on SPM M/EEG objects)

## Installation

```matlab
addpath('/path/to/ERW');
addpath('/path/to/ERW/sim_tools');   % optional: simulation and figure tools
```

## Pipeline

Each step takes a single input structure `S`. See each function's help text for the full list of options.

| Step | Function | Purpose |
|------|----------|---------|
| 1 | `spm_eeg_era_design` | Define the trial design from within-trial events (onsets and their distributions) in an epoched MEEG file |
| 2 | `spm_eeg_era_template` | Build temporal warping templates for each subject, or a group template |
| 3 | `spm_eeg_era_estimate` | Estimate warping trajectories between the design and the template (`S.method`: `'linear'`, `'hermite'`, `'parabolic'`, `'sinc'`) |
| 4 | `spm_eeg_era_warp` | Apply the estimated warps to the M/EEG time series |
| 5 | `spm_eeg_era_average` | Average over trials, with optional robust, trimmed or barycentric weighting |

Minimal example:

```matlab
% 1. Design
S            = [];
S.D          = 'epoched_data.mat';
S.designdef(1).conditionlabel       = 'go';
S.designdef(1).onset(1).eventtype   = 'cue';     % first onset is trial-locked (t = 0)
S.designdef(1).onset(1).eventvalue  = 1;
S.designdef(1).onset(2).eventtype   = 'target';
S.designdef(1).onset(2).eventvalue  = 2;
ERA          = spm_eeg_era_design(S);

% 2. Template
S            = [];
S.ERA        = ERA;
ERA          = spm_eeg_era_template(S);

% 3. Estimate warps
S            = [];
S.ERA        = ERA;
S.method     = 'parabolic';
ERA          = spm_eeg_era_estimate(S);

% 4. Warp data
S            = [];
S.D          = 'epoched_data.mat';
S.ERA        = ERA;
ERA          = spm_eeg_era_warp(S);

% 5. Average
S            = [];
S.D          = 'warped_data.mat';                 % output of step 4
S.ERA        = ERA;
D            = spm_eeg_era_average(S);
```

Event types and values, and the output file names, depend on your data.

## Validation

- `spm_eeg_era_validate_erp`: RMSE between ERW-aligned averages and a reference ERP
- `spm_eeg_era_validate_xcorr_corrected`: checks that inter-channel lags are preserved (ERW-aligned vs. reference averages). This replaces `spm_eeg_era_validate_xcorr`, whose channel-pair mask kept the weakly correlated pairs it was meant to remove. See the header of the corrected file for details.

## Simulation tools (`sim_tools/`)

Functions for generating synthetic transient signals and noise, simulating jittered trials, evaluating alignment metrics, and making summary figures (e.g. `spm_era_gen_transient`, `spm_era_gen_noise`, `spm_era_sim_transient`, `spm_era_sim_evalmetrics`, `spm_era_sim_results_figure`).

## Acknowledgements

- `spm_era_estimate_ttw` / `spm_era_estimate_fltw` are adapted from the Trainable Time Warping code by Soheil Khorram: Khorram, S., McInnis, M. G., & Provost, E. M. (2019). *Trainable Time Warping: Aligning Time-Series in the Continuous-Time Domain.* arXiv:1903.09245.
- `spm_robust_barycentric` is adapted from SPM's `spm_robust_average` (James Kilner, Wellcome Trust Centre for Neuroimaging).

## Contact

Andrew D. Levy
