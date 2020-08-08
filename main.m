



clear all; % close all; % Clear all previous variables and close all previous windows.

% parameter for signals
opts_syndata.Input_Datalength   = 200;
opts_syndata.SNR                      = 5;
opts_syndata.Input_Periods      = {[3,7,11]};

% extra settings
opts_syndata.incomplete          = 0; % off: no, 1: on
opts_syndata.ratio_incomplete    = 0.3;
opts_syndata.missing_window_size = 1;
opts_syndata.visual_incomplete   = 0;
opts_syndata.visual                     = 0;
%% Step 0: Selection of Dictionary Parameters

Pmax            = [90,90]; %The largest period spanned by the NPDs
Dictionary_pool = {'Ramanujan','NaturalBasis','random' };%Type of the dictionary
Dictionary_type = Dictionary_pool{1};

opts.Dictionary_type = Dictionary_type;

opts.Pmax_LP     = Pmax(1);

%% step  1: genertate signal

x =  data_syn_generator(opts_syndata);

%% parameter settings

opts.Dictionary_type = Dictionary_type;
opts.Pmax            = Pmax;
opts.lambda_0        = 1;
opts.lambda_1        = 0.001;
opts.lambda_2        = 0.001;

opts.rho             = 1e-3;
opts.lp_show          = 0;
opts.max_iter        = 50;
opts.DIPS = 1;

%% main

tstart = tic;

[completed_x,beta_output,periods_vector_ours] = PIE(x,opts);

telapsed_lp = toc(tstart);

if (opts.DIPS)
    figure, stem(periods_vector_ours,'linewidth',3,'color',[0 0 0]);
    title('LP');
    xlabel('Period');
    ylabel('Strength');
end





















