

function x =  get_periodic_signal(opts_syndata)

rng('default')
% genertate perfect signal with input periods
% opts_syndata.Input_Datalength = 200;


%%


% Step 0: Selection of Dictionary Parameters
%
% Pmax = 50; %The largest period spanned by the NPDs
% Dictionary_type = 'Ramanujan'; %Type of the dictionary
%


% --- signal parameters ++++++++++++++
if ~isfield(opts_syndata, 'Input_Datalength'),  opts_syndata.Input_Datalength   = 200; end  %
if ~isfield(opts_syndata, 'Input_Periods'),     opts_syndata.Input_Periods      = {[3,7,11]}; end
% if ~isfield(opts_syndata, 'SNR'),               opts_syndata.SNR                = 10; end  % Signal to noise ratio of the input in dB
% if ~isfield(opts_syndata, 'visual_incomplete'), opts_syndata.visual_incomplete  = 0; end
% if ~isfield(opts_syndata, 'visual_signal'),     opts_syndata.visual_signal      = 0; end


%%
% x   = generate_periodic_signal(opts_syndata.Input_Periods{1},opts_syndata.Input_Datalength,opts_syndata.SNR);
% x = generate_periodic_signal(Input_Periods,Input_Datalength,SNR,visual)

% Step 1: Selection of Input Signal Paparemters
% Input_Periods    = [3,5,7];  % Component Periods of the input signal
% Input_Datalength = 200;      % Data length of the input signal
% SNR              = 100;      % Signal to noise ratio of the input in dB


Input_Periods    = opts_syndata.Input_Periods{1};
Input_Datalength = opts_syndata.Input_Datalength;

% Step 2: Generating a random input periodic signal
np = length(Input_Periods);    % number of base period
x  = zeros(Input_Datalength,1);

for i = 1:np
    x_temp = abs(randn(Input_Periods(i),1));  % for each period
    x_temp = repmat(x_temp,ceil(Input_Datalength/Input_Periods(i)),1); % repeated the number of interval for the current period
    x_temp = x_temp(1:Input_Datalength); % truncate to make the length is the same as the desired lenght
    x_temp = x_temp./norm(x_temp,2);
    x      = x + x_temp;                 % add signal with different periods together
end
% x = x./norm(x,2) + 1; % make sure all positive
x = x./norm(x,2);

if opts_syndata.visual
    figure;plot(x);
    title('Input Periodic Signal without Noise');
    xlabel('time index');
    ylabel('signal amplitude');
end


