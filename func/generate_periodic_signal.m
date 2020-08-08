
function x = generate_periodic_signal(Input_Periods,Input_Datalength,SNR,visual)

if nargin < 3
    SNR = [];
    visual = [];
end

if nargin < 4
    visual = [];
end
rng('default')

% Step 1: Selection of Input Signal Paparemters
% Input_Periods    = [3,5,7];  % Component Periods of the input signal
% Input_Datalength = 200;      % Data length of the input signal
% SNR              = 100;      % Signal to noise ratio of the input in dB


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

if visual
    figure;plot(x);
    title('Input Signal without Noise');
    xlabel('time index');
    ylabel('signal amplitude');
end



%if SNR
% Step 3: Adding Noise to the Input
lambda = 0.5;
ns = abs(randn(Input_Datalength,1));
ns = ns./norm(ns,2);
Desired_Noise_Power = 10^(-1*SNR/20);
ns = Desired_Noise_Power.*ns;
x  = lambda*x + (1-lambda)*ns; % signal + noise

if visual
    figure;plot(x);
    title('Input Signal With Noise');
    xlabel('time index');
    ylabel('signal amplitude');
end
%end




