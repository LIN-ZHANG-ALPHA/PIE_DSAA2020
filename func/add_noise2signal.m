

function x = add_noise2signal(x,opts_syndata)

% add noise to periodic signal

SNR =  opts_syndata.SNR;

Input_Datalength = opts_syndata.Input_Datalength;
% Step 3: Adding Noise to the Input
ns = abs(randn(Input_Datalength,1));
% ns = randn(Input_Datalength,1);
ns = ns./norm(ns,2);
Desired_Noise_Power = 10^(-1*SNR/20);
ns = Desired_Noise_Power.*ns;
x  = x + ns; % signal + noise

if opts_syndata.visual
    figure;plot(x);
    title('Input Signal With Noise');
    xlabel('time index');
    ylabel('signal amplitude');
end
%end