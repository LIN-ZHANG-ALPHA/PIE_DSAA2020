 
 
function x  = get_incomplete_data(x,ratio_incomplete,missing_method)
% remove part of values from signal x
rng('default')
num_signal = length(x);
num_missing =  floor(ratio_incomplete*num_signal);
 
switch missing_method
    case 'random'
        missing_idx = randsample(num_signal,num_missing);
        x(missing_idx) = 0;
    case 'uniform'
        pace = floor(num_signal/num_missing);
        
        x(1:pace:end)=0;
        
    otherwise
        warning('Warning: No method specified!!!!.')
        
end
 
 
 
 