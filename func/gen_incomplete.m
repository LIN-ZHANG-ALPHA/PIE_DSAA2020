
function x_raw = gen_incomplete(x_raw, len_signal,ratio_burst, missing_window_size)

if nargin< 3
    missing_window_size = 1;
end

rng('default')

%%

% x_tmp =  zeros(len_signal,1);

len_incomplete = floor(ratio_burst*len_signal);
num_window     = round(len_incomplete/missing_window_size);

num_candidates = floor(len_signal/missing_window_size);

bursty_idx   = randsample([1:num_candidates],num_window);

for i = 1: num_window
    cur_end_position   = min(bursty_idx(i)*missing_window_size,len_signal); % make sure within x length
    x_raw(cur_end_position-missing_window_size+1:cur_end_position) =  0;
end




end