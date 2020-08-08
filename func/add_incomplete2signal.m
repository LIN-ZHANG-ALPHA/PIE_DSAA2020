
function  x =  add_incomplete2signal(x,opts_syndata)



% if ~isfield(opts_syndata, 'incomplete'),          opts_syndata.incomplete          = 0; end  %
% if ~isfield(opts_syndata, 'ratio_incomplete'),    opts_syndata.ratio_incomplete    = 0.2; end  %
if ~isfield(opts_syndata, 'missing_window_size'), opts_syndata.missing_window_size = 1; end  % Signal to noise ratio of the input in dB
if ~isfield(opts_syndata, 'visual_incomplete'),   opts_syndata.visual_incomplete   = 0; end  % Signal to noise ratio of the input in dB

%if opts_syndata.incomplete
x = gen_incomplete(x,opts_syndata.Input_Datalength ,opts_syndata.ratio_incomplete, opts_syndata.missing_window_size);
x = x./norm(x,2);
if opts_syndata.visual_incomplete
    figure, plot(x);
    title('Input Signal With Missing data');
    xlabel('time index');
    ylabel('signal amplitude');
    %  set(h.a2, 'ButtonDownFcn',@buttonDownCallback)
end

