
function [fc_thresh, fc_mtf10, fc_mtf50] = mtf_cutoff_linear(freq, mtf, thresh)
if nargin<3, thresh = 1e-3; end
freq = freq(:); mtf = mtf(:);
[~,imax] = max(mtf);
x = freq(imax:end); y = mtf(imax:end);
fc_thresh = crossing(x,y,thresh);
fc_mtf10  = crossing(x,y,0.10);
fc_mtf50  = crossing(x,y,0.50);
end

function fc = crossing(x,y,level)
fc = NaN;
idx = find(y<=level,1,'first');
if isempty(idx) || idx<2, return; end
fc = interp1(y([idx-1 idx]), x([idx-1 idx]), level, 'linear');
end
