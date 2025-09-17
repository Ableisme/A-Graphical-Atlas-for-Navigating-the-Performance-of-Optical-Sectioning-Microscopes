
function [FR, prof, centers] = radial_profile(fx, fy, M, nBins)
%RADIAL_PROFILE Radial average of 2D map M over frequency axes fx, fy
if nargin<4, nBins = 200; end
[FX,FY] = meshgrid(fx, fy);
FR = sqrt(FX.^2 + FY.^2);
fr_min = 0; fr_max = min(max(abs(fx)), max(abs(fy)));
edges = linspace(fr_min, fr_max, nBins+1);
centers = 0.5*(edges(1:end-1)+edges(2:end));
prof = zeros(numel(centers),1);
for b = 1:numel(centers)
    mask = (FR>=edges(b)) & (FR<edges(b+1));
    vals = M(mask);
    if isempty(vals)
        prof(b) = NaN;
    else
        prof(b) = mean(vals(:),'omitnan');
    end
end
% fill NaNs
nanid = isnan(prof);
if any(nanid)
    prof(nanid) = interp1(find(~nanid), prof(~nanid), find(nanid), 'linear', 'extrap');
end
end
