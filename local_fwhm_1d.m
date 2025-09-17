
function w = local_fwhm_1d(x, y)
x = x(:); y = y(:);
[ym,im] = max(y);
if isempty(im) || ~isfinite(ym) || ym<=0, w = NaN; return; end
half = ym/2;
% left crossing
iL = find(y(1:im) <= half, 1, 'last');
if isempty(iL) || iL==im, xL = x(1);
else, xL = interp1(y([iL iL+1]), x([iL iL+1]), half, 'linear'); end
% right crossing
iRrel = find(y(im:end) <= half, 1, 'first');
if isempty(iRrel), xR = x(end);
else
    iR = im-1+iRrel;
    if iR<=1, xR = x(end);
    else, xR = interp1(y([iR-1 iR]), x([iR-1 iR]), half, 'linear'); end
end
w = abs(xR - xL);
end
