function inds = findBreaks(x, y, threshPercentile, plt)
%FINDBREAKS detect discontinuties in an increasing linear function

assert( ...
    issorted(x), ...
    'sync:brokenStickInterp:unsortedX', ...
    'Input "x" must have ascending sorted values.');

if nargin < 3 || isempty(threshPercentile), threshPercentile = 99; end
if nargin < 4 || isempty(plt), plt = false; end

MIN_KINK_SEPARATION = 50;   % Minimum number of points between two kinks
MIN_THRESH = 1e-6;          % Minimum gradient threshold for kink detection
N_MOV_MEAN = 20;            % moving average for gradient calculation

npoints = numel(x);
mmy = movmean(y(:), N_MOV_MEAN);
mmx = movmean(x(:), N_MOV_MEAN);

nhalf = floor(N_MOV_MEAN/2);
ic = N_MOV_MEAN+1 : npoints-N_MOV_MEAN;
i1 = ic - nhalf - 1;
i2 = ic + nhalf;
m = nan(npoints, 1);

for n = 1:numel(ic)
    y1 = mmy(i1(n));
    y2 = mmy(i2(n));
    x1 = mmx(i1(n));
    x2 = mmx(i2(n));
    m(ic(n)) = (y2-y1) ./ (x2-x1);
end

absm  = abs(m);
mabsm = nanmedian(absm);
absdevPc = prctile(abs(m-nanmedian(m)), threshPercentile);
absdabsm = abs(absm - mabsm);
thresh = max(5*absdevPc, MIN_THRESH);

if plt
   figure();
   subplot(1, 2, 1);
   histogram(m);
   title('Piecewise gradients');
   xlabel('Gradient');
   ylabel('Count');
   
   ax = subplot(1, 2, 2);
   histogram(ax, absdabsm);
   title(ax, 'Absolute gradients');
   xlabel(ax, 'Abs. gradient');
   ylabel(ax, 'Count');
   hold(ax, 'on');
   xpc = mabsm + absdevPc;
   ax.YLim = ax.YLim;
   h(1) = plot(ax, [1, 1]*xpc, ax.YLim, 'k--');
   h(2) = plot(ax, [1, 1]*thresh, ax.YLim, 'r--');
   xmx = max(absdabsm);
   ymx = ax.YLim(1)+0.1*mean(ax.YLim);
   ymx2 = ax.YLim(1)+0.13*mean(ax.YLim);
   h(3) = plot(ax, xmx, ymx, 'kv', 'markerFaceColor', 'k');
   text(ax, xmx, ymx2, string(xmx), ...
       'verticalAlignment', 'bottom', ...
       'horizontalAlignment', 'center')
   nAbove = sum(absdabsm > thresh);
   strs = { ...
       sprintf('%.0fth percentile', threshPercentile), ...
       sprintf('Threshold (%u above)', nAbove), ...
       sprintf('Maximum value') };
   legend(h, strs);
end

if any(absdabsm > thresh)
    [~, xKinks] = findpeaks(absdabsm, x, ...
        'MinPeakDistance', MIN_KINK_SEPARATION, ...
        'MinPeakHeight', thresh);
else
    inds = [];
    return;
end

nKinks = numel(xKinks);

for k = 1:nKinks
    if k <= nKinks
        inds(k, 1) = find(x < xKinks(k), 1, 'last');
    else
        inds(k, 1) = numel(x);
    end
end

end

