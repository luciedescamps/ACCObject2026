function [yi, resid, mdls] = brokenStickInterp(x, y, xi)
% Model y with piecewise linear spline, so that we can find the changes in
% gradient

iKinks = sync.findBreaks(x, y);

nKinks = numel(iKinks);

yi = nan(size(xi));

resid = nan(size(x));

iBefKink0 = 0;
% iBefKinkI0 = 0;
iAftKinkI0 = 1;

for k = 1:nKinks+1
    
    if k <= nKinks
        % Current region ends with a kink: stop at the last sample before
        % the kink
        iBefKink = iKinks(k) - 1;
        xBefKink = x(iBefKink);
    else
        % Last region: handle all remaining data
        iBefKink = numel(x);
        xBefKink = inf;
    end
    inds = iBefKink0+1 : iBefKink;
    
    mdl = fitlm(x(inds), y(inds), 'RobustOpts', 'on');
    mdls{k} = mdl;
    warning('on', 'stats:LinearModel:RankDefDesignMat');
    resid(inds) = mdl.Residuals(:, 1).Variables;
    coeffs = mdl.Coefficients(:, 1).Variables;
    b = coeffs(2);
    c = coeffs(1);
    
    iBefKinkI = find(xi <= xBefKink, 1, 'last');
%     indsI = iBefKinkI0+1 : iKinkI;
    indsI = iAftKinkI0 : iBefKinkI;
    yi(indsI) = b*xi(indsI) + c;
    iBefKink0 = iBefKink;
%     iBefKinkI0 = iBefKinkI;
    
    if k <= nKinks
        iAftKinkI0 = find(xi >= x(iBefKink+1), 1, 'first');
    end
end

assert(~any(isnan(resid)));
% assert(~any(isnan(yi)));

end