classdef BaseAligner < sync.Base & matlab.mixin.CustomDisplay
    %ALIGNER alignment tool for a random sync pulse train recorded in
    %two timeseries
    
    properties
        % Pulse trains to be aligned
        pulseTrainA sync.PulseTrain
        pulseTrainB sync.PulseTrain
        
        % Sync interpolation parameters
        extrapLimit double = []     % maximum time to extrapolate beyond aligned pulses
    end
    
    properties (SetAccess = protected)
        iMatchAToB          % indices of matched tB pulses for each tA pulse
        iMatchBToA          % indices of matched tA pulses for each tB pulse 
    end
    
    properties(Dependent)
        % Convenience dependent props
        tA
        tB
        nA
        nB
        name
        
        tAMatched           % values of tA matched to tB
        tBMatched           % values of tB matched to tA
        
        pulseTrainAMatched
        pulseTrainBMatched
    end
    
    properties (Constant, Hidden)
        DEFAULT_ANCHOR = "tCenter";
    end
    
    methods
        
        function self = BaseAligner(varargin)
            %ALIGNER constructor
            %
            % OBJ = ALIGNER(PTA, PTB) creates a new Aligner object
            % for the pair of PulseTrain objects PTA and PTB.
            if nargin
                self.pulseTrainA = varargin{1};
                self.pulseTrainB = varargin{2};
            end
        end
        
        function self = reverse(self)
            %REVERSE switch places of timeseries A and B
            [self.pulseTrainB, self.pulseTrainA] = deal(self.pulseTrainA, self.pulseTrainB);
        end
        
        function [tB, validInterp] = aToB(self, tA, varargin)
            [tB, validInterp] = self.interpSub("aToB", tA, varargin{:});
        end
        
        function [tA, validInterp] = bToA(self, tB, varargin)
            [tA, validInterp] = self.interpSub("bToA", tB, varargin{:});
        end
        
        function plot(self, varargin)
            chain = sync.Chain(self);
            chain.plot(varargin{:});
        end
        
        function val = get.tA(self)
            if isempty(self.pulseTrainA)
                val = [];
            else
                fd = self.onGetAnchor();
                val = self.pulseTrainA.(fd);
            end
        end
        
        function val = get.tB(self)
            if isempty(self.pulseTrainA)
                val = [];
            else
                fd = self.onGetAnchor();
                val = self.pulseTrainB.(fd);
            end
        end
        
        function val = get.name(self)
            val = sprintf("%s <-> %s", self.pulseTrainA.name, self.pulseTrainB.name);
        end
        
        function chain = plus(self, obj)
            if isa(obj, 'sync.BaseAligner')
                chain = sync.Chain(self, obj);
            elseif isa(obj, 'sync.Chain')
                chain = obj;
                chain.append(self);
            else
                error('operator "+" can only be used with arguments of type "Chain" or "Aligner"');
            end
        end
        
        function val = get.pulseTrainAMatched(self)
            iBA = self.iMatchBToA;
            v = ~isnan(iBA);
            inds = iBA(v);
            pt = self.pulseTrainA.subsample("inds", inds);
            pt.name = pt.name + " (matched)";
            val = pt;
        end
        
        function val = get.pulseTrainBMatched(self)
            iAB = self.iMatchAToB;
            v = ~isnan(iAB);
            inds = iAB(v);
            pt = self.pulseTrainB.subsample("inds", inds);
            pt.name = pt.name + " (matched)";
            val = pt;
        end
        
        function val = get.tAMatched(self)
            iBA = self.iMatchBToA;
            v = ~isnan(iBA);
            val = self.tA(iBA(v));
        end
        
        function val = get.tBMatched(self)
            iAB = self.iMatchAToB;
            v = ~isnan(iAB);
            val = self.tB(iAB(v));
        end
        
        function val = get.nA(self)
            val = numel(self.tA);
        end
        
        function val = get.nB(self)
            val = numel(self.tB);
        end
    end
    
    methods (Hidden)
        function [t2I, validInterp] = interpSub(self, flag, t1I, method)
            if nargin < 4 || method == "", method = "brokenstick"; end
            
            GRAD_TOL = 0.001;
            method = lower(method);
            
            % Handle input times as a vector, reshape later.
            szI = size(t1I);
            t1I = t1I(:);
            
            % Create pulse trains without unmatched pulses
            train1 = self.pulseTrainAMatched;
            train2 = self.pulseTrainBMatched;
            
            if flag == "bToA"
                [train2, train1] = deal(train1, train2);
            elseif flag ~= "aToB"
                error('Invalid flag value "%s". Must be "bToA" or "aToB"', flag);
            end
            
            anchor = self.onGetAnchor;
            t1 = train1.(anchor);
            t2 = train2.(anchor);
            
            % If we have fewer than two points, interpolation is
            % impossible, so return NaN.
            if numel(t1) < 2
                t2I = nan(szI);
                validInterp = false(szI);
                return;
            end
            
            % Sort query times so we can easily check for problems
            doSort = ~issorted(t1I);
            if doSort
                [t1I, isort] = sort(t1I);
            end
                
            if method == "piecewise"
                t2I = interp1(t1, t2, t1I, 'linear', 'extrap');
                
            elseif any(method == ["regression", "brokenstick"])
                if method == "brokenstick"
                    [t2I, resid, mdls] = sync.brokenStickInterp(t1, t2, t1I);
                    nKinks = numel(mdls) - 1;
                    b = cellfun(@(mdl) mdl.Coefficients(2, 1).Variables, mdls);
                    if nKinks
                        fprintf('Found %u kinks in t1 vs. t2 relationship.\n', nKinks);
                    end
                elseif method == "regression"
                    % Fit robust regression
                    [b, c, resid] = robustReg(t1, t2);
                    t2I = b*t1I + c;
                end
                
                mse = mean(resid.^2);
                tolMse = (1e-5 * diff(t2([1, end]))) .^ 2;
                tolB = 1e-4; % assume clock doesn't drift by more than 1 second in 10,000
                if mse > tolMse
                    warning('Aligner:interpSub:badResid', ...
                        'Mean residual of regression fit (%.6g) is above maximum tolerance (%.6g).', ...
                        mse, tolMse);
                end
                if any(abs(b-1) > tolB)
                    warning('Aligner:interpSub:badGradient', ...
                        'Regression gradient (%.6g) exceeds tolerance (%.6g).', ...
                        max(abs(b)), tolB);
                end
            else
                error('Interpolation mode must be one of "spline", "piecewise" or "regression"');
            end
            
            % Check for sudden discontinuities in gradient
            m = diff(t2I)./diff(t1I);
            absdm = abs(m-nanmedian(m));
            thresh = max(GRAD_TOL,  5*median(absdm));
            badGrad = absdm > thresh;
            i1 = find(badGrad, 1, 'first');
            i2 = find(badGrad, 1, 'last');
            validInterp = true(size(t1I));
            validInterp(i1:i2+1) = false;
            if any(badGrad)
                warning('Aligner:interpSub:unstableGrad', 'Interpolated times have unstable gradient.');
            end
            tExtrap = self.extrapLimit;
            if ~isempty(tExtrap)
                validExtrap = t1I >= t1(1)-tExtrap & t1I <= t1(end)+tExtrap;
                validInterp = validInterp & validExtrap;
            end
            %t2I(~validInterp) = NaN;
            
            % Reverse the sorting
            if doSort
                [~, isortrev] = ismember(1:numel(t1I), isort);
                t2I = t2I(isortrev);
            end
            t2I = reshape(t2I, szI);
        end
        
    end
    
    methods (Access = protected)
        function grps = getPropertyGroups(self)
            import matlab.mixin.util.PropertyGroup
            grps(1) = PropertyGroup({'pulseTrainA', 'pulseTrainB', 'tA', 'tB', 'nA', 'nB', 'name'}, "PULSE TRAINS");
            grps(2) = PropertyGroup({'tAMatched', 'tBMatched', 'iMatchAToB', 'iMatchBToA'}, "ALIGNED DATA");
            grps(3) = PropertyGroup({'extrapLimit'}, "INTERPOLATION PARAMETERS");
        end
    end
    
    methods (Abstract, Access = protected)
        anchor = onGetAnchor(self)
    end
    
end

function [b, c, resid, mdl] = robustReg(x, y)
mdl = fitlm(x, y, 'RobustOpts', 'on');
coeffs = mdl.Coefficients(:, 1).Variables;
b = coeffs(2);
c = coeffs(1);
resid = mdl.Residuals(:, 1).Variables;
end