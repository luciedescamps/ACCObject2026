classdef PulseTrain < sync.Base
    %PULSETRAIN describes a recording of synchronization pulses
    
    properties(SetAccess = protected)
        tRise       (:,1)   double  % pulse onset times
        tFall       (:,1)   double  % pulse offset times
    end
    
    properties
        name        (1,1)   string
    end
    
    properties (Dependent)
        tCenter
        nPulses
        pulseLengths
        meanIpi
        maxIpi
        cvIpi
        medianPulseLength
        tLims
    end
    
    methods
        function self = PulseTrain(tRise, tFall, name)
            %PULSETRAIN constructor
            %
            %PT = PULSETRAIN(TON, TOFF) creates a PulseTrain object from
            %the onset and onset times specified by TON and TOFF.
            if nargin
                checkTimes(tRise, tFall)
                self.tRise = tRise;
                self.tFall = tFall;
                if nargin == 3
                    self.name = name;
                end
            end
        end
        
        function [self, inds] = subsample(self, subsampleMode, param)
            smode = lower(subsampleMode);
            if smode == "proportion"
                n = round(self.nPulses * param);
                inds = randsample((1:self.nPulses)', n);
            elseif smode == "number"
                n = param;
                inds = randsample((1:self.nPulses)', n);
            elseif smode == "inds"
                inds = param;
            else
                error('Invalid subsampling mode "%s".', subsampleMode);
            end
            inds = sort(inds);
            self.tRise = self.tRise(inds);
            self.tFall = self.tFall(inds);
        end
        
        function [inds, times] = findFlankingPulses(self, anchor, t)
            tPulse = self.(anchor);
            v = t > tPulse(1) & t < tPulse(end);
            iBef = rg.helpers.binarySearch(t, tPulse, -1);
            iAft = rg.helpers.binarySearch(t, tPulse, 1);
            iBef(~v) = NaN;
            iAft(~v) = NaN;
            inds = [iBef, iAft];
            times = nan(size(inds));
            times(v, :) = tPulse(inds(v, :));
        end
        
        function plot(self, varargin)
            inp = inputParser();
            inp.addParameter('axes', gca());
            inp.addParameter('color', 'k');
            inp.addParameter('yPos', [0, 1]);
            inp.addParameter('yLim', [-1, 2]);
            inp.addParameter('lineWidth', 0.5);
            inp.addParameter('xOffset', 0);
            inp.parse(varargin{:});
            P = inp.Results;
            
            xOff = P.xOffset;
            xR = repmat(self.tRise, 1, 2);
            yR = zeros(size(xR)) + P.yPos;
            
            xF = repmat(self.tFall, 1, 2);
            yF = zeros(size(xF)) + P.yPos([2, 1]);
            
            x = [xR, xF]' + xOff;
            y = [yR, yF]' ;
            
            line(P.axes, x(:), y(:), 'color', P.color, 'lineWidth', P.lineWidth);
            P.axes.YLim = P.yLim;
        end
        
        function val = get.tCenter(self)
            val = mean([self.tRise, self.tFall], 2);
        end
        
        function val = get.nPulses(self)
            val = numel(self.tRise);
        end
        
        function val = get.pulseLengths(self)
            val = self.tFall-self.tRise;
        end
        
        function val = get.meanIpi(self)
            if self.nPulses > 1
                val = mean(self.tFall(2:end) - self.tRise(1:end));
            else
                val = NaN;
            end
        end
        
        function val = get.medianPulseLength(self)
            val = median(self.pulseLengths);
        end
        
        function val = get.cvIpi(self)
            ipi = self.tFall(2:end) - self.tRise(1:end-1);
            val = std(ipi) ./ mean(ipi);
        end
        
        function val = get.maxIpi(self)
            ipi = mean(self.tFall(2:end) - self.tRise(1:end-1));
            val = max(ipi);
        end
        
        function val = get.tLims(self)
            if self.nPulses == 0
                val = [];
            else
                val(1) = self.tRise(1);
                val(2) = self.tFall(end);
            end
        end
        
        function chain = gt(self, arg)
            if isa(arg, 'sync.PulseTrain')
                chain = sync.Chain([self; arg(:)]);
            elseif isa(arg, 'sync.Chain')
                chain = arg;
                chain.append(self, "left");
            end
        end
        
    end
    
    methods (Static)
        
        function [pulseTrains, info] = simulate(varargin)
            inp = inputParser();
            inp.addParameter('nTrains', 2);
            inp.addParameter('nPulses', 1000);
            inp.addParameter('intervalRange', [0.3, 1.5]);
            inp.addParameter('sigmaDetect', 0.01);
            inp.addParameter('sigmaDrift', 0);
            inp.addParameter('duration', 0.05);
            inp.addParameter('sigmaDuration', 0.0001);
            inp.addParameter('detectRate', 1);
            inp.addParameter('offsets', 0);
            inp.addParameter('jumpPos', 0);
            inp.addParameter('jumpIncrement', 0);
            inp.parse(varargin{:});
            P = inp.Results;
            P = checkScalarParams(P, {'detectRate', 'sigmaDrift', 'sigmaDetect', 'offsets', 'duration', 'sigmaDuration', 'jumpPos', 'jumpIncrement'});
            
            
            % Calculate the "raw" train of pulses, without any detection or drift noise
            rng = P.intervalRange;
            x = (rand(P.nPulses, 1) .* diff(rng)) + rng(1);
            
            tRaw = cumsum(x);
            info.tRaw = tRaw;
            info.tRawMax = max(tRaw);
            info.P = P;
            
            MIN_PULSE_DURATION = 1e-6;
            
            for n = 1:P.nTrains
                
                % Generate simulated detection of the pulse train:
                nDetect = ceil(P.detectRate(n)*P.nPulses);
                
                % 1) Add random-walk drift noise
                errDrift = cumsum(randn(nDetect, 1)*P.sigmaDrift(n));
                
                % 2) Add pulse detection time noise
                errDetect = randn(nDetect, 1)*P.sigmaDetect(n);
                
                % 3) Add pulse duration noise
                errDuration = randn(nDetect, 1)*P.sigmaDuration(n);
                
                % 4) Add discontinuties
                jPos = P.jumpPos(n);
                errJump = zeros(nDetect, 1);
                if jPos ~= 0
                    for j = 1:numel(jPos)
                        inds = jPos : P.nPulses;
                        errJump(inds) = errJump(inds) + P.jumpIncrement(n);
                    end
                end
                
                % 5) Randomly subsample at specified detection rate
                indsDetect = sort(randperm(P.nPulses, nDetect))';
                
                tRise = tRaw(indsDetect) + P.offsets(n) + errDrift + errJump + errDetect;
                tFall = tRise + max(MIN_PULSE_DURATION, P.duration(n) + errDuration); % don't allow negative durations
                name = "PulseTrain " + char('A'+n-1);
                pulseTrains(n) = sync.PulseTrain(tRise, tFall, name);
                
                info.errDrift{n, 1} = errDrift;
                info.errDetect{n, 1} = errDetect;
                info.tRawInds{n, 1} = indsDetect;
            end
        end
        
        function pulseTrain = fromSignal(t, y)
            %FROMSIGNAL create a 
        end
        
    end
    
end

function P = checkScalarParams(P, params)
nParams = numel(params);
for n = 1:nParams
    param = params{n};
    val = P.(param);
    if isscalar(val)
        P.(param) = repmat(val, P.nTrains, 1);
    else
        assert( ...
            numel(val)==P.nTrains, ...
            'simulatePulseTrains:badParamSize', ...
            'Parameter "%s" should either be a scalar or have the same number of elements as the number of pulse trains.', ...
            param);
    end
end
end

function checkTimes(tRise, tFall)
tAll = [tRise(:), tFall(:)];
assert(isequal(size(tRise), size(tFall)), ...
    'PulseTrain:checkTimes:unequalSizes', ...
    'Sizes of "tRise" and "tFall" and must be equal');

    assert(issorted(tAll'), ...
        'PulseTrain:checkTimes:tRiseFallNotSorted', ...
        'Each value of "tFall" must be greater than the corresponding value of "tRise"');
    assert(issorted(tAll), ...
        'PulseTrain:checkTimes:tNotSorted', ...
        'Intput vectors "tRise" and "tFall" must be sorted ascending');
end