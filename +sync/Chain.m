classdef Chain < sync.Base
    %CHAIN represents a chain of Aligner objects
    
    properties
        name            (1,1)   string = ""     % descriptive string
        alignerClass    (1,1)   string = "PatternAligner"
    end
    
    properties(SetAccess = protected, Hidden)
        aligners
    end
    
    properties (Dependent)
        pulseTrains
    end
    
    properties (Dependent, Hidden)
        nChain      % number of Aligner objects in chain
    end
    
    methods
        
        function self = Chain(arg)
            if nargin
                if isa(arg, 'sync.BaseAligner')
                    self.aligners = arg;
                    cls = strsplit(class(arg), '.');
                    self.alignerClass = cls{2};
                elseif isa(arg, 'sync.PulseTrain')
                    self = self.append(arg);
                end
                self.name = string(self);
            end
        end
        
        function self = append(self, pulseTrains, side)
            % APPEND add new PulseTrains to chain RHS
            if nargin < 3 || isempty(side), side = "right"; end
            assert( ...
                isa(pulseTrains, 'sync.PulseTrain'), ...
                'Chain:append:invalidPulseTrain', ...
                'Input must be a PulseTrain object');
            pulseTrains = pulseTrains(:);
            if self.nChain > 0
                if side == "right"
                    trainsAll = [self.pulseTrains(end); pulseTrains];
                elseif side == "left"
                    trainsAll = [pulseTrains; self.pulseTrains(1)];
                end
            else
                trainsAll = pulseTrains;
                if numel(trainsAll) < 2
                    error('A Chain must contain two or more PulseTrain objects');
                end
            end
            
            cls = "sync." + self.alignerClass;
            newAligners = feval(cls + ".empty");
            for n = 1:numel(trainsAll)-1
                newAligners(n, 1) = feval(cls, trainsAll(n), trainsAll(n+1));
            end
            
            if side == "right"
                self.aligners = [self.aligners; newAligners];
            elseif side == "left"
                self.aligners = [newAligners; self.aligners];
            end
        end
        
        function tR = mapForward(self, tL, varargin)
            %ATOB convert times from the first to the last timeframe
            for n = 1:self.nChain-1
                aligner = self.aligners(n);
                tL = aligner.aToB(tL, varargin{:});
            end
            tR = tL;
        end
        
        function tL = mapReverse(self, tR, varargin)
            %BTOA convert times from the last to the first timeframe
            for n = 1:self.nChain-1
                aligner = self.aligners(n);
                tR = aligner.bToA(tR, varargin{:});
            end
            tL = tR;
        end
        
        function plot(self, varargin)
            inp = inputParser();
            inp.addParameter('axes', gca());
            inp.addParameter('color', 'k');
            inp.addParameter('yPos', 0);
            inp.addParameter('lineWidth', 0.5);
            inp.parse(varargin{:});
            P = inp.Results;
            
            ax = P.axes;
            hold(ax, 'on');
            
            trains = self.pulseTrains;
            xOffset0 = -trains(1).tRise(1);
            yPos0 = P.yPos;
            
            for tr = 1:self.nChain
                train = trains(tr);
                yPos = yPos0 + [0, 1];
                if tr > 1
                    aligner = self.aligners(tr-1);
                    xOffsetTr = aligner.aToB(0);
                else
                    xOffsetTr = 0;
                end
                xOffset = xOffset0 - xOffsetTr;
                train.plot( ...
                    'axes', ax, ...
                    'color', P.color, ...
                    'yPos', yPos, ...
                    'xOffset', xOffset);
                
                if tr > 1
                    tA = aligner.tAMatched + xOffset0;
                    tB = aligner.tBMatched + xOffset;

                    sz = size(tA);
                    x = [tA, tB, nan(sz)]';
                    y = [(yPos0-1)*ones(sz), yPos(1)*ones(sz), nan(sz)]';
                    line(ax, x(:), y(:), ...
                        'color', 'r', ...
                        'lineStyle', '-', ...
                        'lineWidth', P.lineWidth);
                end
                
                yPos0 = yPos(1) + 2;
                xOffset0 = xOffset;
            
            end
            
            train = trains(1);
            xOffset = -train.tRise(1);
            t = train.tRise;
            x = [t, t, nan(size(t))]' + xOffset;
            y = [P.yPos-1; yPos(2)+1; nan] + zeros(size(x));
            line(ax, x(:), y(:), 'color', 'k', 'lineStyle', ':');
            yt = 0.5 + 2*((1:self.nChain) - 1);
            yticks(yt);
            yticklabels([trains.name]);
            
            ax.YLim = [P.yPos, yPos(2)] + [-1, 1];
            ax.XLim = [0, 100];
%             ax.YTick = [];
            fig = ax.Parent;
            xlabel(ax, sprintf('Time (%s)', trains(1).name));
            zoom(fig, 'xon');
            pan(fig, 'xon');
            
        end
        
        function str = string(self)
            name = self.pulseTrains(1).name;
            for n = 2:self.nChain
                name = name + sprintf(" -> %s",  self.pulseTrains(n).name);
            end
            str = name;
        end
        
        function val = get.pulseTrains(self)
            aligners = self.aligners;
            if isempty(aligners)
                val = sync.PulseTrain.empty();
            else
                val = aligners(1).pulseTrainA;
                for n = 1:numel(aligners)
                    val = [val; aligners(n).pulseTrainB];
                end
            end
        end
        
        function val = get.nChain(self)
            val = numel(self.pulseTrains);
        end       
        
        function result = gt(self, arg)
            assert(isa(arg, 'sync.PulseTrain'))
            result = self;
            self = self.append(arg, "right");
        end
        
    end
    
end
