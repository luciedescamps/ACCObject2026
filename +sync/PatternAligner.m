classdef PatternAligner < sync.BaseAligner
    %ALIGNER alignment tool for a random sync pulse train recorded in
    %two timeseries

    properties(Dependent)
        % Alignment results
        nChunks
        percentMatched
        medianError
        meanSquaredError
    end
    
    properties (SetAccess = protected)
        % Alignment parameters
        alignmentParams (1,1) sync.AlignmentParams
        
        % Alignment results
        indsChunkA
        indsChunkB
        chunkErrors
        chunkValid
        logMse
        separation
        chunkErrorGmm
    end
    
    methods
        
        function self = PatternAligner(varargin)
            %ALIGNER constructor
            %
            % OBJ = ALIGNER(PTA, PTB) creates a new Aligner object
            % for the pair of PulseTrain objects PTA and PTB.
            %
            % OBJ = ALIGNER(PTA, PTB, PARAMS) specifies the alignment
            % parameters via the argument PARAMS. PARAMS must be an
            % AlignmentParams object.
            
            self = self@sync.BaseAligner(varargin{:});
            if nargin
                % With inputs, perform alignment
                if nargin == 3
                    params = varargin{3};
                    self.alignmentParams = params;
                end
                self = self.align();
            end
        end
        
        function self = reverse(self)
            %REVERSE switch places of timeseries A and B
            self = self.reverse@sync.BaseAligner();
            [self.iMatchAToB, self.iMatchBToA] = deal(self.iMatchBToA, self.iMatchAToB);
            [self.indsChunkB, self.indsChunkA] = deal(self.indsChunkA, self.indsChunkB);
        end
        
        function val = get.medianError(self)
            tA = self.tAMatched;
            tB = self.tBMatched;
            if numel(tA) == numel(tB)
                [~, ~, resid] = regress(tA, [ones(size(tB)), tB]);
                val = median(abs(resid));
            else
                val = nan;
            end
        end
        
        function val = get.meanSquaredError(self)
            tA = self.tAMatched;
            tB = self.tBMatched;
            if numel(tA) == numel(tB)
                [~, ~, resid] = regress(tA, [ones(size(tB)), tB]);
                val = mean(resid.^2);
            else
                val = nan;
            end
        end
        
        function val = get.nChunks(self)
            val = numel(self.chunkValid);
        end
        
        function val = get.percentMatched(self)
            val(1) = numel(self.tAMatched) ./ numel(self.tA);
            val(2) = numel(self.tBMatched) ./ numel(self.tB);
            val = val * 100;
        end
        
    end
    
    methods (Access = protected)
        function anchor = onGetAnchor(self)
            anchor = self.alignmentParams.anchor;
        end
        
        function self = align(self)
            params = self.alignmentParams;
            [iAB, iBA, info] = sync.alignSyncTrains(self.tA, self.tB, params);
            
            self.chunkValid     = info.chunkValid;
            self.indsChunkA     = info.chunkIndsA;
            self.indsChunkB     = info.chunkIndsB;
            self.chunkErrors    = info.chunkErrors(:, 1);
            self.logMse         = info.logMse;
            self.chunkErrorGmm  = info.gmm;
            self.separation     = info.separation;
            
            self.iMatchAToB = iAB;
            self.iMatchBToA = iBA;
        end
        
        function grps = getPropertyGroups(self)
            grps = self.getPropertyGroups@sync.BaseAligner();
            grps(end+1) = matlab.mixin.util.PropertyGroup( ...
                {'nChunks', 'indsChunkA', 'indsChunkB', 'chunkErrors1', 'chunkErrors2', 'chunkValid', 'percentMatched', 'medianError'}, ...
                "PATTERN ALIGNMENT DATA");
        end
    end
    
    methods (Static)
        function obj = loadobj(obj)
            if isstruct(obj)
                s = obj;
                if isfield(s, 'aligned')
                    s = rmfield(s, 'aligned');
                end
                obj = sync.PatternAligner;
                fds = fieldnames(s);
                for f = 1:numel(fds)
                    fd = fds{f};
                    obj.(fd) = s.(fd);
                end
            end
        end
    end
    
    
end