classdef SimpleAligner < sync.BaseAligner
%SIMPLEALIGNER handles alignment of two pre-matched pulse trains

    properties
        % ANCHOR pulse reference point used for alignment
        % 
        % Alignment can be defined in relation to the pulse rise times
        % "tRise", or fall times "tFall", or to the pulse mid-point, 
        % "tCenter" (default). "tCenter" may help to reduce alignment
        % errors if sampling quantization is a significant factor.
        anchor                  (1,1)   string  {mustBeMember(anchor, ["tRise", "tFall", "tCenter"])} = sync.BaseAligner.DEFAULT_ANCHOR;
    end
    
    methods
        function self = SimpleAligner(varargin)
            % SIMPLEALIGNER constructor
            self = self@sync.BaseAligner(varargin{:});
            if nargin
                ptA = varargin{1};
                ptB = varargin{2};
                n = ptA.nPulses;
                assert(ptB.nPulses == n, ...
                    "sync:SimpleAligner:SimpleAligner:unmatchedLength", ...
                    "The lengths of the two pulse trains must match.");
                self.iMatchAToB = (1:n)';
                self.iMatchBToA = (1:n)';
            end
        end
    end
    
    methods (Access = protected)
        function anchor = onGetAnchor(self)
            anchor = self.anchor;
        end
    end
    
end