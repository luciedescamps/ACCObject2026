classdef Aligner
    %ALIGNER stub class, serves conversion from Aligner to PatternAligner
    
    methods
        function self = Aligner(varargin)
            error("The Aligner class is not functional and only exists for legacy reasons. Please don't try to use it!");
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if isfield(s, 'aligned')
                s = rmfield(s, 'aligned');
            end
            obj = sync.PatternAligner.loadobj(s);
        end
    end
end

