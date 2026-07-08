classdef Network < sync.Base
    %NETWORK represents a graph of synchronized recording sytems
    
    % TODO extrapLimit value is lost when resetting Aligners - should fix
    properties
        name            (1,1)   string = ""
        alignerClass    (1,1)   string = "PatternAligner"
    end
    
    properties (SetAccess = protected)
        pulseTrains (:,1) sync.PulseTrain
    end
    
    properties (SetAccess = protected, Hidden)
        graphObj
        aligners = struct();
        masterIndex
        alignmentParamsImpl (1,1) sync.AlignmentParams
    end
    
    properties (Dependent)
        extrapLimit
        masterName
    end
    
    properties (Dependent, SetAccess = protected)
        nPulseTrains
        nEdges
        alignmentParams 
    end
        
    
    methods
        
        function self = Network(pulseTrains, alignmentParams)
            if nargin
                self = self.setPulseTrainsSub(pulseTrains, false);
                if nargin == 2
                    self.alignmentParamsImpl = alignmentParams;
                end
                self = self.align();
            end
        end
        
        function self = setPulseTrains(self, pulseTrains)
            % SETPULSETRAINS set PulseTrains objects
            self = self.setPulseTrainsSub(pulseTrains, true);
        end
        
        function self = appendPulseTrains(self, pulseTrains)
            % APPENDPULSETRAIN append one or more PulseTrain objects
           pts = [self.pulseTrains; pulseTrains(:)];
           self = self.setPulseTrainsSub(pts, false);
        end
        
        function t2 = mapTimes(self, node1, node2, t1, varargin)
            % MAPTIMES map from one timebase to another
            if isequal(node1, node2)
                t2 = t1;
            else
                chain = self.getChain(node1, node2);
                t2 = chain.mapForward(t1, varargin{:});
            end
        end
        
        function t2 = mapTimesToMaster(self, node1, t1, varargin)
            % MAPTIMESTOMASTER map a timebase to the "master" timebase
            node1 = self.getNodes(node1);
            node2 = self.getNodes(self.masterName);
            t2 = self.mapTimes(node1, node2, t1, varargin{:});
        end
        
        function [pulseTrains, nodeInds] = getNodes(self, varargin)
            % GETNODES retrieve nodes (PulseTrain objs) of the graph
            % PT = NET.GETNODES(STR) retrieves the PulseTrain object
            % PT specified by its name STR. STR may be of type char or a
            % cell array of char, or alternatively it may be a
            % string or array of strings.
            %
            % PT = NET.GETNODES(IDX) retrives a PulseTrain object by
            % its node index IDX. To retrieve multiple PulseTrain objects,
            % IDX may be an array of indices.
            %
            % PT = NET.GETNODES(PTQ) retrieves PT using the input PulseTrain
            % object PTQ. The name of PTQ is compared with the names of
            % PulseTrains stored in NET, and the matching PulseTrain is
            % returned.
            %
            % [PT, IDX] = NET.GETNODES(...) additionally returns the
            % node indices of the returned PulseTrains in vector IDX.
            
            haveNodes = false;
            
            % Determine the input type
            % 1: PulseTrain
            if nargin > 1
                if isa(varargin{1}, 'sync.PulseTrain')
                    pulseTrains = [varargin{:}];
                    names = [pulseTrains.name];
                % 2: numeric array of indices
                elseif isnumeric(varargin{1})
                    nodeInds = [varargin{:}];
                    haveNodes = true;
                % 3: char or string
                else
                    % Cell array of chars
                    if numel(varargin) > 1
                        arg = varargin;
                    % Char or string array
                    else
                        arg = varargin{1};
                    end
                    names = string(arg);
                end
            else
                nodeInds = (1:self.nPulseTrains)';
                haveNodes = true;
            end
            
            if ~haveNodes
                allNames = self.graphObj.Nodes.Variables;
                [nameExists, nodeInds] = ismember(names, allNames);
                if ~all(nameExists)
                    error('Network:getNode:nameNotFound', ...
                        'Could not find PulseTrain objects with requested name(s).');
                end
            end
            pulseTrains = self.pulseTrains(nodeInds);
        end
        
        function h = plot(self, ax, nodes)
            if nargin < 2 || isempty(ax), ax = gca; end
%                 n = self.nPulseTrains;
%                 nodeFontWeight = repmat({'normal'}, n, 1); % requires r2019a
                nodeLabel = [self.pulseTrains.name];
                nodeLabel = nodeLabel.cellstr();
                idx = self.masterIndex;
                if ~isempty(idx)
%                     nodeFontWeight{idx} = 'bold';
                    nodeLabel{idx} = ['(M) ', nodeLabel{idx}];
                end
            args = { ...
                'lineWidth', 2, ...
                'nodeLabel', nodeLabel};
            if nargin == 3
                [~, nodeInds] = self.getNodes(nodes);
                [iNode, iEdge] = self.shortestPath(nodeInds(1), nodeInds(2));
                lineWidths = ones(self.nEdges, 1);
                lineWidths(iEdge) = 3;
                args = [args, {'lineWidth', lineWidths}];
            end
            h = plot(ax, self.graphObj, args{:});
        end
        
        function val = get.nPulseTrains(self)
            val = numel(self.pulseTrains);
        end
        
        function val = get.nEdges(self)
            val = size(self.graphObj.Edges, 1);
        end
        
        function self = set.extrapLimit(self, val)
            assert(numel(val)<=1 && ~any(val<0), ...
                'sync:Network:setExtrapLimit:invalidSize', ...
                'Value of parameter "extrapLimit" must be empty or a nonnegative scalar.');
            [aligners, keys] = self.getAllAligners();
            for a = 1:numel(aligners)
                aligner = aligners{a};
                aligner.extrapLimit = val;
                self.aligners.(keys{a}) = aligner;
            end
        end
        
        function val = get.extrapLimit(self)
            if self.nEdges == 0
                val = NaN;
            else
                aligners = self.getAllAligners();
                val = aligners(1).extrapLimit;
            end
        end
        
        function val = get.masterName(self)
            idx = self.masterIndex;
            if isempty(idx)
                val = "";
            else
                val = self.pulseTrains(idx).name;
            end
        end
        
        function self = set.masterName(self, val)
            names = [self.pulseTrains.name];
            matches = names == val;
            assert(any(matches), ...
                'sync:Network:setMasterName', ...
                'Value of property "masterName" must specify a valid PulseTrain name.');
            self.masterIndex = find(matches);
        end
        
        function val = get.alignmentParams(self)
            val = self.alignmentParamsImpl;
        end
        
    end
    
    methods (Hidden)
        
        function self = reset(self)
            % Reset alignment
            self.aligners = struct();
            self = self.align();
        end
        
        function [self, modified] = align(self)
            %Generate all Aligners
            naligners = numel(fieldnames(self.aligners));
            pts = self.pulseTrains;
            npts = numel(pts);
            combs = combnk(1:npts, 2);
            ncombs = size(combs, 1);
            verbose = self.alignmentParams.verbose;
            if verbose
                fprintf('Calculating all %u alignments between %u PulseTrains...', ncombs, npts);
            end
            for c = 1:ncombs
                comb = combs(c, :);
                pt1 = pts(comb(1));
                pt2 = pts(comb(2));
                [~, self] = self.getAligner(pt1, pt2);
                if verbose
                    fprintf('%u ', c);
                end
            end
            if verbose
                fprintf('done.\n');
            end
            modified = numel(fieldnames(self.aligners)) > naligners;
        end

        function chain = getChain(self, ptStart, ptEnd)
            ptInds = self.shortestPath(ptStart, ptEnd);
            pts = self.pulseTrains(ptInds);
            nChain = numel(pts);
            for c = 1:nChain-1
                [aligners(c), self] = self.getAligner(pts(c), pts(c+1));
            end
            chain = sync.Chain(aligners);
        end
        
        function [iNode, iEdge] = shortestPath(self, ptStart, ptEnd)
            [~, inds] = getNodes(self, ptStart, ptEnd);
            assert( ...
               inds(1) ~= inds(2), ...
                'Network:getChain:invalidNames', ...
                'The input PulseTrain names must be different.');
            iStart = inds(1);
            iEnd = inds(2);
            [iNode, ~, iEdge] = self.graphObj.shortestpath(iStart, iEnd);
        end
        
        function [aligner, self] = getAligner(self, pt1, pt2)
           % GETALIGNER retrieve Aligner for a pair of PulseTrains
           % N.B. needs to return self because this method caches any newly
           % created aligners
           pt1 = self.getNodes(pt1);
           pt2 = self.getNodes(pt2);
           aligners = self.aligners;
           
           % Generate key string for Aligner, direction determined
           % alphabetically
           [names, isort] = sort([pt1.name, pt2.name]);
           doReverse = isort(1) == 2;
           pts = [pt1, pt2];
           pts = pts(isort);
           key = names(1) + "_" + names(2);
           key = char(key);
           % Get Aligner with sorted direction
           
           if isfield(self.aligners, key)
                aligner = self.aligners.(key);
           else
               cls = "sync." + self.alignerClass;
               p = self.alignmentParamsImpl;
               aligner = feval(cls, pts(1), pts(2), p);
               self.aligners.(key) = aligner;
           end
           
           % Restore requested direction
           if doReverse
               aligner = aligner.reverse();
           end
        end
        
        function [aligners, keys] = getAllAligners(self)
            aligners = struct2cell(self.aligners);
            keys = fieldnames(self.aligners);
        end
        
    end
    
    methods (Access = protected)
        
        function self = setPulseTrainsSub(self, pulseTrains, reset)
            nTrains = numel(pulseTrains);
            A = ones(nTrains) - eye(nTrains);
            names = [pulseTrains.name];
            assert( ...
                ~any(cellfun(@isempty, names)), ...
                'Network:setPulseTrainsSub:unnamedPulseTrain', ...
                'All input PulseTrain objects must have a valid name.');
            assert( ...
                numel(unique(names)) == numel(names), ...
                'Network:setPulseTrainsSub:namesNotUnique', ...
                'All input PulseTrain names must be unique.');
            self.graphObj = graph(A, cellstr(names));
            self.pulseTrains = pulseTrains;
            if reset
                self = self.reset();
            end
        end
        
        
    end
    
    methods (Static)
       
        function obj = loadobj(obj)
            if obj.version < 2
                % bugfix: for versions <2, aligner objs were stored in a
                % containers.Map object, which caused weird problems
                % because of its 'handle' behavior. This upgrade converts
                % the old storage container to a simple struct instead,
                % which will copy its contents whenever the Network object
                % is reassigned.
                if isstruct(obj)
                    s = obj;
                    aligners = s.aligners.values;
                    p = aligners{1}.alignmentParams;
                    obj = sync.Network(s.pulseTrains, p);
                    obj.name = s.name;
                    obj.masterIndex = s.masterIndex;
                else
                    obj = obj.reset();
                end
                
                warning("sync:Network:upgraded", "sync.Network object has been upgraded from v%u to v%u", obj.version, sync.version);
                obj.version = sync.version;
                obj.upgraded = true;
            end
        end
        
    end
    
end
