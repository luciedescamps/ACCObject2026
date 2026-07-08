%% Classify neuron object responses.

% This function classifies object responsiveness based on neuronal activity
% within recording arena. It is a different method than the one used in
% calculate.objScoringRand, and it?s based on the shuffling of spikes/Ca2+
% events. 

% USAGE
%   [objectRateScoring, allShuffledValues] = score.objScoringShuff(posData, cellSpikeT, cellRateMap, objMask)
%
%   posData                    matrix with timestamps and XY coordinates of a selected reference point - usually scope.
%   cellSpikeT                  cell array with spike timestamps for each cell.
%   cellRateMap              cell array containing each cell's rate map.
%   objMask                    matrix indicating object zones.
%
% OUTPUT
%
%   objectRateScoring     matrix with firing rate and average firing rate of each cell in each object zone, firing rate ratio between in each 
%                                     object zone and rest of arena, and statistics between objecteach  zone firing rate and rest of arena.
%
%   allShuffledValues       matrix with all shuffled values for all cells.
%
% SEE ALSO
%   analyses.map
%   scripts.shuffling
%
%   Written by @miguelvmc 2020

function [objectRateScoring, allShuffledValues] = objScoringShuff_circshift(posData, cellTraces_filtered_3std, objMask, threshold)

% define user analysis settings.
global hippoGlobe
smoothing = hippoGlobe.smoothing;
binWidth =  hippoGlobe.binWidth;
minTime = hippoGlobe.minTime;
mapLimits = hippoGlobe.mapLimits;
shuffOffset = hippoGlobe.shuffOffset;
nShuffles = hippoGlobe.nShuffles;

disp('Starting Object Coding Analysis. This is going to take a while...go stretch your legs and drink some water.');
%% object firing rate analyses
numObj= unique(objMask(objMask >= 1));
% define range of possible random shuffle start for the session.
shuffEndIndex = find(posData(:, 1) >= (posData(end, 1) - shuffOffset), 1);
shuffStartIndex = find(posData(:, 1) <= shuffOffset, 1, 'last');

allShuffledValues = [];
objectRateScoring = cell(size(cellTraces_filtered_3std, 1), length(numObj)*4);

% run through cells
numCells = size(cellTraces_filtered_3std, 1);
 objects = length(numObj);
parfor iCell = 1:numCells

    celltrace = cellTraces_filtered_3std(iCell,:);
    cellMap = calculate.rateMap(posData,cellTraces_filtered_3std(iCell, :)', 'caFeature', 'amplitude', 'smooth', smoothing, 'binWidth', binWidth, 'minTime',minTime, 'limits', mapLimits); % generate rate map.;
    cellMap = cellMap.z;
    
    % generate a random start of the Ca2+ event sequence for each cell shuffle
    shuffStartIndices = randi([shuffStartIndex shuffEndIndex], 1, nShuffles);
    
    %do a for loop and use circhift 500 times, use shuffstartindex as shuff
    %value
    %end variable is shuffledCellEvents
    shuffledCellEvents = zeros(length(celltrace),nShuffles);
    
     shuffObjRate = zeros(size(shuffledCellEvents, 2), objects );
    for i = 1:nShuffles
        shuffledCellEvents(:,i) = circshift(celltrace, shuffStartIndices(1,i))';
    end
    
    % perform the shuffling
   
    for s= 1: size(shuffledCellEvents, 2)
        shuffTrace = shuffledCellEvents (:, s);
        shuffMap = calculate.rateMap(posData, shuffTrace, 'caFeature', 'amplitude', 'smooth', smoothing, 'binWidth', binWidth, 'minTime',minTime, 'limits', mapLimits); % generate rate map.
        
        % determine for each shuffle the average firing rate of each object zone.
        % NOTE: make sure to select the intended rate map.
        for iobj = 1: objects
            shuffObjRate(s, iobj) = nanmean(shuffMap.z(objMask == numObj(iobj)));
        end
    end
    allShuffledValues = vertcat(allShuffledValues, shuffObjRate);
    
    % gather data for each cell firing rate in each object zone and compare
    % observed mean firing rate to shuffled mean firing rate.
    cellShuff = {};
    for iobj = 1: length(numObj)
        objRateBins = cellMap(objMask == numObj(iobj));
        objMeanRate = nanmean(objRateBins);
        objShuff = prctile(shuffObjRate(:, iobj), threshold);
        objPassed = double(objMeanRate > objShuff);
        cellShuff = [cellShuff objRateBins objMeanRate objShuff objPassed];
    end
    objectRateScoring(iCell, :) = cellShuff;
    fprintf('Object Coding Analysis: Processed cell %d out of %d\n', iCell, numCells);
end

labels = {};
for lb = 1: length(numObj)
    if threshold == 99
        objLabels = {sprintf('Obj %d Bins Rate', lb), sprintf('Obj %d Mean Rate', lb), sprintf('Obj %d P99', lb), sprintf('Obj %d Passed', lb)};
    elseif threshold ==95
        objLabels = {sprintf('Obj %d Bins Rate', lb), sprintf('Obj %d Mean Rate', lb), sprintf('Obj %d P95', lb), sprintf('Obj %d Passed', lb)};
    end
    labels = [labels objLabels];
end
objectRateScoring = vertcat(labels, objectRateScoring);
empties = cellfun('isempty', objectRateScoring);
objectRateScoring(empties) = {NaN};

% count numbers of object cells
for i =1:objects
objScore(:,i) = cell2mat(objectRateScoring(2:end, 4*(i)));
end

% objScore1 = cell2mat(objectRateScoring(2:end, 4));
% objScore2 = cell2mat(objectRateScoring(2:end, 8));
%objScore = [objScore1 objScore2];

objCounter = zeros(1, length(numObj)+1);
for iCell = 1:numCells
    objPass =  find(objScore(iCell, :) == 1);
    if isempty(objPass)
        continue
    elseif numel(objPass) == 1
        objCounter(1, objPass) = objCounter(1, objPass)+1;
    elseif numel(objPass) == 2
        objCounter(1, objPass) = objCounter(1, objPass)+1;
        elseif numel(objPass) == 3
        objCounter(1, objPass) = objCounter(1, objPass)+1;
        elseif numel(objPass) == 4
        objCounter(1, objPass) = objCounter(1, objPass)+1;
    end
end
objCounter = [numCells, objCounter];

objCounterLabels = {'# Cells'};
for lb = 1: length(numObj)
    objLabels = {sprintf('#Obj %d Cells', lb)};
    objCounterLabels = [objCounterLabels objLabels];
end
objCounterLabels{1, end+1} = ('# Mult Obj Cells');

for add = 1:length(objCounter)
    objectRateScoring(1, end+1) = objCounterLabels(add);
    objectRateScoring(2, end) = num2cell(objCounter(add));
end

end