
% This script uses: 
%
% 1 - Tracking data extracted from DeepLabCut
% 2 - Single cell Ca2+ events extracted with calciumImagingAnalysis (https://github.com/bahanonu/calciumImagingAnalysis)
%
% Analyses Performed:
%
% 1 - Compute activity maps (rate maps)
% 2 - Identification of object-coding cells
%
%It needs the BNT toolbox (https://github.com/brkanter/BNT)
%
% Written by Lucie Descamps and Miguel Carvalho 2020

%% Identify directory for subject and session and define various experimental and analysis settings.

% get work directory from user. Data should be stored as folder(Animal)/folder(Session).
path = uigetdir (''); % change to user data directory. 
split_path = split(path, "\");
session = split_path(end);
subject = split_path(end-1);

% generate output folder.
outputFolder = fullfile(path, 'Results');
mkdir(outputFolder);

% define user analysis settings.
global hippoGlobe
fNames = fieldnames(hippoGlobe);
for iField = 1:length(fNames)
   eval(sprintf('%s = hippoGlobe.%s;',fNames{iField},fNames{iField}));
end


%% Tracking Data %%
% Upload tracking data from DeepLabCut and perform some cleaning.
[syncMethod] = helper.identifySyncMethod;
disp('Processing Tracking Data');
[trackingData_raw]= upload.importDLC(path, param); 

% Rescale tracking data.
[trackingData, originalTrackingLimits] = helper.rescalePath(path, param, trackingData_raw);
disp('Done!');

%User inputs which sync method to use. If using nvista LED, identifies when the nVista LED is first on. If using
% the random LED blinking, will run the Aligner script.

if syncMethod == 1
    [first_frame] = helper.identifyLED(path);
    trackingData_Masked = structfun(@(x) x(first_frame:end,:), trackingData,  'UniformOutput', false);
    trackingData = trackingData_Masked;
    disp('Done!');
end

if syncMethod == 2 && ~exist ('aligner')
    [aligner, pulseTrains, remove] = sync.createVectors;
end

if syncMethod == 2
    for i = 1:length(trackingData.(param.refBdPt))
        behaviour_ts(i,1) = i/60; %Because tracking data is sampled 60Hz
    end
    trackingData = structfun(@(x) horzcat(behaviour_ts(:, 1), x), trackingData, 'UniformOutput', false);
    posData = trackingData.(param.refBdPt);
    disp('Done');
end

 %% Object Locations %%
% Upload tracking video and identify object coordinates
[objCoordinates, trackingVideo, ~, objtype] = upload.obj_Locations(path);
numObj = size(objCoordinates, 2);
disp('Done');

%% Calcium Event Data %%
% Upload Ca2+ event file from Ca2+ imaging data.
% If needed, this step crops Ca2+ data and downsamples behavior tracking data to match Ca2+ sampling rate.
% See ca2_Import documentation for details
disp('Processing Calcium Event Data.');
 if syncMethod == 1 
        posData = [];
        aligner = [];
 end

[cellTraces, cellTraces_filtered_3std, signalPeaks, posData, match_Im_to_Bh] = upload.ca2_Import_Amplitude(path, posData, trackingData, syncMethod, aligner);
      
numCells = size(cellTraces,1);

% Save all cell tracking data and Ca2+ event information as .mat file in the output folder.
save(fullfile(outputFolder, 'posData.mat'), 'posData');
save(fullfile(outputFolder, 'cellData_filteredTraces.mat'), 'cellTraces');
disp('Done!');

%% RATE MAPS %%
cellAmplitudeMapsSmoothed = cell(1,numCells);
for iCell = 1: numCells
[map, edges] = calculate.rateMap(posData, cellTraces(iCell, :)', 'caFeature', 'amplitude', 'smooth', smoothing, 'binWidth', binWidth, 'minTime',minTime, 'limits', mapLimits); % generate rate map.
cellAmplitudeMapsSmoothed{iCell} = map.z;
cellMaps{iCell} = map;
end

%% Object scoring

xMin = originalTrackingLimits(1);
xMax = originalTrackingLimits(2);
yMin = originalTrackingLimits(3);
yMax = originalTrackingLimits(4);


xObjPos(:, 1) = abs(objCoordinates(:, 1) - xMax);
yObjPos(:, 1) = abs(objCoordinates(:, 2)- yMax);
xRange = xMax - xMin;
yRange = yMax - yMin;
xObjRescale = param.arenaDim/2 - (xObjPos*param.arenaDim)/xRange;
yObjRescale = param.arenaDim/2 - (yObjPos*param.arenaDim)/yRange;
objCoordinates = horzcat(xObjRescale, yObjRescale);

limitsX = [mapLimits(1) mapLimits(2)];
limitsY = [mapLimits(3) mapLimits(4)];
nBinsX = ceil((limitsX(2) - limitsX(1)) / binWidth);
nBinsY = ceil((limitsY(2) - limitsY(1)) / binWidth);

xEdges = limitsX(1):binWidth:limitsX(1) + binWidth*nBinsX; % create nBinsX bins
yEdges = limitsY(1):binWidth:limitsY(1) + binWidth*nBinsY; % create nBinsY bins

[~, ~, objXBin] = histcounts(round(objCoordinates(:,1)), xEdges);
[~, ~, objYBin] = histcounts(round(objCoordinates(:,2)), yEdges);
objectBins = (objectBins-1)/2; % we demark the object zone as an even number of bins left and right of the bin that has the centre of the object

objMask = zeros(nBinsX, nBinsY);


for obj = 1 : numObj
    xx = objXBin(obj, :)-objectBins:1: objXBin(obj, :) + objectBins;
    yy = objYBin(obj, :)-objectBins:1: objYBin(obj, :) + objectBins;
    if any(xx<1) || any(yy<1)
        outBinsIndex = [find(xx >=1, 1) find(yy >=1, 1)];
        xx = xx(1, outBinsIndex(1):end);
        yy = yy(1, outBinsIndex(2):end);
    elseif any(xx>nBinsX) || any(yy>nBinsY)
        outBinsIndex = [find(xx <=nBinsX, 1, 'last') find(yy <=nBinsX, 1, 'last')];
        xx = xx(1, 1:outBinsIndex(1));
        yy = yy(1, 1:outBinsIndex(2));
    end
    objMask(yy, xx) = obj;
end

% object scoring based on shuffling of Ca2+ events along animal's path
threshold = 95; %percentile threshold
 [objectRateScoring, allShuffledValues] = score.objScoringShuff_circshift(posData, cellTraces_filtered_3std, objMask, threshold);

objScore1 = cell2mat(objectRateScoring(2:end, 4));
objScore2 = cell2mat(objectRateScoring(2:end, 8));
objScore = [objScore1 objScore2];

% object temporal occupancy analyses
% temporal occupancy of each object zone.
objectTimeScoring = {};
objectTimeComp = {};
occupancyMap = map.time;
% identify time spent in each object zone.
for iobj = 1: numObj
    objTimeBins = occupancyMap(objMask == iobj);
    objectTimeComp{1, iobj} = objTimeBins;
    objMeanTime = nanmean(objTimeBins);
    objectTimeScoring = [objectTimeScoring objTimeBins objMeanTime];
end

% run simple stats.
skipStats = find(cellfun(@(x) numel(x)==1 & any(isnan(x)), objectTimeScoring));
if length(objectTimeComp) < 3 % stats for only 2 objects
    %If the animal hasn't visited the object bins, then ObjectTimeComp contains only nans, so we should skip the stat calculation.
    if ~isempty(skipStats)
        disp('At least one object has not been explored, no statistical comparison of object occupancy will be computed.');
        objTimeP = NaN;
        objTimeStats = NaN;
        objectTimeScoring = [objectTimeScoring objTimeP objTimeStats];
    else
        [objTimeP, ~, objTimeStats] = ranksum(objectTimeComp{1, 1}, objectTimeComp{1, 2});
        objectTimeScoring = [objectTimeScoring objTimeP objTimeStats];
    end
    
    timeLabels = {};
    for lb = 1: numObj
        objLabels = {sprintf('Obj %d Bins Time', lb), sprintf('Obj %d Mean Time', lb)};
        timeLabels = [timeLabels objLabels];
    end
    statsLabels = {sprintf('Wilcoxon p-value'), sprintf('Wilcoxon stats')};
    timeLabels = [timeLabels statsLabels];
    
else % stats for more than 2 objects
    if ~isempty(skipStats)
        disp('At least one object has not been explored, no statistical comparison will be computed.');
        p = NaN;
        stats = NaN;
        multComp = NaN
        objectTimeScoring = [objectTimeScoring p multComp(:, 6)];
    else
        [p, ~, stats] = kruskalwallis(cell2mat(objectTimeComp), [], 'off');
        multComp = multcompare(stats, 'CType', 'bonferroni');
        objectTimeScoring = [objectTimeScoring p multComp(:, 6)];
    end
    
    timeLabels = {};
    for lb = 1: numObj
        objLabels = {sprintf('Obj %d Bins Time', lb), sprintf('Obj %d Mean Time', lb)};
        timeLabels = [timeLabels objLabels];
    end
    statsLabels = {sprintf('K-Wallis p-value'), sprintf('Mult Comp p-value')};
    timeLabels = [timeLabels statsLabels];
end
objectTimeScoring = vertcat(timeLabels, objectTimeScoring);
objTime = objectTimeScoring(2,find(cellfun(@(x) numel(x) >2, objectTimeScoring(2,:))));

% save shuffling results
save(fullfile(outputFolder, 'cellData_ShuffleResults.mat'), 'allShuffledValues');
writematrix(allShuffledValues, fullfile(outputFolder, 'cellData_ShuffleResults.xlsx'));

% save all cell object rate and scoring and session object time and scoring as .mat file in the output folder.
save(fullfile(outputFolder, 'cellData_ObjRateScoring.mat'), 'objectRateScoring');
temp= objectRateScoring(:,find(cellfun(@(x) numel(x)<2, objectRateScoring(2,:))));
xlswrite(fullfile(outputFolder, 'cellData_ObjRateScoring.xlsx'),temp);

save(fullfile(outputFolder, 'session_ObjTimeScoring.mat'), 'objectTimeScoring');
temp= objectTimeScoring(:,find(cellfun(@(x) isnumeric(x) && numel(x)<2, objectTimeScoring(2,:))));
writecell(temp, fullfile(outputFolder, 'session_ObjTimeScoring.xlsx'));

disp('Object Coding Analysis: Done!');