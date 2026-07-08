%% Script to obtain calcium imaging peaks following cell extraction and manual sorting

% Used in combination with calciumImagingAnalysis pipeline (https://bahanonu.github.io/calciumImagingAnalysis/)

% Lucie Descamps and Miguel Carvalho 2021 

%% Data Upload
% upload outputs of signal extraction algorithm and manual cell sorting.
[file,path] = uigetfile('.mat', 'Select cell extraction and manual sorting outputs', 'Multiselect', 'on', path);
for ffile = 1:2
    load(char(fullfile(path, file(ffile))));
end

%% Select valid manually sorted cells.
if exist('validCNMFE', 'var')
    cell_idx = find(validCNMFE == 1);
    goodTraces= cnmfeAnalysisOutput.extractedSignals(cell_idx, :);
    goodfilters =cnmfeAnalysisOutput.extractedImages(:, :, cell_idx);
elseif exist('validEXTRACT', 'var')
    cell_idx = find(validEXTRACT == 1);
    goodTraces = extractAnalysisOutput.traces(cell_idx, :);
    goodfilters =extractAnalysisOutput.filters(:, :, cell_idx);
end


%% Filter the traces to only keep the significant transients.

cellTraces = goodTraces;
numCells = size(cellTraces,1);
cellTraces_filtered_3std = zeros(numCells,length(cellTraces));
mask_std = cellTraces > 3.*std(cellTraces, [], 2);
cellTraces_filtered_3std(mask_std) = cellTraces(mask_std);
%Finds the peaks from the filtered traces
[signalPeaks, signalPeaksArray, ~] = computeSignalPeaks(cellTraces_filtered_3std);
%Normalizes the filtered traces with zscore - Could be changed to a matfun
%line!
for i = 1:numCells
    cellTraces_filtered_3std_zscored(i,:) = zscore(cellTraces_filtered_3std(i,:));
end

%% Save final results.
processedCellData = struct();

processedCellData.cellIdx = cell_idx;
processedCellData.cellTracesZscored = cellTraces_filtered_3std_zscored;
processedCellData.cellTraces = cellTraces_filtered_3std;
%save non zscored one!
processedCellData.cellPeaks = signalPeaks;
processedCellData.cellFootPrints = goodfilters;

save(fullfile(path, 'processedCellData_filteredTraces.mat'), 'processedCellData');