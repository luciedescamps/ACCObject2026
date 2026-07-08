%% Upload Ca2+ event file 

% User identifies the file with Ca2+ events, extracted with EXTRACT through Ciatah


% NOTES: 
%           (1) Following intial experiemnts of Lucie behavior videos recorded with ANYmaze were being cropped, thus endeing up
%                having a shorter duration than the Ca2+ event file. Ca2+ event files need to be cropped accordingly.
%           (2) For behavior videos recorded with ANYmaze the sampling rate does not match
%                Ca2+ sampling rate. The script downsamples tracking data, with one of two methods, so that they match.
%           (3) If syncing with the random LED pattern, we will modify the
%                 timestamps of the behaviour and calcium vector so that
%                 they are aligned, and will create a new trackingData
%                 vector aligned with the calcium data

% USAGE
%   [cellTraces, trackingData] = upload.ca2_Import(path, trackingVideo, trackingData, syncMethod, aligner)
%   path                                 path to subject/session folder containing RAW Ca2+ data extracted with OASIS  
%   trackingVideo                  behavior video used to track body parts 
%   posData                           structure with tracking data the reference body part
%   syncMethod                    Identifies which technique was used: nVista sync led, random led pattern
%   aligner                             Patter aligner for behaviour and calcium vectors if using the random led sync method 
%
% OUTPUT
%   cellTraces                       matrix with fluorescent trace filtered to only keep significant events (>3std) and normalized (zscore) for each cell
%   posData                          new structure with tracking data in the event of the need to downsample behavior data or align with the cell traces

% Written by Lucie Descamps and Miguel Carvalho 2020

%% 

function [cellTraces, cellTraces_filtered_3std, signalPeaks, posData, match_Im_to_Bh] = ca2_Import_Amplitude_WIP(path, posData, trackingData, syncMethod, aligner)

answer = questdlg('Do you need to filter the traces, normalize them and extract the peaks?', ...
	'Peaks Fever', ...
	'Yes','No', 'Maybe');

% Handle response
switch answer
    case 'Yes'
        disp(['Processing traces now'])
        %Computes processedCellData
        helper.peaksFever;
        %Fetches the cell Traces and filter them
        cellTraces_filtered_3std_zscored = processedCellData.cellTracesZscored;
        cellTraces_filtered_3std = processedCellData.cellTraces;
        signalPeaks = processedCellData.cellPeaks;
        numCells = size( cellTraces_filtered_3std_zscored,1);
        
    case 'No' %Need to work on this option!
        disp([answer 'Okidoki, then select the Ca2+ event file'])
        [file, path] = uigetfile('.mat', 'Select Ca2+ event file', path);
        selectedFile = fullfile(path,file);
        caEvents = readmatrix(selectedFile);
        caEvents = caEvents';
end

%% processing of Ca2+ data
imagingFR = 20; % nVista FR is ALWAYS 20Hz
session_duration_ca = length(cellTraces_filtered_3std_zscored);

%Here we need to check and adapt the data to one of these 3 case:
% (1) Tracking data acquired with AnyMaze: it will be shorter than the
%       calcium data, so we need to chop off the end of the calcium data
% (2)  Sync with the nVista Sync port: Need to chop off the start of the
%        tracking data
% (3) Sync with random led pattern: we will use the aligner to convert the
%       timestamps of the imaging data into the tracking referential, then use
%       knn to match the timestamps between the 2, then finally create a new
%       tracking vector that only comprises those matched with the imaging. 

% (1) Checks if tracking data was acquired using AnyMaze
answer = questdlg('Are you using tracking data that was acquired with AnyMaze and analysed with DLC?', ...
	'AnyMaze videos', ...
	'Yes','No', 'Maybe');

% Handle response
switch answer
    case 'Yes'
        anymaze = 1;
    case 'No'
        anymaze = 0;
end

if anymaze == 1
    
    [file, path] = uigetfile('.mp4', 'Select the DLC labeled video file', path);
    selectedFile = fullfile(path,file);
    video = VideoReader(selectedFile);
    
    new_imaging_duration = round(video.Duration*20);
    cellTraces_filtered_3std_zscored =  cellTraces_filtered_3std_zscored(:, 1:new_imaging_duration);
    cellTraces_filtered_3std = cellTraces_filtered_3std(:, 1:new_imaging_duration);
    signalPeaks = signalPeaks(:, 1:new_imaging_duration);
end

% (2) Checks if syncing with nVista sync port
%I am actually not sure the code below is correct, so need to double check it. 
answer = questdlg('Are you syncing data with the LED connected to the SYNC port of the nVista DAQ?', ...
	'Synchronization method', ...
	'Yes','No', 'Maybe');

% Handle response
switch answer
    case 'Yes'
        syncled = 1;
    case 'No'
        syncled = 0;
end

if syncled == 1
    newtrackingVideo_duration = caEvents_duration * 60 %60 is frame rate of behavioural video
    posData = structfun(@(x) x(1:newtrackingVideo_duration,:), posData,  'UniformOutput', false); %??
end

% (3) Checks if syncing done with random led pattern
if syncMethod == 2
    newvector = zeros(session_duration_ca,1); %Pre assignment
    for i = 2:session_duration_ca
        newvector(i,1) = i/20;%20 because nVista recordings are at 20Hz
    end
    newvector = newvector'; %This is a vector of seconds corresponding to frames indices
    cellTraces_filtered_3std_zscored_ts = [newvector;cellTraces_filtered_3std_zscored]; %Concatenating the timestamps to the filtered traces
    cellTraces_filtered_3std_zscored_ts = cellTraces_filtered_3std_zscored_ts'; %Flipping the axis so it's frames x cells (nb: the first "cell" is actually the timestamps, so there are cells-1 cells in the recording)
    %We are using the aligner function to convert the second indices in imaging referential into behaviour referential
    cellTraces_filtered_3std_zscored_ts_aligned(:,1)= aligner.aToB(cellTraces_filtered_3std_zscored_ts(:,1), 'piecewise');
    %We can then use a knnsearch to find the corresponding timestamps between behaviour and imaging
    match_Im_to_Bh = knnsearch(posData(:, 1), cellTraces_filtered_3std_zscored_ts_aligned(:,1));
    match_Bh_to_Im = knnsearch(cellTraces_filtered_3std_zscored_ts_aligned(:,1), posData(:, 1));
    posData = posData(match_Im_to_Bh,:);
    cellTraces = cellTraces_filtered_3std_zscored;
end


%% downsample tracking data to match Ca2+ sampling rate
% Apply downsampling procedure to behavior data.
% Two possible methods to choose from.
%Needs to be adjusted to posData.

if  anymaze == 1%We only want this to happen for the anymaze tracking data
    cellTraces = cellTraces_filtered_3std_zscored;
    trackingData_cell = struct2cell(trackingData);
    nPos = size(trackingData_cell{1},1);
   
    % method A. Selects specific indexes based on a vector with
    % same length of Ca2+ events
%     ds_inds = round(linspace(1, nPos, size(caEvents,1)));
%     trackingData_ds = cellfun(@(F) F(ds_inds,:), trackingData_cell, 'UniformOutput', 0);
    
    % method B. Interpolates positions based on interpolation vector with
    % same length of Ca2+ events. Allows smoothing of sudden position jumps.
    interpVec = linspace(1,nPos,size(cellTraces,2));    % Interpolation Vector
    trackingData_ds = cellfun(@(x) interp1(1:size(x,1), x, interpVec, 'linear'), trackingData_cell, 'UniformOutput', false);
    
    trackPoints = fieldnames(trackingData);
    trackingData = cell2struct(trackingData_ds, trackPoints, 1);
    posData = trackingData.scope;
    for i = 1:length(posData);
        behaviour_ts(i,1) = i/20; %Because tracking data is now sampled at 20Hz
    end
    posData = [behaviour_ts posData];
    
    earR = trackingData.earR;
    earL = trackingData.earL;
    earsTracking = [behaviour_ts earL earR];
    
    
    match_Im_to_Bh = [];
    
end


end



