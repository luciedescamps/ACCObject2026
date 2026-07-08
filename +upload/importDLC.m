%% Upload tracking data from DeepLabCut and perform some cleaning

% User identifies the file with DLC tracking data. The files comes in the form of a .csv file
% with timestamps and coordinates for several tracked points. Spurious
% detection errors in DLC are cleaned
%
% USAGE
%   [trackingData_raw, varargout] = upload.importDLC(path, param)
%   path                                 path to subject/session folder containing RAW tracking data extracted from DeepLabCut.
%   param                             user parameters
%
% OUTPUT
%   trackingData_raw            structure with rescaled tracking data for each body part identified with DeepLabCut
%   varargout                         optional output of behavior timestamps (only for ucla scopes)

% Written by Lucie Descamps and Miguel Carvalho 2020

%%
function [trackingData_raw, varargout] = importDLC(path, param)

% get parameters used.
system = param.system;

% upload behavior timestamp file
if ~strcmp(param.system, 'nvista')
    posT = fullfile(path,'\behaviour\timeStamps.csv');
    posT = readcell(posT);
    posT = cell2mat(posT(2:end,:));
    
    %Trim some portion of the start of the recording (ex if mouse is not in the
    %arena at the start)
    if param.trim_tracking == 1
        tsDiff = round(mean(diff(posT(:, 2))), -1);
        dlc_frameRate = (1/tsDiff)*1000;
        
        frame_idx_trim = dlc_frameRate*param.TrimSec;
        
        posT(1:frame_idx_trim+1,:) = [];
    end

varargout{1}= posT(:, 1:2);
end
% upload tracking .csv file and import data
[file, path] = uigetfile('.csv', 'Select DeepLabCut tracking file', path);
selectedFile = fullfile(path,file);
opts = detectImportOptions(selectedFile);
opts.VariableNamesLine = 2;
dlcDataRaw = readtable(selectedFile, opts);

trackingInfo = dlcDataRaw.Properties.VariableNames(2:end);
bodyparts = trackingInfo(1:3:end);

dlcData = table2array(dlcDataRaw);

if param.trim_tracking == 1
    dlcData =  dlcData(frame_idx_trim+2:end, 2:end);
else
    dlcData =  dlcData(:, 2:end);
end

%% get tracking data for different bodyparts and correct tracking jumps.


% detect tracking jumps
jumpThreshold = 19; % in pixels
for bp = 1: length(bodyparts)
    disp(sprintf('Processing tracking data for %s', bodyparts{bp}));
    indx = double(contains(trackingInfo(1 , :), bodyparts(bp)));
    bodypartData = dlcData(:, indx==1);
    jumpX = find(abs(diff(bodypartData(:, 1)))> jumpThreshold)+1;
    jumpY = find(abs(diff(bodypartData(:, 2)))> jumpThreshold)+1;
    negTracking = find(bodypartData(:, 1:2) < 0);
    
    if exist ('negTracking') &&  ~isempty(negTracking)
        for i = size(negTracking,1):-1:1
            if negTracking(i) > size(bodypartData, 1)
                negTracking(i) = [];
            end
        end    
    end
    
    jumpPos = unique([jumpX; jumpY; negTracking]);
    
    if ~isempty (jumpPos)
        % clean X positions
        for i = 1: length(jumpPos)
            if isnan(bodypartData(jumpPos(i), 1:2))
                % this accounts for previously cleaned positions as part of back and forth jumps.
                continue
                bodypartData(jumpPos(i), 1:2) = NaN;
            elseif i == length(jumpPos) || jumpPos(i+1) -  jumpPos(i) > 20
                % the rest of tracking errors end up being sudden jumps
                % that don't return to original position. We clean only the first frame.
                bodypartData(jumpPos(i), 1:2) = NaN;
            elseif jumpPos(i+1) -  jumpPos(i) <= 20
                % look for short jumps and return to previous position. Clean positions within that interval.
                bodypartData(jumpPos(i): jumpPos(i+1), 1:2) = NaN;
            end
        end
        
        % interpolate position data to remove NaNs.
        bodypartData = interp1(1: size(bodypartData, 1), bodypartData, 1: size(bodypartData, 1), 'pchip');
    end
    writematrix(bodypartData, fullfile(path, sprintf('Tracking data - %s.xlsx', bodyparts{bp})));
    trackingData.(bodyparts{bp}) = bodypartData;
end

%% Save path plots and position data

% save a sample figure with raw tracking data for each body part.
fig = figure('Position', [10, 50, 1400, 750]);
for bp = 1 : length(bodyparts)
    indx = double(contains(trackingInfo(1 , :), bodyparts(bp)));
    rawPos = dlcData(:, indx==1);
    subplot(ceil(length(bodyparts)/2), 4, 2*bp-1);
    plot(rawPos(:,1), rawPos(:,2),'color', [0.7 0.4 0.4], 'LineWidth', 0.1);
    title(sprintf('Raw Tracking for %s', bodyparts{bp}));
    
    cleanPos = trackingData.(bodyparts{bp});
    subplot(ceil(length(bodyparts)/2), 4, 2*bp);
    plot(cleanPos(:,1), cleanPos(:,2),'color', [0.4 0.7 0.4],'LineWidth', 0.1);
    title(sprintf('Clean Tracking for %s', bodyparts{bp}));
end
w = waitforbuttonpress;
%saveas(gcf, fullfile(path, sprintf('Sample Tracking')), 'png');
close(fig);

trackingData_raw = structfun(@(x) x(:,1:2), trackingData,  'UniformOutput', false);
end