%% Rescaling of tracking data.

% Function rescales tracking data according to arena dimensions and 
% saves data for each tracked point (body parts defined by user with DLC) 
% to the same folder as DLC output. This function also detects potential
% sessions with poor coverage. In cases of poor coverage the animals may not visit some edges of the
% arena. This will mess up the scaling of the data and is therefore
% better use the coordinates of the arena to define max and min coordinates. 

% USAGE
%   [trackingData, originalTrackingLimits]= upload.dlc_Import(path, param, trackingData_raw) 
%
%   param                             user parameters
%   trackingData_raw          DLC tracking data, cleaned.
%
% OUTPUT
%   trackingData                   structure with rescaled tracking data for each body part identified with DeepLabCut
%   originalTrackingLimits   limits of original tracking data.

% Written by Lucie Descamps and Miguel Carvalho 2020

function [trackingData, originalTrackingLimits] = rescalePath(path, param, trackingData_raw)

% get parameters used.
system = param.system;
refBdPt = param.refBdPt;
poorCovThr  = param.poorCovThr;
arenaDim = param.arenaDim; % arena side dimensions in cm.
mapLimits = param.mapLimits; % arena X and Y limits in cm.

%% detect potential sessions with poor coverage and define arena limits to be used in data rescaling.
originalTrackingLimits = [min(trackingData_raw.(refBdPt)(:, 1)) max(trackingData_raw.(refBdPt)(:, 1)) ...
    min(trackingData_raw.(refBdPt)(:, 2)) max(trackingData_raw.(refBdPt)(:, 2))];

if  originalTrackingLimits(2) - originalTrackingLimits(1) < poorCovThr || originalTrackingLimits(4) - originalTrackingLimits(3) <poorCovThr
    opts = struct('WindowStyle','modal', 'Interpreter','tex');
    warn = warndlg('Potential session with poor coverage!','Warning', opts);
    uiwait(warn);
    
    % upload video frame and plot path for visual inspection.
    [file, path] = uigetfile('.avi', 'Select behavior video file', path);
    selectedFile = fullfile(path,file);
    video = VideoReader(selectedFile);
    fig = figure;
    set(fig,'WindowStyle','normal')
    refFrame = read(video, 5);
    imshow(refFrame);
    set(gcf, 'MenuBar', 'none');
    set(gca, 'Visible', 'on');
    hold on
    plot(trackingData_raw.(refBdPt)(:, 1),trackingData_raw.(refBdPt)(:, 2),'color',[1 0 0]);
    
    % ask user which set of coordinates to use in data rescaling.
    answer = questdlg('What do you want to do?', 'Poor Coverage?', ...
        'Use Tracking Coordinates','Define Arena Coordinates','Define Arena Coordinates');
    
    % handle response
    switch answer
        case 'Use Tracking Coordinates'
            % % rescale tracking data based on tracked positions.
            xRescale = structfun(@(x) rescale(x(:,1), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
            yRescale =  structfun(@(x) rescale(x(:,2), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
            trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
            
        case 'Define Arena Coordinates'
            disp('Define arena limits');
            arenaShape = strsplit(param.arena);
            if ismember('cylinder', arenaShape)
                arenaROI = drawcircle;
                pos = customWait(arenaROI);
            else
                arenaROI = drawrectangle;
                pos = customWait(arenaROI);
            end
            xMin = min(arenaROI.Vertices(:,1));
            xMax = max(arenaROI.Vertices(:,1));
            yMin = min(arenaROI.Vertices(:,2));
            yMax = max(arenaROI.Vertices(:,2));
            originalTrackingLimits = [xMin xMax yMin yMax];
            close(fig);
            % rescale tracking data based on arena limits drawn by user.
            oldRange= xMax - xMin;
            newRange = arenaDim;
            xRescale = structfun(@(x) ((newRange * (x(:,1) - min(x(:, 1))))./ oldRange) + mapLimits(1), trackingData_raw, 'UniformOutput', false);
            yRescale = structfun(@(x) ((newRange * (x(:,2) - min(x(:, 2))))./ oldRange) + mapLimits(3), trackingData_raw, 'UniformOutput', false);
            trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
    end
else
    % rescale tracking data based on tracked positions.
    xRescale = structfun(@(x) rescale(x(:,1), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
    yRescale =  structfun(@(x) rescale(x(:,2), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
    trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
    
end

end

%% additional functions used to draw circles
function pos = customWait(temp)
% Listen for mouse clicks on the ROI
l = addlistener(temp,'ROIClicked',@clickCallback);
% Block program execution
uiwait;
% Remove listener
delete(l);
% Return the current position
pos = temp.Vertices;
end

function clickCallback(~,evt)
if strcmp(evt.SelectionType,'double')
    uiresume;
end
end







% %% Rescaling of tracking data.
% 
% % Function rescales tracking data according to arena dimensions and 
% % saves data for each tracked point (body parts defined by user with DLC) 
% % to the same folder as DLC output. This function also detects potential
% % sessions with poor coverage. In cases of poor coverage the animals may not visit some edges of the
% % arena. This will mess up the scaling of the data and is therefore
% % better use the coordinates of the arena to define max and min coordinates. 
% 
% % USAGE
% %   trackingData = upload.dlc_Import(trackingData_raw, arena, arenaDim, binWidth, mapLimits) 
% %
% %   trackingData_raw            DLC tracking data, cleaned.
% %   arena                               character array with information about the shape of the arena.
% %   arenaDim                        dimensions of the arena. Used to rescale tracked positions. 
% %   binWidth                         size of arena bins.
% %   mapLimits                       real arena limits in cm.
% %
% % OUTPUT
% %   trackingData                   structure with rescaled tracking data for each body part identified with DeepLabCut
% %   xMax
% %   xMin
% %   yMax
% %   yMin
% 
% % Written by Lucie Descamps and Miguel Carvalho 2020
% 
% function [trackingData, xMax, xMin, yMax, yMin] = rescalePath(path, trackingData_raw, arena, arenaDim, binWidth, mapLimits)
% 
% 
% %% detect potential sessions with poor coverage and define arena limits to be used in data rescaling. 
% trackingLimits = [min(trackingData_raw.Scope(:, 1)) max(trackingData_raw.Scope(:, 1)) ...
%     min(trackingData_raw.Scope(:, 2)) max(trackingData_raw.Scope(:, 2))];
% 
% if  trackingLimits(2) - trackingLimits(1) < 350 || trackingLimits(4) - trackingLimits(3) <350 
%     % the value of 350 is roughly the size of the 50cm circle in pixels, based on a sesssion. Only used as reference for potentially suspicious sessions.
%     opts = struct('WindowStyle','modal', 'Interpreter','tex');
%     warn = warndlg('Potential session with poor coverage!','Warning', opts);
%     uiwait(warn);
%     
%     % upload video frame and plot path for visual inspection.
%     [file, path] = uigetfile('.avi', 'Select behavior video file', path);
%     selectedFile = fullfile(path,file);
%     video = VideoReader(selectedFile);
%     fig = figure;
%     set(fig,'WindowStyle','normal')
%     refFrame = read(video, 5);
%     imshow(refFrame);
%     set(gcf, 'MenuBar', 'none');
%     set(gca, 'Visible', 'on');
%     hold on
%     plot(trackingData_raw.Scope(:, 1),trackingData_raw.Scope(:, 2),'color',[1 0 0]);
%     
%     % ask user which set of coordinates to use in data rescaling.
%     answer = questdlg('What do you want to do?', 'Poor Coverage?', ...
% 	'Use Tracking Coordinates','Define Arena Coordinates','Define Arena Coordinates');
% 
%     % handle response
%     switch answer
%         case 'Use Tracking Coordinates'
%             % % rescale tracking data based on tracked positions.
%             xRescale = structfun(@(x) rescale(x(:,1), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
%             yRescale =  structfun(@(x) rescale(x(:,2), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
%             trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
%             
%             xMax = max(trackingData_raw.Scope(:, 1));
%             xMin = min(trackingData_raw.Scope(:, 1));
%             
%             yMax = max(trackingData_raw.Scope(:, 2));
%             yMin = min(trackingData_raw.Scope(:, 2));
%             
%         case 'Define Arena Coordinates'
%              disp('Define arena limits');
%              arenaShape = strsplit(arena);
%              if ismember('cylinder', arenaShape)
%                  arenaROI = drawcircle;
%                  pos = customWait(arenaROI);
%              else
%                  arenaROI = drawrectangle;
%                  pos = customWait(arenaROI);
%              end
%              xMin = min(arenaROI.Vertices(:,1));
%              xMax = max(arenaROI.Vertices(:,1));
%              yMin = min(arenaROI.Vertices(:,2));
%              yMax = max(arenaROI.Vertices(:,2));
%              arenaLimits = [xMin xMax yMin yMax];
%              close(fig);
%              % rescale tracking data based on arena limits drawn by user.
%              oldRange= xMax - xMin;
%              newRange = arenaDim;
%              xRescale = structfun(@(x) ((newRange * (x(:,1) - min(x(:, 1))))./ oldRange) + mapLimits(1), trackingData_raw, 'UniformOutput', false);
%              yRescale = structfun(@(x) ((newRange * (x(:,2) - min(x(:, 2))))./ oldRange) + mapLimits(3), trackingData_raw, 'UniformOutput', false);
%              trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
%     end
% else
%     % rescale tracking data based on tracked positions.
%     xRescale = structfun(@(x) rescale(x(:,1), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
%     yRescale =  structfun(@(x) rescale(x(:,2), -arenaDim/2, arenaDim/2), trackingData_raw, 'UniformOutput', false);
%     trackingData = cell2struct(cellfun(@horzcat,struct2cell(xRescale),struct2cell(yRescale),'uni',0),fieldnames(xRescale),1);
%     
%     xMax = max(trackingData_raw.Scope(:, 1));
%     xMin = min(trackingData_raw.Scope(:, 1));
%     
%     yMax = max(trackingData_raw.Scope(:, 2));
%     yMin = min(trackingData_raw.Scope(:, 2));
% end
% 
% end
% 
% %% additional functions used to draw circles
% function pos = customWait(temp)
% % Listen for mouse clicks on the ROI
% l = addlistener(temp,'ROIClicked',@clickCallback);
% % Block program execution
% uiwait;
% % Remove listener
% delete(l);
% % Return the current position
% pos = temp.Vertices;
% end
% 
% function clickCallback(~,evt)
% if strcmp(evt.SelectionType,'double')
%     uiresume;
% end
% end