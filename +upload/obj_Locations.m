%% Upload tracking video, identify area of the recording arena and object locations

% User identifies the tracking video file. The function loads a video frame from input behaviour video
% where user can identify the centre of objects or other features of the
% environment.

% USAGE
%   [objCoordinates, video, ROI] = upload.obj_Locations(path)
%   path                                 path to subject/session folder containing RAW tracking video
%
% OUTPUT
%   objCoordinates               matrix with X and Y coordinates for each object marked
%   video                               tracking video
%   ROI                                 cell array with description of objects

% Written by Lucie Descamps and Miguel Carvalho 2020

%%

function [objCoordinates, video, ROI] = obj_Locations(path)

% upload video file and select a video frame
[file, path] = uigetfile('.avi', 'Select behavior video file', path);
selectedFile = fullfile(path,file);
video = VideoReader(selectedFile);


% load video frame.
fig = figure;
set(fig,'WindowStyle','normal')
refFrame = read(video, 5);
imshow(refFrame)
set(gcf)

% identify regions of interest such as objects.
list = {'Object 1','Object 2', 'Object 3', 'Object 4'};
[roi,tf] = listdlg('ListString', list, 'PromptString','Select ROIs:', 'ListSize', [250,300]);
ROI = cell(1, length(roi));

% identify the centre of object locations for each object identified.
for r = 1 : length(roi)
    ROI{r} = list{r};
    disp(sprintf('Mark the centre of %s ', list{r}))
    [xObj(r), yObj(r)] = ginput(1); hold on
    plot(xObj(r), yObj(r),'o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', [0.9882 0.4941 0.8588]);
end

close(fig);
objCoordinates = horzcat(xObj', yObj');

end