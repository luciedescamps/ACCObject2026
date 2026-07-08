%% Script to find when is the Miniscope LED first on.

%When doing calcium imaging experiments with the miniscope, we first start
%the behavioural cameras, then the imaging scope

%A small red LED is flashing everytime the scope is taking a frame. We need
%to identify when is the first LED flash happening to know when the scope
%was first active, in order to sync both imaging and behavioural data. 

% Written by Lucie Descamps 2020



function [first_frame] = identifyLED(path)
%% Detect threshold 
%Do this every time you are changing the lighting conditions in the rig. If
%the lighting conditions are stable you can just do it once.

    [filename pathname] = uigetfile('*.avi', 'Select the behavioural video', path);
    moviename = fullfile(pathname,filename);
    vidObj = VideoReader(moviename);
    numFrames = ceil(vidObj.FrameRate*vidObj.Duration);
    
    
    refFrame_LED_OFF = read(vidObj,5);
    refFrame_LED_ON = read(vidObj, 700);
    
    f1 = figure;
    imagesc(refFrame_LED_OFF);
    disp('Draw an outline around the LED - LED should be OFF');
    truesize;
    set(gcf, 'MenuBar', 'none')  % gcf = get current frame
    set(gca, 'DataAspectRatioMode', 'auto') %gca = get current axes
    set(gca, 'Position', [0 0 1 1 ])
    drawnow
    led_off = roipoly(refFrame_LED_OFF);
    close(f1);

    f1 = figure;
    imagesc(refFrame_LED_ON);
    disp('Draw an outline around the LED - LED should be ON');
    truesize;
    set(gcf, 'MenuBar', 'none')  % gcf = get current frame
    set(gca, 'DataAspectRatioMode', 'auto') %gca = get current axes
    set(gca, 'Position', [0 0 1 1 ])
    drawnow
    led_on = roipoly(refFrame_LED_ON);
    close(f1);
    
    ave_pixel_value_on = mean2(refFrame_LED_ON(logical(led_on)))
    ave_pixel_value_off = mean2(refFrame_LED_OFF(logical(led_off)))
    
    threshold = ave_pixel_value_on-5;



%% Find first frame with LED ON

refFrame= read(vidObj,500);
imagesc(refFrame);
disp('Draw an outline around the LED');
truesize;
set(gcf, 'MenuBar', 'none')  % gcf = get current frame 
set(gca, 'DataAspectRatioMode', 'auto') %gca = get current axes
set(gca, 'Position', [0 0 1 1 ])
drawnow
bw_led = roipoly(refFrame);

for frame_idx = 1:numFrames
    frame = read(vidObj,frame_idx);
    ave_pixel_value = mean2(frame(logical(bw_led)));
    
    if ave_pixel_value > threshold
        sprintf('LED is on for the first time at frame number %d\n', frame_idx')
        first_frame = frame_idx;
        break 
    end
    frame_idx = frame_idx+1;
end


