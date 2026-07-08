%%This script is to create the input vectors needed to align the imaging
%%and behaviour data.
%Vecbrilho section heavily relies on code from Rafael Pedrosa,
%Battaglia lab.
%Alignment section heavily relies on code from Rich Gardner, Moser group.

%% Create pulse vector from the behavioural video
function  [aligner, pulseTrains, remove] = createVectors
[filename pathname] = uigetfile('*.avi', 'Select the video');
moviename = fullfile(pathname,filename);
vidObj = VideoReader(moviename);

refFrame= read(vidObj,800); %Change this number to a frame number where you can clearly identify the LED position
imagesc(refFrame);
truesize;
set(gcf)  % gcf = get current frame
set(gca, 'DataAspectRatioMode', 'auto') %gca = get current axes
set(gca, 'Position', [0 0 1 1 ])

clear vetorbrilho

NumberOfFrames = ceil(vidObj.FrameRate*vidObj.Duration);

answer = questdlg('Which environment are you using?', ...
	'Environment', ...
	'Circular (Lucie)','Rectangular (Elise)', 'Maybe');

% Handle response
switch answer
    case 'Circular (Lucie)'
        envtype = 1;
    case 'Rectangular (Elise)'
        envtype = 2;
end

for i = 1:NumberOfFrames
    i
    video = read(vidObj,i);
    
    %X and Y values of the LED in the reference frame. Change if needed (depend on the video)
    if envtype ==1
       %xled = 649:653;
       xled= 660:665; %both should be 5x5 matrices
       %xled = 660:665;
       %yled = 882:886;
      %yled = 904:909;
       yled=920:925;
       %yled = 910:915;
    end
    
    if envtype ==2
        xled=237:239;
        yled=465:469;
    end
    
    %Your blinking pattern
    vetorbrilho(i) = squeeze(mean(mean(mean(video(yled,xled,:)))));
end

%Indexing LED blinking
vecbrilho = vetorbrilho;
vecbrilho = vecbrilho-min(vecbrilho);
blinking_vector=(vecbrilho > 220) ; %Usually 50! had to change for bubbles H1 
plot(blinking_vector,'o') %Visually check if the pattern looks correct

FrameRate = 60; %Input here the frame rate of your behavioural camera
rise_index = strfind(blinking_vector,[0 1]); %finds when it's switching from 0 to 1
fall_index = strfind(blinking_vector,[1 0]); %finds when it's switching from 1 to 0

pulseTracking.rise = rise_index';
pulseTracking.fall = fall_index';

pulseTimes.tracking.rise = pulseTracking.rise / FrameRate; %converts frame index to second index
pulseTimes.tracking.fall = pulseTracking.fall / FrameRate; %converts frame index to second index

%% Create pulse vector from the Saleae file

%To do before running this section: you need to open the logic file from the session you want to align, then export it as a csv for only the LEDs
%channel. The input will be a csv of 2 columns, the first column is the
%timestamps in seconds and the second column is the status of the LED
%channel. There is a new row every time the status is changing (ie. either
%rising or falling).

[filename pathname] = uigetfile('*.csv', 'Select Saleae sync file');
syncchannel = fullfile(pathname,filename);
syncchannel = csvread(syncchannel);

%Find the timestamp of the nVista frame.
start_nvista = find(syncchannel(1:end,2),1,'first');
timestampFirstFrame = syncchannel(start_nvista,1);

%Load Saleae LED file
[filename pathname] = uigetfile('*.csv', 'Select Saleae LED file');
saleae = fullfile(pathname,filename);
saleae = csvread(saleae);

%Find first saleae pulse after timestamp first frame!
first_saleae_pulse_after_nvista = find( saleae > timestampFirstFrame, 1 );

%So how many nvista frame will we need to ditch?
remove = saleae(first_saleae_pulse_after_nvista,1) - timestampFirstFrame;


%Change t=0 on the saleae file so it starts when nVista starts
for i = 1:length(saleae)
    saleae(i,1) = (saleae(i,1) - timestampFirstFrame);
end

indices = find(saleae(:,1) < 0);
saleae(indices,:) = [];

%create pulse vectors for saleae
saleae_rise = saleae(find(saleae(:,2) == 1));
saleae_fall = saleae(find(saleae(:,2) == 0));

pulseTimes.saleae.rise = saleae_rise;
pulseTimes.saleae.fall = saleae_fall;

%The 2 rise and fall vectors need to be the same length, and the fall value
%always needs to be greater than the rise value at the corresponding index.
%However one of the 2 vectors will always be bigger than the other since it
%will contain the status of the channel at timestamp = 0s. So we are just
%going to chop chop the one that contains this value :-)
%But that also means that the t0 for the nvista time stamps is going to be
%modified, so we need to calculate which frame it is at the first saleae
%event; then convert this frame index in s, substract it from caEvents

if size(pulseTimes.saleae.fall, 1) > size(pulseTimes.saleae.rise, 1)
    pulseTimes.saleae.fall(1,:) = [];
end

if size(pulseTimes.tracking.fall, 1) > size(pulseTimes.tracking.rise, 1)
    pulseTimes.tracking.fall(1,:) = [];
end

if size(pulseTimes.tracking.rise, 1) > size(pulseTimes.tracking.fall, 1)
    pulseTimes.tracking.rise(end,:) = [];
end

if size(pulseTimes.saleae.rise, 1) > size(pulseTimes.saleae.fall, 1)
    pulseTimes.saleae.rise(end,:) = [];
end


%% Visualize the vectors

figure();

signalNames = ["Saleae", "tracking LED pulses"];
pulseData{1} = pulseTimes.saleae;
pulseData{2} = pulseTimes.tracking;


for n = 1:2
    
    name = signalNames(n);
    t = pulseData{n};
    
    % Create the PulseTrain object.
    %
    % We create a PulseTrain by calling the class constructor method
    % ("sync.PulseTrain"), with the pulse rise and fall times as the two
    % arguments. This creates a PulseTrain, which we will call "pt".
    pt = sync.PulseTrain(t.rise, t.fall);
    pt.name = name;
    pulseTrains(n) = pt;
    
    % Plot the pulse times.
    %
    % It's possible to plot a PulseTrain object with the simple command
    % "plot(pt)".
    subplot(4, 1, n);
    plot(pt);
    xlim(pt.tRise(1) + [-1 30]);
    xlabel("Time / s");
    yticks([0, 1]);
    yticklabels(["OFF", "ON"]);
    title(name);
end

%%Align the vectors

% ALIGN THE SYNC PULSES
%
% Since the sequence of pulses has random intervals, the temporal pattern
% of intervals in any given chunk of the recording should be unique. This
% is the principle which we use to align the two recordings. We divide the
% first recording into chunks of 10 pulses, and for each chunk we search
% for the most similar chunk in the second recording.

% Create a PatternAligner object
%
% The PatternAligner class is used to handle the process of aligning two
% PulseTrains to each other by matching the unique temporal sequences of
% intervals between the pulses. We create a PatternAligner by calling its
% class constructor with the two PulseTrain objects as the arguments.
aligner = sync.PatternAligner(pulseTrains(1), pulseTrains(2));

% Plot the aligned pulses.
%
% We can create a convenient plot of the aligned pulses directly from a
% PatternAligner object, using the command "plot(aligner)". Here we will
% see the two pulse trains plotted in parallel, with the second pulse train
% (tracking) shifted in time such that it is aligned with the first pulse
% train (spikes).
subplot(2, 1, 2);
plot(aligner);
xlim([0, 30]);

% We can check that the matching has worked correctly by making a scatter
% plot of the matched pulse times. They should form a near-perfect straight
% line (we expect a small amount of error due to quantization of the sample
% timestamps).
%
% The field 'tAMatched' contains a vector of timestamps of matched pulses,
% in time frame A (which is 'ephys' in this tutorial). tBMatched contains
% the corresponding pulse times in time frame B ('tracking' in this
% tutorial).
figure();
x = aligner.tAMatched;
y = aligner.tBMatched;
plot(x, y, 'k.');
xlabel('t_{saleae} / s', 'interpreter', 'tex');
ylabel('t_{tracking} / s', 'interpreter', 'tex');
title('Matched pulse times');

%Save the aligner and pulseTrains
% save(fullfile(path, 'aligner'), 'aligner');
% save(fullfile(path, 'pulseTrains'), 'pulseTrains');

end
%% Align cell events data to tracking data  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAP TIMES BETWEEN THE TWO TIME FRAMES
%
% The purpose of aligning the sync pulses is so that we can relate data 
% recorded in one time frame to data recorded in the other time frame.
%
% For example, we have the spike times of a neuron recorded in the "ephys"
% (saleae)
% time frame, and we want to relate that to the animal's position in the 
% "tracking" time frame. We can do this by interpolating any arbitrary 
% time point between time frames, using the matched sync pulse times.

% First let's simulate some 2D tracking data
% nTracking = 100*25;
% trackingTimeLims = [1770, 1850];
% tracking.t = linspace(trackingTimeLims(1), trackingTimeLims(2), nTracking);
% tracking.x = smooth(cumsum(randn(nTracking, 1)), 50);
% tracking.y = smooth(cumsum(randn(nTracking, 1)), 50);


%Import tracking data from imaGUI


% Now we'll generate some random spike times that we will match up with the
% tracking data
% nSpikes = 30;
% spikeTimes = sort(rand(nSpikes, 1)*50);

% We can map the spike times from their original time frame to the tracking
% time frame using the PatternAligner "aToB" method. This method takes a 
% vector of times specified in time frame A (events/saleae), and interpolates them 
% across to time frame B (tracking).
%
% (If you're interested, this interpolation is calculated via a linear 
% regression fitted to the two series of sync pulses.)
%spikeTimesTracking = aligner.aToB(spikeTimes);
% spikeTimesTracking = aligner.aToB(spikeTimes, 'piecewise') %% Use this one to be more precise!!

%save it! 

% For both timeseries (spikes and tracking), plot the sync pulses and the
% simulated data
% figure();
% subplot(4, 1, 1);
% plot(pulseTrains(1));
% 
% x = spikeTimes;
% x = [x, x, nan(nSpikes, 1)]';
% x = x(:);
% y = nan(3, nSpikes);
% y(1, :) = 0;
% y(2, :) = 1;
% y(3, :) = nan;
% hold on
% plot(x(:), y(:)+2, 'r', "lineWidth", 1);
% xlim([0, 50]);
% ylim([-0.5, 3.5]);
% yticks([0.5, 2.5]);
% yticklabels(["TTL pulses", "Spike times"]);
% xlabel("t_{ephys} / s");
% title("Ephys data")
% 
% subplot(4, 1, [3, 4]);
% plot(pulseTrains(2));
% x = tracking.t;
% y = zscore([tracking.x, tracking.y]);
% hold on
% h = plot(x, y + 3);
% legend(h, ["position X", "position Y"]);
% ylim([-0.5 6]);
% xlim(trackingTimeLims);
% yticks([0.5, 3]);
% yticklabels(["LED pulses", "Position"]);
% xlabel("t_{tracking} / s");
% title("Tracking data")
% 
% x = tracking.x;
% y = tracking.y;
% t = tracking.t;
% spikePos = interp1(t, [x, y], spikeTimesTracking);
% 
% clear h
% figure
% h(1) = plot(x, y, 'k');
% hold on
% h(2) = plot(spikePos(:, 1), spikePos(:, 2), 'r.', 'markerSize', 10);
% title("Path and spike positions");
% axis equal tight
% title("Tracking data");
% xlabel("x");
% ylabel("y");
% legend(h, ["Animal path", "Spike position"]);