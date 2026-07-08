%%First need to add timestamps to trackingData (the 108k frames one!!! Not downsampled one! )

for i = 2:session_duration_tracking
    newvector(i,1) = i/60; %Because tracking data is 60Hz
end
trackingData = structfun(@(x) horzcat(newvector(:, 1), x), trackingData, 'UniformOutput', false);


%Then filter the calcium traces to keep only the significant transiants
%(see amplitude_draft)
%Then we need to add timestamps to these guys too: 

session_duration_ca = length(cellTraces_filtered_3std_zscored); %How long is the imaging data
newvector = zeros(session_duration_ca,1); %Pre assignment

for i = 2:session_duration_ca
    newvector(i,1) = i/20;%20 because it was at 20Hz 
end
newvector = newvector'; %This is a vector of seconds corresponding to frames indices 

cellTraces_filtered_3std_zscored_ts = [newvector;cellTraces_filtered_3std_zscored]; %Concatenating the timestamps to the filtered traces
cellTraces_filtered_3std_zscored_ts = cellTraces_filtered_3std_zscored_ts'; %Flipping the axis so it's frames x cells (nb: the first "cell" is actually the timestamps, so there are cells-1 cells in the recording)
cellTraces_filtered_3std_zscored_ts_aligned(:,1)= aligner.aToB(cellTraces_filtered_3std_zscored_ts(:,1), 'piecewise') %We are using the aligner function to convert the second indices in imaging referential into behaviour referential , Need to double check this is the correct way around


%We can then use a knnsearch to find the corresponding timestamps between
%behaviour and imaginng
match_Im_to_Bh = knnsearch(trackingData.Scope(:, 1), cellTraces_filtered_3std_zscored_ts_aligned(:,1));
match_Bh_to_Im = knnsearch(cellTraces_filtered_3std_zscored_ts_aligned(:,1), trackingData.Scope(:, 1));



%Or we can ds the behaviour from 60 to 20 hz, then match knn 

 test_ds_behaviour = downsample(trackingData,Scope, 3); %3 is the downsampling factor because 60/20 = 3
 

session_duration_tracking = length(test_ds_behaviour );
newvector = zeros(session_duration_tracking,1);

for i = 2:session_duration_tracking
    newvector(i,1) = i/20;
end

test_ds_behaviour = [newvector test_ds_behaviour];
match_Im_to_Bh2 = knnsearch(test_ds_behaviour(:, 1), cellTraces_filtered_3std_zscored_ts_aligned(:,1));
match_Bh_to_Im2 = knnsearch(cellTraces_filtered_3std_zscored_ts_aligned(:,1), test_ds_behaviour(:, 1));

%Then once you have the knn I would use all these indexes to create
%a new behav vector that only contains the frames index of match_im_to_bh
%Might be a bit more precise to not ds the behaviour vector 

new_track_data =trackingData.Scope(match_Im_to_Bh,:); %And voila!!! :D 