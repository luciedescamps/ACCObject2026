%Script to run CLEAN a bit more easiliy! 

%% Input data to put in the CLEAN obj 

for i = 1:11
    display(i);
    [file,path] = uigetfile('.mat', 'Select cell extraction and manual sorting outputs', 'Multiselect', 'on');
    for ffile = 1:2
        load(char(fullfile(path, file(ffile))));
    end
    
    inputImages{i} = extractAnalysisOutput.filters;
    inputSignals{i} = extractAnalysisOutput.traces;
    inputTargets{i} = validEXTRACT;
    
    [filename pathname] = uigetfile('*.h5', 'Select the dfof video');
    moviename = fullfile(pathname,filename);
    inputMovieList{i} = moviename;
    
    clear validEXTRACT extractAnalysisOutput moviename
end

%% Initialize clean object and train the model
cleanObj = ciapkg.clean();

% Input training sessions and [optional] classification sessions.
cleanObj.inputImages = inputImages; % cell array {1 nSessions} of [x y nSignals] matrices
cleanObj.inputSignals = inputSignals; % cell array {1 nSessions} of [nSignals frames] matrices
cleanObj.inputTargets = inputTargets; % cell array {1 nSessions} of [1 nSignals] vectors
cleanObj.inputMovieList = inputMovieList; % cell array {1 nSessions} of char strings
cleanObj.readMovieChunks = 1; % Read movies from disk.
cleanObj.inputDatasetName = '/1'; % HDF5 dataset name, change as needed.

% Index of sessions to use for training (assuming we input 11 sessions)
cleanObj.trainIdx = [1 2 3 4 5 6 7 8 9 10 11];

% Compute features for all training and test data. Allows faster training and testing using different session combinations.
cleanObj.computeFeatures;

% Create CLEAN model based on training set
cleanObj = train(cleanObj,cleanObj.inputImages,cleanObj.inputSignals,cleanObj.inputTargets);

% Export the CLEAN object and trained model
cleanObj.export;

%% Import data you want to classify


for i = 1:15
    display(i);
    [file,path] = uigetfile('.mat', 'Select cell extraction output');
    load(fullfile(path,file));
    
    testImgs{i} = extractAnalysisOutput.filters;
    testSignals{i} = extractAnalysisOutput.traces;
    testMovies{i} = extractAnalysisOutput.file;
    
    clear extractAnalysisOutput 
end

%% Classify new data

for i = 1:15
    %cleanObj.classify(testImgs,testSignals,inputMovieList,testMovies);
    cleanObj = classify(cleanObj,testImgs{i},testSignals{i});
    % Get majority vote and SVM consensus classifications for that folder in.
    majorityVoteLabels{i} = cleanObj.classifyStruct.classifications;
    svmLabels{i} = cleanObj.classifyStruct.classificationsModelDecisions;
end

%% Visualize the sorted cells - majorityVoteLabels

for i = 1:15
    cell_idx = find(majorityVoteLabels{i} == 1);
    goodfilters =testImgs{i}(:, :, cell_idx);
    
    ind = find(goodfilters);
    [i1, i2, i3] = ind2sub(size(goodfilters), ind);
    figure();
    scatter(i1, i2, 5, i3);
    
    clear cell_idx goodfilters
end

%% Visualize the sorted cells - svm labels

for i = 1:15
    cell_idx = find(svmLabels{i} == 1);
    goodfilters =testImgs{i}(:, :, cell_idx);
    
    ind = find(goodfilters);
    [i1, i2, i3] = ind2sub(size(goodfilters), ind);
    figure();
    scatter(i1, i2, 5, i3);
    
    clear cell_idx goodfilters
end

%% Save the sorting results 

k = 1;
for i = 1:15
    
    extractAnalysisSorted_Clean = majorityVoteLabels{i};
    extractAnalysisSorted_Clean = extractAnalysisSorted_Clean';
    extractAnalysisSorted_Clean_svm = svmLabels{i};
    extractAnalysisSorted_Clean_svm = extractAnalysisSorted_Clean_svm';
    fname = sprintf('extractAnalysisSorted_Clean_C%d.mat', k);
    filename = sprintf('extractAnalysisSorted_Clean_svm_C%d.mat', k);
    save(fname, 'extractAnalysisSorted_Clean');
    save(filename, 'extractAnalysisSorted_Clean_svm');
    clear extracttAnalysisSorted_Clean fname extractAnalysisSorted_Clean_svm filename;
    k = k +2;
end

