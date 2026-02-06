% This code should do the tree classifier using position labels and any
% coarse graining or regular data. Needs position path, rois of interest,
% coarse-grained data and single-cell data.
%Written by Wesley Clawson with modifications from Lucie Descamps


% This first section loads in data
clear
addpath(genpath('Code'));
disp('Loading Coarse Data...');

% Which mouse file
mouse = 'Bubbles';

% Where the coarse grained data is
coarsepath = ['W:\luciede\NTNU\Miniscope-rig\2022\Tracker results ACC mice\' mouse '\'];

% Path to the generated ROIs
roipath = ['W:\luciede\NTNU\Codes\Wes\' mouse '\roi2\'];

% The  big post imaGUI file
postgui = ['W:\luciede\NTNU\Miniscope-rig\2022\Tracker results ACC mice\' mouse '\post_ImaGui_95_HalfSessions.mat'];

% Load in needed data
 load([coarsepath 'Coarse_Session.mat']);  %Bubbles
 %load([coarsepath 'Coarse_Session_S1S10.mat']); %other mice

% What sessions
sessions = [1:3 5:10];

% Load in position data
positionpath =  ['W:\luciede\NTNU\Codes\Wes\' mouse '\'];
load([positionpath 'PosData' mouse '.mat']);

%% - Create label curves
%   You don't need to run this if you've run it before, and saved it. Just
%   proceed to the next section, where you load them in. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI CURVES %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numsesh = length(sessions);
label_curves = cell(1, numsesh);

% Build ROI curves for every file in each animal ID folder
for dt = 1:numsesh
    % Select the dataset (absent, baseline, etc)
    sesh_id = sessions(dt);

    % Load the reclen - /full recording
    if dt == 1 | dt ==  2 | dt == 3 
    rec_len = size(Coarse_Session.(['S' num2str(sesh_id)]).K_1.Data, 1);
    end


    %Load the reclen - /half sess
    if   dt == 4  ||  dt ==  5 || dt == 6 || dt == 7 || dt == 8 || dt ==9 
    rec_len = floor(size(Coarse_Session.(['S' num2str(sesh_id)]).K_1.Data, 1)/2);
    end 


    % Load the position data in
    posdata = AllPosData{1, sesh_id};
    if dt == 1 | dt ==  2 | dt == 3 
        pos_len = size(posdata, 1);
    end

    if dt == 4  | dt ==  5 || dt == 6 || dt == 7 || dt == 8 || dt ==9 
        pos_len = floor(size(posdata, 1)/2);
    end

    xpos = posdata(1:pos_len,1);
    ypos = posdata(1:pos_len,2);

    % Create labels of ROIs
    label_curve = cell(size(xpos, 1), 1);       % init
    label_curve(:, 1) = {'none'};               % fill with 'none' for now

    % We need to build a ROI curve for each case
    for zone = [1, 2, 3, 4, 6] %1:6 %Removed 5 so it doesnt merge left and right

        if dt == 1 || 2
            sesh_id = sessions(3);
        end

        switch zone
            case 1 % LEFT
                tempvar = 'blurmaskL';
                blurmask = load([roipath '\S' num2str(sesh_id) '_left.mat'], tempvar).(tempvar);
            case 2 % RIGHT
                tempvar = 'blurmaskR';
                blurmask = load([roipath '\S' num2str(sesh_id) '_right.mat'], tempvar).(tempvar);
            case 3 % LEFT EDGE
                tempvar = 'blurmaskleftedge';
                blurmask = load([roipath '\S' num2str(sesh_id) '_leftedge.mat'], tempvar).(tempvar);
            case 4 % RIGHT EDGE
                tempvar = 'blurmaskrightedge';
                blurmask = load([roipath '\S' num2str(sesh_id) '_rightedge.mat'], tempvar).(tempvar);
            case 5 % BOTH
                tempvar = 'blurmaskboth';
                blurmask = load([roipath '\S' num2str(sesh_id) '_both.mat'], tempvar).(tempvar);
            case 6 % MAZE EDGE
                tempvar = 'blurmaskmazeedge';
                %blurmask = load([roipath '\S' num2str(sesh_id) '_mazeedge.mat'], tempvar).(tempvar);
                tempx = round(xpos);
                tempy = round(ypos);
                tblur = zeros([max(tempy), max(tempx)]);
                mask1 = false(size(tblur));
                mask2 = false(size(tblur));
                mask3 = false(size(tblur));
                for i = 1:length(tempx)
                    tblur(tempy(i), tempx(i)) = 1;
                end

                figure('Name', 'MazeE');
                imagesc(tblur);
                title('Draw ellipse that touches the edges');
                h = drawcircle;
                pause;
                % mask1 = mask1 | createMask(h);
                % h2 = h;
                % h2.SemiAxes = h2.SemiAxes+(h2.SemiAxes*0.1);
                % mask2 = mask2 | createMask(h2);
                % h3 = h;
                % h3.SemiAxes = h3.SemiAxes-(h3.SemiAxes*0.1);
                % mask3 = mask3 | createMask(h3);

                mask1 = mask1 | createMask(h);
                h2 = h;
                h2.Radius = h2.Radius+(h2.Radius*0.1);
                mask2 = mask2 | createMask(h2);
                h3 = h;
                h3.Radius = h3.Radius-(h3.Radius*0.1);
                mask3 = mask3 | createMask(h3);

                edgemask = xor(mask2, mask3);
                blurmask = imgaussfilt(double(edgemask), 10);
                blurmask = flipud(blurmask);
                close('MazeE');
        end

        % Now, make the ROI curve %%%%%%%%%%%%%%%%%%%%
        % use 'griddata.m' to make a grid using the the same space as
        % blurmask, with the same values. Then, interpolate what x and y
        % should be as if blurmask was a 2D surface.
        x = 1:size(blurmask, 1);
        y = 1:size(blurmask, 2);
        v = flipud(blurmask);
        xq = griddata(y, x, v, xpos, ypos);
        roicurve = xq;

        label_curve(roicurve >= 0.75) = {tempvar};
    end

    label_curves(1, dt) = {label_curve};
end

disp("Don't forget to save");
save([positionpath 'Label_curves_separated_S1S10_1515.mat'], 'label_curves')

%% - Load in the label curves if you have them, if not, run the section above.
% which are time long and at each time point, a label of where the mouse
% is. 

% Load in needed data, should be 1 x sessions (10), long. 
load([positionpath 'Label_curves_separated_S1S10_1515.mat']);
%% - This will do decoding within an individual session.
tic

maxK =11;  

ObjSeshIDs = [1 2 3 5 6 7 8 9 10]; % <--- this is the important dude for decode
sessions = ObjSeshIDs;
num_sessions = length(sessions);
num_iter = 1;  % < ----------- this is how many times you'll run the trees
TreeResults_CV = [];
progressbar('Sessions', 'K');
for SI = 1:num_sessions
    % Grab the session ID
    sesh_id = sessions(SI);
    % Get the label curves
    lb1 = label_curves{1, SI};
    % Loop over K
    for K = 1:maxK
        if K > length(fieldnames((Coarse_Session.(['S' num2str(sesh_id)])))) -1 %If this K iteration  doesn't exist, continue, then eventually loops through SI
            continue
        end

        if sesh_id== 1 || sesh_id== 2 || sesh_id== 3 %whole sess because it's 15min anyway
            data1 = Coarse_Session.(['S' num2str(sesh_id)]).(['K_' num2str(K)]).Data;
        end
        if sesh_id ==  5 || sesh_id== 6 || sesh_id== 7 || sesh_id== 8 || sesh_id== 9 || sesh_id== 10
            %half sess so it's 15min too
            rec_len = floor(size(Coarse_Session.(['S' num2str(sesh_id)]).K_1.Data, 1)/2);
            data1 = Coarse_Session.(['S' num2str(sesh_id)]).(['K_' num2str(K)]).Data(1:rec_len,:);
        end

        classnames = {'blurmaskL', 'blurmaskR', 'blurmaskleftedge', 'blurmaskrightedge', ...
            'blurmaskboth', 'blurmaskmazeedge', 'none'};
        
            lbr = lb1(randperm(size(lb1, 1)), 1);
            cv1 = cvpartition(lb1, 'HoldOut', 0.2, 'Stratify', true);
            trainInds = training(cv1);
            testInds = test(cv1);
            cvr = cvpartition(lb1, 'HoldOut', 0.2, 'Stratify', true);
            trainIndsR = training(cvr);
            testIndsR = test(cvr);
            XTrain = data1(trainInds, :);
            YTrain = lb1(trainInds);
            XTest = data1(testInds, :);
            YTest = lb1(testInds);
            XTrainR = data1(randperm(size(data1, 1)),  :);
            XTrainR = XTrainR(trainInds, :);
            YTrainR = lbr(trainInds);
            YTestR = lbr(testInds);

            % Regular model  (you can turn Optimize off)
            %Mdl = fitctree(XTrain, YTrain, 'OptimizeHyperparameters','auto');
            Mdl = fitctree(XTrain, YTrain, 'CrossVal', 'on', ClassNames=classnames);
            % Random model (you can turn Optimize off)
            MdlR = fitctree(XTrain, YTrainR, 'CrossVal', 'on', ClassNames=classnames);
            num_samp = size(XTrain, 2);

            TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Data = Mdl;
            TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Random = MdlR;
            TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).AccessoriesTrain =struct('XTrain', XTrain, 'YTrain', YTrain, 'XTrainR', XTrainR, 'YTrainR', YTrainR);
            TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).AccessoriesTest= struct('XTest', XTest, 'YTest', YTest, 'YTestR', YTestR);

            clear Mdl MdlR
            progressbar([], K/maxK);
    end

    progressbar(SI/num_sessions, []);
end
toc


%% Extract scores and depth for each session/K 

for SI = 1:num_sessions
    sesh_id = sessions(SI)
    for K = 1:maxK
        if K > length(fieldnames(TreeResults_CV.(['S_' num2str(sesh_id)])))  %If this K iteration  doesn't exist, continue, then eventually loops through SI
            continue
        end
        for i_cv = 1:10 %Cross validated 10 folds
            XTest =  TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).AccessoriesTest.XTest;
            Mdl = TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Data.Trained{i_cv};
            MdlR = TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Random.Trained{i_cv};

            [predlabel,score1,~] = predict(Mdl, XTest);
            [predlabelR,score2,~] = predict(MdlR, XTest);


            YTest = {TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).AccessoriesTest.YTest}';
            YTestR =  {TreeResults_CV.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).AccessoriesTest.YTestR}';

            scores(i_cv,1) = 1 - sum(cellfun(@isequal, YTest, predlabel)) / length(YTest);
            scoresR(i_cv,1) = 1 - sum(cellfun(@isequal, YTestR, predlabelR)) / length(YTest);
            tree_depth(i_cv,1) = getDepth(Mdl);
            tree_depthR(i_cv,1) = getDepth(MdlR);

            clear XTest Mdl MdlR predlabel predlabelR Ytest YTestR
        end
        %Save the distributions for each SI/K, this should be much
        %more manageable
        TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Score = scores;
        TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).ScoreR = scoresR;
        TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Depth = tree_depth;
        TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).DepthR = tree_depthR;
    end
    clear scores scoresR tree_depth tree_depthR
end
save([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat'], 'TreeResults_CV_Compact')

%Keep one value of score and depth per SIxK
for SI = 1:num_sessions
    sesh_id = sessions(SI);
    for K = 1:maxK
        if K > length(fieldnames(TreeResults_CV.(['S_' num2str(sesh_id)])))  %If this K iteration  doesn't exist, continue, then eventually loops through SI
            score(sesh_id,K) = NaN;
            scoreR(sesh_id,K) = NaN;
            depth(sesh_id, K) = NaN;
            depthR(sesh_id, K) = NaN;
            continue
        end

        score(sesh_id,K) = median(TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Score);
        scoreR(sesh_id,K) = median(TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).ScoreR);
        depth(sesh_id, K) = median(TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).Depth);
        depthR(sesh_id, K) = median(TreeResults_CV_Compact.(['S_' num2str(sesh_id)]).(['K' num2str(K)]).DepthR);
    end
end
%Append it to TreeResults_CV_Compact
TreeResults_CV_Compact.Scores = score;
TreeResults_CV_Compact.ScoresR = scoreR;
TreeResults_CV_Compact.Depths = depth;
TreeResults_CV_Compact.DepthsR =depthR;
save([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat'], 'TreeResults_CV_Compact')
