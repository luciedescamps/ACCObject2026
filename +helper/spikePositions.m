% Get spike positions in 2D
%
% This function converts spike timestamps to x/y-coordinates.
%
%  USAGE
%   [spkPos, spkInd, rejected, posCounts] = data.getSpikePositions(spikeTs, pos, <posCol>)
%   spikeTs     Vector of spike timestamps, size Sx1.
%   pos         NxK matrix of position samples in form [t x y ...]. Can have 2 or more columns.
%               The first column is the most important as it must contain time samples.
%   posCol      Optional positions validity columns, vector 1xI. If provided, contains indices of
%               columns that are checked for having NaN values when a decision about rejecting a spike
%               is made. If omitted, then all non-time columns of pos are used.
%               See also `rejected` output argument.
%
%   spkPos      Matrix MxK of spike positions. 1 column is original spike time, the rest are
%               the non-time columns from pos matrix. M is the number of valid spikes, can be less than S.
%   spkInd      Linear indices of spikes in pos matrix, size Mx1.
%   rejected    Indices of rejected spikes. These indices correspond to input vector spikeTs.
%               length(rejected) + length(spkInd) = length(spikeTs);
%               Spike is rejected if any (posCol) of non-time pos columns for that spike is NaN,
%               i.e. any(isnan(pos(spikeIndex, 2:end))).
%   posCounts   Vector Nx1. Each element indicates how many spikes correspond to position sample
%               given by element index. posCounts(3) = 2 indicates that 2 spikes correspond to position
%               sample 2. And posCounts(5) = 0 indicates that no spikes correspond to position sample 5.
%               If every spike falls to a distinct position sample, then this vector will consists of [0 1] only.
%
%  Copied from BNT data.getSpikePositions

function [spkPos, spkInd, rejected, posCounts] = spikePositions(spikeTs, pos, validityColumns)
    if nargin < 1
        error('Incorrect number of input arguments (type ''help <a href="matlab:help data.getSpikePositions">data.getSpikePositions</a>'' for details).');
    end

    if isempty(spikeTs) || isempty(pos)
        spkPos = [];
        spkInd = [];
        rejected = [];
        posCounts = [];
        return;
    end
    if nargin < 3
        validityColumns = 2:size(pos, 2);
    else
        if nanmax(validityColumns) > size(pos, 2)
            error('BNT:arg', 'Argument ''posCol'' is invalid. It contains reference to column that is greater than number of columns in ''pos'' matrix.');
        end
    end

    post = pos(:, bntConstants.PosT);
    posValues = pos;
    % remove time and leave only coordinates or whatever is present in pos
    posValues(:, bntConstants.PosT) = [];
    posCounts = zeros(size(post));

    sampleTime = mean(diff(post));
    if isnan(sampleTime)
        sampleTime = data.sampleTime('sec');
    end
    minTime = min(post);
    maxTime = max(post) + sampleTime;
    rejected = [];
    numRejected = 0;

    % make sure spikeTs is a column vector
    spikeTs = spikeTs(:);

    N = size(spikeTs, 1);
    spkPos = zeros(N, size(pos, 2));
    spkInd = zeros(N, 1);

    count = 0;
    tmpSpkInd = knnsearch(post, spikeTs); % another option:
                                          % >> tri = delaunayn(post');
                                          % >> dsearchn(post', tri, spikeTs)

    for i = 1:length(tmpSpkInd)
        ind = tmpSpkInd(i);
        if spikeTs(i) < minTime || spikeTs(i) > maxTime
            numRejected = numRejected + 1;
            rejected(numRejected, 1) = i; %#ok<*AGROW>
            continue;
        end
        if any(isnan(pos(ind, validityColumns)))
            numRejected = numRejected + 1;
            rejected(numRejected, 1) = i;
        else
            count = count + 1;
            spkPos(count, 1) = spikeTs(i);
            spkPos(count, 2:end) = posValues(ind, :);
            spkInd(count) = ind;
            posCounts(ind) = posCounts(ind) + 1;
        end
    end

    spkPos = spkPos(1:count, :);
    spkInd = spkInd(1:count);
end
