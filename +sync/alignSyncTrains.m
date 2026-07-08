function [iAB, iBA, info] = alignSyncTrains(tA, tB, params)
%ALIGNSYNCTRAINS aligns two irregular sync event trains
%
% [INDS, DISTS] = ALIGNSYNCTRAINS(T1, T2), where T1 and T2 are vectors
% containing sync event times, finds the closest match in T2 for every
% event in T1, by using local pattern matching.
%
% [INDS, DISTS] = ALIGNSYNCTRAINS(T1, T2, PLEN), when PLEN is a scalar
% integer, sets the length of the search pattern (default 5).

if nargin < 3 || isempty(params), params = sync.AlignmentParams(); end
P = params;

nA = numel(tA);
nB = numel(tB);

chsz = P.chunkSize;

if P.alignmentMethod == "rg"
    tBChunkMn = movmean(tB, chsz);
    tBChunkMn(1:floor(chsz/2)) = [];
    chunkStartIndsA = (1 : chsz : (nA-chsz+1))';
elseif P.alignmentMethod == "xcorr"
    chunkStartIndsA = (1 : chsz : (nA-chsz))';
    intA = diff(tA);
    intB = diff(tB);
    intA2 = intA.^2;
    intB2 = intB.^2;
    ch1 = ones(chsz, 1);
    xcB2 = xcvalid(intB2, ch1);
end

nChunks = numel(chunkStartIndsA);
chunkStartIndsB = nan(nChunks, 1);
chunkStartIndsB2 = nan(nChunks, 1);
chunkErrs = nan(nChunks, 2);

if P.verbose
    fprintf('matching pulses (n = %u, %u) in %u chunks, method "%s"...', ...
        nA, nB, nChunks, P.alignmentMethod);
    tic();
    fprintf('    ');
end

for c = 1:nChunks
    
    % Extract the local t1 pattern
    chunkStartIdx = chunkStartIndsA(c);
    chunkIndsA = chunkStartIdx + (0:chsz-1);
    
    if P.alignmentMethod == "rg"
        tAChunk = tA(chunkIndsA);
        [inds, errs] = localPatternSearchRg(tAChunk, tB, tBChunkMn);
        chunkErrs(c, :) = errs;
    elseif P.alignmentMethod == "xcorr"
        chA = intA(chunkIndsA);
        sumChA2 = sum(intA2(chunkIndsA));
        err = xcB2 + sumChA2 - 2*xcvalid(intB, chA);
        [chunkErrs(c, :), inds] = mink(err, 2);
    end
    
    chunkStartIndsB(c) = inds(1);
    chunkStartIndsB2(c) = inds(2);
    
    % If final chunk only contains one possible alignment, there is no
    % second-best error.
    if c == nChunks && numel(inds) == 1
        chunkErrs(c, 2) = nan;
    end
    
    if P.verbose
        pcDone0 = floor((c-1)/nChunks * 100);
        pcDone = floor(c/nChunks * 100);
        if pcDone > pcDone0
            fprintf('\b\b\b\b%3u%%', pcDone);
        end
    end
end

validChunks = true(nChunks, 1);

% Truncate pulse trains such that ends match; this option should be used
% when the recorded pulse trains do not fully overlap. Enabling this option
% will truncate unmatched ends of either pulse train, meaning that these
% unmatched pulses will not affect the alignment quality statistics.
if P.truncate
    
    % chA/cbB = chunk in pulseTrain A/B
    diChB = diff(chunkStartIndsB);
    diChB1 = diChB(1:end-1);
    diChB2 = diChB(2:end);
    v = diChB1==chsz & diChB2==chsz;
    
    % Find the first and last consecutive chunk pairs
    iChAFirst = find(v, 1, 'first'); % first valid chunk
    iChALast = find(v, 1, 'last') + 2; % last valid chunk
    
    % Define
    validChunks(1:iChAFirst-1) = false;
    validChunks(iChALast+1:end) = false;
    
end

chunkErrs = chunkErrs ./ chsz;
chunkErrs(chunkErrs<=0) = eps; % avoid -Inf in logMse

if P.verbose, fprintf('... done in %.3f seconds.\n', toc()); end

logMse = log(chunkErrs);
logMseDiff = diff(logMse, [], 2);
if P.minChunkSeparation
    validMseDiff = logMseDiff > P.minChunkSeparation;
    pc = sum(validMseDiff) / nChunks * 100;
    if P.verbose
        fprintf('%.1f%% of chunks pass the separation criterion (%.1f).\n', ...
            pc, P.minChunkSeparation);
    end
    validChunks = validChunks & validMseDiff;
end

logMseValid = logMse(validChunks, :);
if any(validChunks)
    gmm = gmdistribution.fit(logMseValid(:), 2);
    iCluValid = gmm.cluster(logMseValid(:, 1));
    iClu = nan(nChunks, 1);
    iClu(validChunks) = iCluValid;
    [~, imn] = min(gmm.mu);
    validChunks = validChunks & iClu == imn;
    iValid = find(validChunks);
    mu = gmm.mu;
    sigma = gmm.Sigma;
else
    gmm = [];
    iValid = [];
    mu = nan(1, 2);
    sigma = nan(1, 2);
end

iAB = nan(nA, 1);
iBA = nan(nB, 1);

for n = 1:numel(iValid)
    % get inds first pulse in current chunk
    ich = iValid(n);
    iA = chunkStartIndsA(ich); % inds in tA
    iB = chunkStartIndsB(ich); % inds in tB
    
    % calculate indices of all pulses in current chunk
    i1 = iA:(iA+chsz-1); % inds in tA
    i2 = iB:(iB+chsz-1); % inds in tBs
    
    % Check whether the chunk overlaps with any previous chunks. This may
    % happen in rare cases when chunks are poorly matched. Overlapping
    % chunks will ruin the indexing of the matched pulses, so we'll ignore
    % any chunk which overlaps with a previous chunk.
    chunkOverlap = any(~isnan(iAB(i1))) | any(~isnan(iBA(i2)));
    
    if chunkOverlap
        validChunks(n) = false;
        warning("Chunk %u will be discarded because it overlaps with a previous chunk.", n);
    else
        iAB(i1) = i2;
        iBA(i2) = i1;
    end
end

nvAB = sum(~isnan(iAB));
nvBA = sum(~isnan(iBA));
assert(nvAB == nvBA, "Unequal lengths of matched pulse trains (%u vs %u)", nvAB, nvBA)


% Check separation
separation = abs(diff(mu)) / sum(sqrt(sigma));
okSeparation =  separation >= P.minTotalSeparation;

% Check sorted
inds = iBA;
inds(isnan(inds)) = [];
tAMatch = tA(inds);

inds = iAB;
inds(isnan(inds)) = [];
tBMatch = tB(inds);

if ~okSeparation
    warn('alignSyncTrains:badSeparation', 'best-matched pulses are poorly separated', P);
end


if any(validChunks)
    okOrder = nanmin(diff(tAMatch)) > 0 && ...
        nanmin(diff(tBMatch)) > 0;
    if ~okOrder
        warn('alignSyncTrains:badOrder', 'matched pulses are not in ascending temporal order', P);
    end
end

if numel(tAMatch) ~= numel(tBMatch)
    warn('alignSyncTrains:badMatchCount', 'number of matches differs between the two pulse trains', P);
end

info = struct( ...
    'chunkErrors', chunkErrs, ...
    'gmm', gmm, ...
    'chunkIndsA', chunkStartIndsA, ...
    'chunkIndsB', chunkStartIndsB, ...
    'chunkValid', validChunks, ...
    'logMse', mean(logMseValid(:, 1)),  ...
    'separation', separation);

end

function warn(id, msg, P)
warning(id, msg);
end

function [bestInds, bestErrs] = localPatternSearchRg(tAChunk, tB, tBChunkMn)
% Search for the occurrence of a short event pattern in a train.
%
% The index and error of the best and 2nd-best matches are returned.

ntA = numel(tAChunk);
ntB = numel(tB);

tAChunk0 = tAChunk-mean(tAChunk);
nShifts = ntB-ntA+1;
errors = zeros(nShifts, 1);
inds0 = 0:ntA-1;

for sh = 1:nShifts
    tBChunk0 = tB(inds0+sh) - tBChunkMn(sh);
    errors(sh) = sum((tAChunk0-tBChunk0).^2);
end

[bestErrs, bestInds] = mink(errors, 2);

end

function y = xcvalid(x1, x2)
%should replicate the behaviour of 'valid' mode of numpy.correlate
[y, lags] = xcorr(x1, x2);
n1 = numel(x1);
n2 = numel(x2);
idx0 = find(lags==0);
if n1>n2
    idx = idx0 + n1-n2;
    y = y(idx0:idx);
else
    idx = n1;
    y = y(idx:idx0);
end
end