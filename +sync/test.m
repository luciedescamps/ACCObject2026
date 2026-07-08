% Test for accuracy of the two SyncAligner implementations, using simulated
% sync trains
clear

% Set some parameters
nPulses = 10000;
intervalRange = [0.3, 1.5];
sigmaDetect = 0.01;
sigmaDrift = 0;
detectRate = 1;
chunkSz = 20;

tOffsets = [0, 555, -253];
alignmentMethods = {"none", ["xcorr", "rg"]};
alignerClasses = ["SimpleAligner", "PatternAligner"];

[pts, simInfo] = sync.PulseTrain.simulate( ...
    'nTrains', 3, ...
    'nPulses', nPulses, ...
    'intervalRange', intervalRange, ...
    'sigmaDetect', sigmaDetect, ...
    'sigmaDrift', sigmaDrift, ...
    'detectRate', detectRate, ...
    'offsets', tOffsets, ...
    'jumpPos', [0, 0, 0], ...
    'jumpIncr', [0, 0, 0]);

trainA = pts(1);
trainB = pts(2);
trainC = pts(3);

tA = trainA.tRise;
tB = trainB.tRise;
tC = trainC.tRise;

iDetect1 = simInfo.tRawInds{1};
iDetect2 = simInfo.tRawInds{2};
[~, indsABTrue] = ismember(iDetect1, iDetect2);

tMax = simInfo.tRawMax;

for a = 1:2
    
    figure,
    alignerClass = alignerClasses(a);
    suptitle(alignerClass);
    
    subplot(2, 2, 1)
    plot(tA, tB, 'k.');
    xlabel('tpulse_a');
    ylabel('tpulse_b');
    title('Pulse times');
    axis square
    
    subplot(2, 2, 2);
    bar(tB-tA-(tOffsets(2)-tOffsets(1)));
    xlabel('pulse index');
    ylabel('error');
    title('Pulse timing errors');
    axis square
    
    vDetectedAll = {};
    
    subplot(2, 1, 2)
    hold on
    xlabels = alignmentMethods{a};
    set(gca, 'ytick', []);
    title('Pulse detection');
    xlabel('Pulse #');
    ylim([-1, 5]);
    pan('xon');
    zoom('xon');
    
    for s = 1:numel(alignmentMethods{a})
        
        aligner = feval("sync." + alignerClass, trainA, trainB);
        if alignerClass == "PatternAligner"
            p = sync.AlignmentParams();
            p.alignmentMethod = alignmentMethods{a}(s);
            p.verbose = true;
            aligner = sync.PatternAligner(trainA, trainB, p);
        end
        
        % Plot pulse detection errors
        indsABSync = aligner.iMatchAToB;
        indsBASync = aligner.iMatchBToA;
        vMatch = indsABTrue==indsABSync;
        
        vDoubleMatchB = ~isnan(indsBASync);
        v = ~isnan(indsBASync);
        vDoubleMatchB(v) = vDoubleMatchB(v) & vMatch(indsBASync(v));
        
        doubleMatchIndsA = indsBASync(vDoubleMatchB);
        doubleMatchIndsB = find(vDoubleMatchB);
        
        indsMismatch = find(doubleMatchIndsB~=indsABTrue(doubleMatchIndsA));
        vFalsePosB = false(aligner.nA, 1);
        vFalsePosB(indsMismatch) = true;
        vFalsePosB(~isnan(indsABSync)) = false;
        vFalsePosA =  ~vMatch & ~isnan(indsABSync);
        vFalsePos = vFalsePosA | vFalsePosB;
        nFalsePos = sum(vFalsePos);
        
        percentDetected = sum(vMatch) ./ numel(vMatch) * 100;
        
        baseVal = (s-1)*3;
        y = zeros(2, aligner.nA) + baseVal;
        
        % Plot matches
        x = repmat(1:aligner.nA, 2, 1);
        v = vMatch;
        y(2, :) = y(2, :) + v'*2;
        plot(x(:, v), y(:, v), 'k', 'lineWidth', 0.5);
        plot(x(:, v), y(2, v), 'ko');
        
        % Plot non-matches
        v = ~vMatch & ~vFalsePos;
        y(2, :) = y(2, :);
        plot(x(:, v), y(2, v), 'k.');
        
        % Plot false-pos
        v = vFalsePos;
        y(2, :) = y(2, :) + v'*2;
        plot(x(:, v), y(:, v), 'r', 'lineWidth', 0.5);
        plot(x(:, v), y(2, v), 'ro');
        
        str = sprintf('%s\n%.1f%%\n FP = %u', xlabels{s}, percentDetected, nFalsePos);
        text(-0.03, (s-1)*0.6 + 0.2, str, 'sc', 'horizontalAlignment', 'right', ...
            'verticalAlignment', 'middle');
        
        vDetectedAll = v;
        aligners{s} = aligner;
        
    end
    
end

%%

errLim = 0.03*[-1 1];

histBins = linspace(errLim(1), errLim(2), 101);

funcs = {@aToB, @bToA};
funStrs = {'A -> B', 'B -> A'};
syncCols = {'k', 'r'};
syncNames = {'RG', 'PyControl'};
interpMethods = {'piecewise', 'regression', 'brokenStick'};
interpMethodStrs = {'Linear piecewise interp', 'Linear regression', 'Broken-stick regression'};

tExtrap = 100;

% For interp, test all three methods:
% 1) linear-piecewise
% 2) linear regression
% 3) broken-stick linear regression
for i = 1:3
    
    figure('units', 'normalized', 'outerPosition', [0.1, 0.1, 0.4, 0.8]);
    suptitle(interpMethodStrs{i});
    
    % do both forward and backward conversion
    for n = 1:2
        
        axErr = subplot(3, 2, n+2); hold(axErr, 'on');
        axErrHist = subplot(3, 2, n+4); hold(axErrHist, 'on');
        xlabel(axErrHist, 'Error / s');
        ylabel(axErrHist, 'Count');
        
        % Get t1 and t2 with only drift noise (no detection noise)
        i1 = iDetect1;
        i2 = iDetect2;
        t1 = simInfo.tRaw(i1) + simInfo.errDrift{1} + tOffsets(1);
        t2 = simInfo.tRaw(i2) + simInfo.errDrift{2} + tOffsets(2);
        [v, indsABTrue] = ismember(iDetect1, iDetect2);
        t2 = t2(indsABTrue(v));
        t1 = t1(v);
        if n==2
            [t1, t2] = deal(t2, t1);
        end
        
        tLims = [t1(1)-tExtrap, t1(end)+tExtrap];
        tq = linspace(tLims(1), tLims(2), 999)';
        tqIsExtrap = tq < tOffsets(n) | tq > tMax+tOffsets(n);
        t1qValid = tq(~tqIsExtrap);
        
        t2iIdeal = interp1(t1, t2, tq, 'linear', 'extrap'); % query times corrected for drift
        
        for s = 1:2
            
            aligner = aligners{s};
            aligner.extrapLimit = tExtrap;
            
            try
                t2i = funcs{n}(aligner, tq, interpMethods{i});
            catch e
                warning('Interpolation failed: "%s"', e.message);
                continue;
            end
            
            subplot(3, 2, n);
            h(s) = line(tq, t2i, 'color', syncCols{s});
            
            tiErr   = t2i - tq - (tOffsets(2)-tOffsets(1));
            tiErrC  = t2i - t2iIdeal;
            
            bar(axErr, tq, tiErrC, 'faceColor', syncCols{s}, 'faceAlpha', 0.5, ...
                'edgeColor', 'none', 'barWidth', 1);
            histogram(axErrHist, tiErrC(~tqIsExtrap), histBins, 'faceColor', syncCols{s}, 'faceAlpha', 0.5, ...
                'edgeColor', 'none');
            tiAll{s} = t2i;
        end
        
        axVal = subplot(3, 2, n);
        xlabel('t_q', 'interpreter', 'tex');
        ylabel('t_i', 'interpreter', 'tex');
        title(funStrs{n});
        axis('square');
        set(gca, 'ylim', get(gca, 'ylim'));
        
        ax = axErr;
        xlabel(ax, 't_q', 'interpreter', 'tex');
        ylabel(ax, 't_i error', 'interpreter', 'tex');
        axis(ax, 'square');
        ax.YLim = errLim;
        
        x = tExtrap*[0, 0, 1, 1];
        x = {-tExtrap+x, tMax+x};
        axs = [axVal, axErr];
        for a = 1:2
            ax = axs(a);
            y = ax.YLim([1, 2, 2, 1]);
            for n = 1:2
                hp = patch(ax, x{n}, y, 'k', 'faceAlpha', 0.2, 'edgeColor', 'none');
            end
            ax.XLim = [min(x{1}), max(x{2})];
        end
    end
    legend(hp, 'Extrapolated');
    
    try
        leg = legend(h, syncNames);
        title(leg, {'Alignment', 'algorithm'});
    end
    
end

%% Test chaining of SyncAligners using a Chain

% Here we test linking together multiple SyncAligner objects

clear t

chainTrainIds{1} = [1, 2, 1]; % sync A -> B -> A
chainTrainIds{2} = [1, 2, 3, 1]; % sync A -> B -> C -> A
nReps = 2;

for ch = 1:2
    
    trainIds = chainTrainIds{ch};
    
    clear trains
    for s = 1:numel(trainIds)
        trains(s, 1) = pts(trainIds(s));
    end
    trainsRep = repmat(trains(1:end-1), nReps, 1);
    trainsRep(end+1) = trains(end);
    chain = sync.Chain(trains);
    chainRep = sync.Chain(trainsRep);
    
    chainRep.name = sprintf("( %s ) x %u", chain.name, nReps);
    fprintf('Created sync.Chain: "%s"\n', chain.name);
    
    t1Lims = chain.aligners(1).tAMatched([1, end]);
    tq = linspace(t1Lims(1), t1Lims(2), 1000)';
    t2 = chainRep.mapForward(tq);
    tErr = t2-tq;
    
    figure();
    suptitle(chainRep.name);
    
    subplot(1, 2, 1);
    plot(tq, tErr, 'ko');
    xlabel('tq');
    ylabel('Error / s');
    
    subplot(1, 2, 2);
    histogram(tErr);
    xlabel('Error / s');
    ylabel('Count');
    
end

%% Test Network

nTrains = 3;
tOffsets = [0, rand(1, nTrains-1)*10];

[pts, simInfo] = sync.PulseTrain.simulate( ...
    'nTrains', nTrains, ...
    'nPulses', nPulses, ...
    'intervalRange', intervalRange, ...
    'sigmaDetect', 0.0001, ...
    'sigmaDrift', sigmaDrift, ...
    'detectRate', detectRate, ...
    'offsets', tOffsets);

p = sync.AlignmentParams();
p.verbose = true;
syncNet = sync.Network(pts, p);
syncNet.name = sprintf('%s ,', [pts.name]);
fprintf('Created SyncNet: "%s"\n', syncNet.name);

t1Lims = pts(1).tLims + [30, -30];
tq = linspace(t1Lims(1), t1Lims(2), 1000)';

figure();
clear h

for n = 1:nTrains
    
    % Sync via SyncNet
    t2 = syncNet.mapTimes(1, n, tq); % map 1 -> n
    t1 = syncNet.mapTimes(n, 1, tq); % map n -> 1
    tErr1 = t2 - tq - tOffsets(n);
    tErr2 = t1 - tq + tOffsets(n);
    
    % Sync via SyncAligner
    if n > 1
        aligner = sync.PatternAligner(pts(1), pts(n));
        t2A = aligner.aToB(tq);
    else
        t2A = t2;
    end
    tErrA = t2A-t2;
    
    subplot(nTrains, 2, (n-1)*2 + 1);
    hold on
    h(1) = plot(tq, tErr1, 'k.');
    h(2) = plot(tq, tErr2, 'r.');
    h(3) = plot(tq, tErrA, 'b.');
    xlabel('tq');
    ylabel('Error / s');
    str = sprintf('%s <-> %s', pts(1).name, pts(n).name);
    title(str);
    legend(h, {'chain FWD', 'net REV', 'aligner'});
    
    subplot(nTrains, 2, (n-1)*2 + 2);
    histogram(tErr1);
    xlabel('Error / s');
    ylabel('Count');
    title(str);
end

%% Test save/load for each class

pulseTrain = pts(1);
syncAligner = sync.PatternAligner(pts(1), pts(2));
chain = sync.Chain(pts);

objs = {syncNet, chain, syncAligner, pulseTrain};
fn = 'sync_test.mat';
exclClasses = {'containers.Map'};

for o = 1:numel(objs)
    obj = objs{o};
    fprintf('Saving %s obj...', class(obj));
    save(fn, 'obj');
    file = dir(fn);
    fprintf('done. Saved size is %.0f kb.\n', file.bytes/1024);
    tmp = load(fn);
    objLoaded = tmp.obj;
    delete(fn);
    props = sync.getSaveableProperties(obj);
    for p = 1:numel(props)
        prop = props{p};
        try
            val1 = obj.(prop);
            val2 = objLoaded.(prop);
            % Skip
            type = class(val1);
            if ismember(type, exclClasses)
                continue;
            end
        catch
            fprintf('Ignored inaccessible property "%s".\n', prop);
            continue;
        end
        assert(isequal(objLoaded.(prop), obj.(prop)), ...
            'Saved and loaded values for property "%s" are not equal.', prop);
    end
end