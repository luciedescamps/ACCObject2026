%%This code is to reproduce panels from Corder and Ahanonu et al, 2019 DOI: 10.1126/science.aap8586
% It uses the centroids to compute distances between all cells and matched
% cells to perform a post-hoc quality evaluation of the registration
% procedure. The centroids are weighted to their mean as described in
% Corder and Ahanonu et al. Before this you need to run the
% transform_centroids.m script. 

pixelSize = 4.0; % Adjust to your actual micron/pixel ratio
mice = {'Blossom', 'Bubbles', 'Buttercup', 'Daisy', 'Oksana', 'Poppy', 'Villanelle'};
numMice = size(mice,2); 

% Master storage
pooled_within_session_dist = []; 
pooled_matched_pairwise_dist = []; 
pooled_errorX = []; 
pooled_errorY = []; 

for mouseIdx = 1:numMice

    mouse = mice{mouseIdx};

    centroids_path = ['W:\luciede\tracker\july22\Iterative_Cell_Registration-FastAligner\centroids_stats_' mouse '.mat']; %Fetching the centroids after re-running the transformation
    load(centroids_path)

    curr_outputs = outputs;
    curr_global_map =global_map;
    
    % Define session subset based on total sessions available
    totalSessions = length(curr_outputs);
    if totalSessions == 10
        targetSessions = [3, 5:10]; 
    else
        targetSessions = [4, 6:11]; 
    end
    
    fprintf('Mouse %d: Using sessions [%s].', mouseIdx, num2str(targetSessions));

    %% Weighted Centroid Extraction (cf Corder & Ahanonu 2019)
    mouse_centroids = cell(max(targetSessions), 1);
    for si = targetSessions
        if si == 1
            spatial = curr_outputs{si}.spatial_weights;
        else
            spatial = curr_outputs{si}.transformed_spatial_weights;
        end
        
        numCells = size(spatial, 3);
        tmpCents = nan(numCells, 2);
        
        for c = 1:numCells
            spatial_filter = full(spatial(:,:,c)); 
            if any(spatial_filter(:) > 0)
                % binarize at 50% max value
                max_val = max(spatial_filter(:));
                binary_mask = spatial_filter >= (0.5 * max_val);

                %keep pixels connected to max value pixel
                [~, max_idx] = max(spatial_filter(:));
                cc = bwconncomp(binary_mask);
                keep_mask = false(size(binary_mask));
                for i = 1:cc.NumObjects
                    if any(cc.PixelIdxList{i} == max_idx)
                        keep_mask(cc.PixelIdxList{i}) = true;
                        break;
                    end
                end

                %weighted arithmetic mean
                [rows, cols] = find(keep_mask);
                weights = spatial_filter(keep_mask);
                sum_weights = sum(weights);
                centroid_x = sum(cols .* weights) / sum_weights;
                centroid_y = sum(rows .* weights) / sum_weights;
                tmpCents(c, :) = [centroid_x, centroid_y] * pixelSize;
            end
        end
        mouse_centroids{si} = tmpCents;
    end
    
    %% Pairwise Within-session Distances
    for si = targetSessions
        session_cents = mouse_centroids{si};
        session_cents(any(isnan(session_cents), 2), :) = []; 
        if size(session_cents, 1) > 1
            pooled_within_session_dist = [pooled_within_session_dist, pdist(session_cents)];
        end
    end
    
    %% for cells matched in at least 2 sessions
    for g = 1:size(curr_global_map, 1)
        % Identify which target sessions this global cell appears in
        active_target_sessions = intersect(find(curr_global_map(g, :) > 0), targetSessions);
        
        % Minimum 2 sessions requirement
        if length(active_target_sessions) >= 2
            g_cents = [];
            for s_idx = active_target_sessions
                c_idx = curr_global_map(g, s_idx);
                g_cents = [g_cents; mouse_centroids{s_idx}(c_idx, :)];
            end
            
            %pairwise distances for all matches of the same cell
            pooled_matched_pairwise_dist = [pooled_matched_pairwise_dist, pdist(g_cents)];
            
            %calculate offsets from the mean location (alignment error)
            mean_cent = mean(g_cents, 1);
            offsets = g_cents - mean_cent;
            pooled_errorX = [pooled_errorX; offsets(:,1)];
            pooled_errorY = [pooled_errorY; offsets(:,2)];
        end
    end
end

%% Save Variables
save('Pooled_Data_Min2Sessions.mat', 'pooled_within_session_dist', ...
     'pooled_matched_pairwise_dist', 'pooled_errorX', 'pooled_errorY');

%% Viz for pooled data

%S1B
figure
[fD, xD] = ecdf(pooled_within_session_dist);
plot(xD, fD, 'k', 'LineWidth', 2); hold on;
p01 = prctile(pooled_within_session_dist, 0.01);
xline(p01,'-r',{'0.01th percentile'});
xlabel('Within session pairwise cell-cell distances (\mu m)'); ylabel('Cumulative proportion'); title('Panel D');

%S1C
figure
binWidth = 0.5; % Finer bins to match paper style
edges = 0:binWidth:10;
countsD = histcounts(pooled_within_session_dist, edges);
plot(edges(1:end-1) + binWidth/2, countsD, 'k', 'LineWidth', 1.5); hold on;
%xline(5, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Grey threshold line
xline(p01, 'r', 'LineWidth', 1.5); % Red 0.01th percentile line
xlabel('Pairwise cell-cell distances (\mu m)');
ylabel('Number of cell-cell pairs'); title('D Inset (Zoom)');
xlim([0 10]); xticks(0:2.5:10);
box off


% S1D
figure
[fE, xE] = ecdf(pooled_matched_pairwise_dist);
plot(xE, fE, 'k', 'LineWidth', 2); hold on;
xline(p01, 'r', 'LineWidth', 1.5); % Red line consistent with Inset
xlim([0 20]); xlabel('Matched cell-cell pairwise distances (\mu m)'); ylabel('Cumulative proportion'); title('Panel E');

% S1E
figure
% Filter coordinates to strictly +/- 10um to remove accumulation frame
inB = abs(pooled_errorX) < 10 & abs(pooled_errorY) < 10;
[countsF, centersF] = hist3([pooled_errorX(inB), pooled_errorY(inB)], 'Ctrs', {-10:0.1:10, -10:0.1:10});
%[countsF, centersF] = hist3([pooled_errorX, pooled_errorY], 'Ctrs', {-8:0.1:8, -8:0.1:8});
pctF = (countsF / length(pooled_errorX)) * 100;
pctF(pctF == 0) = NaN; % Background white space

h = imagesc(centersF{1}, centersF{2}, pctF'); set(h, 'AlphaData', ~isnan(pctF'));
colormap(turbo); cb = colorbar; ylabel(cb, 'Percentage (%)');
clim([0 0.6]); axis image; set(gca, 'YDir', 'normal');
xlabel('X (\mu m)'); ylabel('Y (\mu m)'); title('Alignment Errors');