%This is to reproduce panel S1F. It computes the distances between centroids for neighbouring cells and matched cells,
% using a function from TRACKER, then plots the resulting distribution. The
% input are the centroids after applying the same FOV registration as
% TRACKER using transform_centroids

mice = {'Blossom', 'Bubbles', 'Buttercup', 'Daisy', 'Oksana', 'Poppy', 'Villanelle'};
numMice = length(mice);

pooled_matched_dist = [];
pooled_nn_dist = [];
pooled_dprime_log = [];

for mouseIdx = 1:numMice

    mouse = mice{mouseIdx};

    centroids_path = ['W:\luciede\tracker\july22\Iterative_Cell_Registration-FastAligner\centroids_stats_' mouse '.mat'];
    load(centroids_path)

    fprintf('Processing Mouse %d of %d...\n', mouseIdx, numMice);

    [matched_dist, nn_dist, dprime_log]  = check_centroid_distances(obj, 4);


    pooled_matched_dist{mouseIdx} = matched_dist;
    pooled_nn_dist{mouseIdx} = nn_dist;
    pooled_dprime_log{mouseIdx} = dprime_log;

    clear matched_dist nn_dist dprime_log
end

dist_per_px = 4; %Change with the actual resolution of your recording

%Matched cells
pooled_matched = vertcat(pooled_matched_dist{:});
pooled_matched  = pooled_matched  .* dist_per_px;

%Neighbouring cells
pooled_nn = vertcat(pooled_nn_dist{:});
pooled_nn = pooled_nn  .* dist_per_px;

length_unit = '\mu m';

figure;
h1 =histogram(pooled_matched, 'Normalization', 'pdf', 'DisplayStyle','stairs');
hold on;
h2 =histogram(pooled_nn, 'Normalization', 'pdf',  'DisplayStyle','stairs');
xlabel(sprintf(' Distance (%s) ', length_unit));
ylabel 'pdf'
title('Matched and NN centroid distances');
legend 'Matched' 'Nearest neighbor' Location best
xlim([0 200])
xlabel('Distance (\mu m)');
ylabel('Probability Density');
legend('Matched cells', 'Nearest neighbor');

lm = log(pooled_matched);
ln = log(pooled_nn);
dprime_log = abs(mean(lm) - mean(ln)) / sqrt( (var(lm) + var(ln))/2 );