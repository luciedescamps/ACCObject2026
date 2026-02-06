%After adding TRACKER to your path, run this to 
% apply the FOV registrations transformations to the spatial weights 
% and then compute the displacements.


for si = 2:size(outputs,2)
    S = full(outputs{si}.spatial_weights);
    ops = W_reg{si}.operations{1};
    S_out = imwarp(S, ops.inref, ops.tform, 'OutputView', ops.outref);
    outputs{si}.transformed_spatial_weights= ndSparse(S_out);
end
