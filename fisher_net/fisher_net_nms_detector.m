function [ scores, sub_scores, ovlps ] = fisher_net_nms_detector( filename, sift, frames, im_size, cell_size, w, boxes, spm )
%FISHER_NET_APPROX_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

height = im_size(1);
width = im_size(2);
cell_height = cell_size(1);
cell_width = cell_size(2);
frames = frames(1:2, :);

if nargin < 8 || isempty(spm)
    npyrs = 1;
    [~, fnet, fnorm, nodes] = fisher_net_norm(filename, sift, frames, [height, width], [cell_height, cell_width], w);
else
    npyrs = sum(prod(spm, 2)) + 1;
    [~, fnet, fnorm, nodes] = fisher_net_norm(filename, sift, frames, [height, width], [cell_height, cell_width], w, spm); 
end

nnodes = numel(nodes);
nmodels = size(w, 2);
nboxes = size(boxes, 1);
scores = zeros(nboxes, nmodels);

% wiht or without spatial pyramids
if npyrs > 1
    
    for i = 1:nboxes
        return;
    end

else
    
    for i = 1:nboxes
        box = boxes(i, :);
        [masks, ~, maskboxes] = computeGridMultiMasks(height, width, [cell_height cell_width], box);
        [scores(i, :), sub_scores, ovlps] = computeNMSFisherScore(fnet, fnorm, box, maskboxes, masks);
    end
    
end

end

