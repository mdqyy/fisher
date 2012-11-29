function [ scores ] = fisher_net_detector( filename, sift, frames, im_size, cell_size, w, boxes, spm )
%FISHER_NET_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

height = im_size(1);
width = im_size(2);
cell_height = cell_size(1);
cell_width = cell_size(2);
frames = frames(1:2, :);

if nargin < 8 || isempty(spm)
    [~, fnet, fnorm, nodes] = fisher_net_norm(filename, sift, frames, [height, width], [cell_height, cell_width], w);
else
    [~, fnet, fnorm, nodes] = fisher_net_norm(filename, sift, frames, [height, width], [cell_height, cell_width], w, spm); 
end

nnodes = numel(nodes);
nmodels = size(w, 2);
npyrs = sum(prod(spm, 2)) + 1;
nboxes = size(boxes, 1);
scores = zeros(nboxes, nmodels);

% wiht or without spatial pyramids
if npyrs > 1
    
    for i = 1:nboxes
        
        box = boxes(i, :);
        score = zeros(npyrs, nmodels);
        
        % 1*1 pyr
        mask = computeGridMask(height, width, [cell_height cell_width], box);
        fn = reshape(fnet(:, 1, :), nnodes, nmodels);
        score(1, :) = computeFisherScore(fn, fnorm, mask);

        ss = 2;  
        for j = 1:size(spm, 1)
           pyr = spm(j, :);
           stepx = (box(3) - box(1) + 1) / pyr(2);
           stepy = (box(4) - box(2) + 1) / pyr(1);

           for p = 1:pyr(2) 
               for q = 1:pyr(1)
                   %box_pyr(1) = ceil(box(1) + (p-1)*stepx);
                   %box_pyr(2) = ceil(box(2) + (q-1)*stepy);
                   %box_pyr(3) = floor(box(1) + p*stepx);
                   %box_pyr(4) = floor(box(2) + q*stepy);       
                   box_pyr(1) = box(1) + (p-1)*stepx;
                   box_pyr(2) = box(2) + (q-1)*stepy;
                   box_pyr(3) = box(1) + p*stepx - 1;
                   box_pyr(4) = box(2) + q*stepy - 1;

                   box_pyr(3) = min([box(3), box_pyr(3)]);
                   box_pyr(4) = min([box(4), box_pyr(4)]);

                   mask_pyr = computeGridMask(height, width, [cell_height cell_width], box_pyr);
                   fn_pyr = reshape(fnet(:, ss, :), nnodes, nmodels);
                   score(ss, :) = computeFisherScore(fn_pyr, fnorm, mask_pyr);
                   ss = ss + 1;
               end
           end

        end
        
        %scores(i, :) = sum(score, 1);
        scores(i, :) = sum(score, 1) ./ sqrt(npyrs);
    end

else
    
    for i = 1:nboxes
        box = boxes(i, :);
        mask = computeGridMask(height, width, [cell_height cell_width], box);
        scores(i, :) = computeFisherScore(fnet, fnorm, mask);
    end
    
end

end

