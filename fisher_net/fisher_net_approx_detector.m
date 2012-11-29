function [ scores, sub_scores ] = fisher_net_approx_detector( filename, sift, frames, im_size, cell_size1, cell_size2, w, boxes, spm )
%FISHER_NET_APPROX_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

height = im_size(1);
width = im_size(2);
frames = frames(1:2, :);

if nargin < 9 || isempty(spm)
    npyrs = 1;
    [fish, fnet, ~, nodes] = fisher_net_nono(filename, sift, frames, [height, width], cell_size1, w);
else
    npyrs = sum(prod(spm, 2)) + 1;
    [fish, fnet, ~, nodes] = fisher_net_nono(filename, sift, frames, [height, width], cell_size1, w, spm);
end

if ~isempty(cell_size2)
    if sum(mod(cell_size2, cell_size1)) == 0
        aggreg_factor = cell_size2 ./ cell_size1;
        fish = aggreg_fisher(fish, size(nodes), aggreg_factor);
    else
        [fish] = fisher_net_nono(filename, sift, frames, [height, width], cell_size2);
    end
else
    cell_size2 = cell_size1;
end

%calc norm matrix
fnorm = fish' * fish;

nnodes = numel(nodes);
nmodels = size(w, 2);
nboxes = size(boxes, 1);
scores = zeros(nboxes, nmodels);
sub_scores = zeros(nboxes, npyrs, nmodels);

% with or without spatial pyramids
if npyrs > 1
    
    for i = 1:nboxes
        
        box = boxes(i, :);
        score = zeros(npyrs, nmodels);
        
        % 1*1 pyr
        [~, imask1] = computeGridMask(height, width, cell_size1, box);
        [~, imask2] = computeGridMask(height, width, cell_size2, box);
        fn = reshape(fnet(:, 1, :), nnodes, nmodels);
        score(1, :) = computeFisherScoreIMask(fn, fnorm, imask1, imask2);
        
        ss = 2;  
        for j = 1:size(spm, 1)
           pyr = spm(j, :);
           stepx = (box(3) - box(1) + 1) / pyr(2);
           stepy = (box(4) - box(2) + 1) / pyr(1);

           for p = 1:pyr(2)
               for q = 1:pyr(1)
                   box_pyr(1) = ceil(box(1) + (p-1)*stepx);
                   box_pyr(2) = ceil(box(2) + (q-1)*stepy);
                   box_pyr(3) = ceil(box(1) + p*stepx - 1);
                   box_pyr(4) = ceil(box(2) + q*stepy - 1);

                   box_pyr(3) = min([box(3), box_pyr(3)]);
                   box_pyr(4) = min([box(4), box_pyr(4)]);

                   [~, imask1_pyr] = computeGridMask(height, width, cell_size1, box_pyr);
                   [~, imask2_pyr] = computeGridMask(height, width, cell_size2, box_pyr);
                   fn_pyr = reshape(fnet(:, ss, :), nnodes, nmodels);
                   score(ss, :) = computeFisherScoreIMask(fn_pyr, fnorm, imask1_pyr, imask2_pyr);
                   ss = ss + 1;
               end
           end
           
        end
        
        %scores(i, :) = sum(score, 1);
        sub_scores(i, :, :) = score ./ sqrt(npyrs);
        scores(i, :) = sum(score, 1) ./ sqrt(npyrs);
    end

else
    
    for i = 1:nboxes
        box = boxes(i, :);
        
        %[masks, ~, maskboxes] = computeGridMultiMasks(height, width, [cell_height cell_width], box);
        %[masks, ~, maskboxes] = computeGridMultiMasks_haha(height, width, [cell_height cell_width], box);
        %scores(i, :) = computeApproxFisherScore(fnet, fnorm, box, maskboxes, masks);

        [~, imask1] = computeGridMask(height, width, cell_size1, box);
        [~, imask2] = computeGridMask(height, width, cell_size2, box);
        scores(i, :) = computeFisherScoreIMask(fnet, fnorm, imask1, imask2);
    end
    
end

end

