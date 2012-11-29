function [ scores ] = fisher_net_deepspm_classifier( filename, sift, frames, im_size, spm, w )
%FISHER_NET_DEEPSPM_CLASSIFIER Summary of this function goes here
%   Detailed explanation goes here

height = im_size(1);
width = im_size(2);
frames = frames(1:2, :);

divs_row = [];
divs_col = [];
for i = 1 : size(spm, 1)
    divs_row = unique([ceil(linspace(1, height+1, spm(i, 1)+1)), divs_row]);
    divs_col = unique([ceil(linspace(1, width+1, spm(i, 2)+1)), divs_col]);
end

npyrs = sum(prod(spm, 2)) + 1;
[~, fnet, fnorm, nodes] = fisher_net(filename, sift, frames, im_size, {divs_row, divs_col}, w, spm);

nnodes = numel(nodes);
nmodels = size(w, 2);
score = zeros(npyrs, nmodels);

% 1*1 pyr
gbox = [1, 1, width, height];
mask = computeIrregularGridMask(height, width, divs_row, divs_col, gbox);
fn = reshape(fnet(:, 1, :), nnodes, nmodels);
score(1, :) = computeFisherScore(fn, fnorm, mask);

ss = 2;
for i = 1:size(spm, 1)
   pyr = spm(i, :);
   stepx = (gbox(3) - gbox(1) + 1) / pyr(2);
   stepy = (gbox(4) - gbox(2) + 1) / pyr(1);

   for p = 1:pyr(2)
       for q = 1:pyr(1)
           gbox_pyr(1) = ceil(gbox(1) + (p-1)*stepx);
           gbox_pyr(2) = ceil(gbox(2) + (q-1)*stepy);
           gbox_pyr(3) = ceil(gbox(1) + p*stepx - 1);
           gbox_pyr(4) = ceil(gbox(2) + q*stepy - 1);

           gbox_pyr(3) = min([gbox(3), gbox_pyr(3)]);
           gbox_pyr(4) = min([gbox(4), gbox_pyr(4)]);

           mask_pyr = computeIrregularGridMask(height, width, divs_row, divs_col, gbox_pyr);
           fn_pyr = reshape(fnet(:, ss, :), nnodes, nmodels);
           score(ss, :) = computeFisherScore(fn_pyr, fnorm, mask_pyr);
           ss = ss + 1;
       end
   end

end

scores = sum(score, 1) ./ sqrt(npyrs);

end
