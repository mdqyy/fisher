function [ scores ] = fisher_approx_detector( filename, sift, frames, im_size, cell_size1, cell_size2, w, spm, boxes )
%FISHER_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

frames = frames(1:2, :);

if isempty(cell_size2)
    scores = mexFisherApproxDetectorEq(sift, int32(frames), int32(im_size), int32(cell_size1), int32(spm), w, int32(boxes), filename);
else
    scores = mexFisherApproxDetector(sift, int32(frames), int32(im_size), int32(cell_size1), int32(spm), w, int32(boxes), filename, int32(cell_size2));
end

end

