function [ scores ] = fisher_detector( filename, sift, frames, im_size, cell_size, w, confs, boxes  )
%FISHER_DETECTOR Summary of this function goes here
%   Detailed explanation goes here

frames = frames(1:2, :);
scores = mexFisherDetector(sift, int32(frames), int32(im_size), int32(cell_size), int32(confs), w, int32(boxes), filename);


end

