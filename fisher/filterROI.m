function [ si, fr ] = filterROI( sift, frames, box )
%filterROI Summary of this function goes here
%   Detailed explanation goes here

idx = frames(1, :) >= box(1) & frames(1, :) <= box(3)...
    & frames(2, :) >= box(2) & frames(2, :) <= box(4);
si = sift(:, idx);
fr = frames(:, idx);

%from absolute to relative location
fr(1, :) = fr(1, :) - box(1) + 1;
fr(2, :) = fr(2, :) - box(2) + 1;

end
