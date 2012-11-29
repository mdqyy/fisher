function [ mask, imask, lc ] = computeGridMask(height, width, cell_size, roi)

cell_height = cell_size(1) ;
cell_width = cell_size(2) ;

rnd = ceil(height / cell_height) ;
cnd = ceil(width / cell_width) ;

cellsx = 1:cell_width:width+1 ;
cellsy = 1:cell_height:height+1 ;
binsx = vl_binsearch(cellsx, roi([1 3])) ;
binsy = vl_binsearch(cellsy, roi([2 4])) ;

mask = zeros(rnd, cnd) ;
mask(binsy(1):binsy(2), binsx(1):binsx(2)) = 1 ;
lc = [binsx(1), binsy(1), binsx(2), binsy(2)] ;

% imask
cell_area = cell_width * cell_height ;
imask = mask ;

% within single cell
if binsx(1) == binsx(2)
    x1 = roi(3) - roi(1) + 1 ;
    x2 = x1 ;
else
    x1 = cellsx(binsx(1)) + cell_width - roi(1) ;
    x2 = roi(3) - cellsx(binsx(2)) + 1 ;
end
if binsy(1) == binsy(2)
    y1 = roi(4) - roi(2) + 1 ;
    y2 = y1 ;
else
    y1 = cellsy(binsy(1)) + cell_height - roi(2) ;
    y2 = roi(4) - cellsy(binsy(2)) + 1 ;
end

imask(binsy(1), binsx(1)) = x1 * y1 / cell_area ;
imask(binsy(1), binsx(2)) = x2 * y1 / cell_area ;
imask(binsy(2), binsx(1)) = x1 * y2 / cell_area ;
imask(binsy(2), binsx(2)) = x2 * y2 / cell_area ;

%
if binsx(2)-binsx(1) > 1
    imask(binsy(1), binsx(1)+1 : binsx(2)-1) = y1 / cell_height ;
    imask(binsy(2), binsx(1)+1 : binsx(2)-1) = y2 / cell_height ;
end
if binsy(2)-binsy(1) > 1
    imask(binsy(1)+1 : binsy(2)-1, binsx(1)) = x1 / cell_width ;
    imask(binsy(1)+1 : binsy(2)-1, binsx(2)) = x2 / cell_width ;
end

%
mask = reshape(mask, 1, numel(mask)) ;
imask = reshape(imask, 1, numel(imask)) ;

end