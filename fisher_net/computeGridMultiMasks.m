function [ masks, lcs, coords ] = computeGridMultiMasks(height, width, cell_size, roi)

cell_height = cell_size(1) ;
cell_width = cell_size(2) ;

rnd = ceil(height / cell_height) ;
cnd = ceil(width / cell_width) ;

cellsx = 1:cell_width:width+1 ;
cellsy = 1:cell_height:height+1 ;
binsx = vl_binsearch(cellsx, roi([1 3])) ;
binsy = vl_binsearch(cellsy, roi([2 4])) ;

near_neighbors = zeros(1, 4) ;
near_neighbors(1) = ~(ismember(roi(1), cellsx) | roi(1)+cell_width > width+1);
near_neighbors(2) = ~(ismember(roi(2), cellsy) | roi(2)+cell_height > height+1);
near_neighbors(3) = ~(ismember(roi(3), cellsx-1) | roi(3)-cell_width < 0);
near_neighbors(4) = ~(ismember(roi(4), cellsy-1) | roi(4)-cell_height < 0);

nmasks = prod(near_neighbors+1, 2);
masks = zeros(nmasks, rnd, cnd) ;
lcs = zeros(nmasks, 4) ;
ss = 1 ;
for x1 = 0 : near_neighbors(1)
    for x2 = 0 : near_neighbors(3)
        for y1 = 0 : near_neighbors(2)
            for y2 = 0 : near_neighbors(4)
                masks(ss, binsy(1)+y1:binsy(2)-y2, binsx(1)+x1:binsx(2)-x2) = 1 ;
                lcs(ss, :) = [binsx(1)+x1, binsy(1)+y1, binsx(2)-x2, binsy(2)-y2] ;
                ss = ss + 1;
            end
        end
    end
end

coords = zeros(nmasks, 4) ;
coords(:, 1) = cellsx(lcs(:, 1))' ;
coords(:, 2) = cellsy(lcs(:, 2))' ;
% coords(:, 3) = cellsx(lcs(:, 3))' + cell_width - 1 ;
% coords(:, 4) = cellsy(lcs(:, 4))' + cell_height - 1 ;
coords(:, 3) = min(cellsx(lcs(:, 3))'+cell_width-1, width) ;
coords(:, 4) = min(cellsy(lcs(:, 4))'+cell_height-1, height) ;

masks = reshape(masks, nmasks, rnd*cnd) ;

end

