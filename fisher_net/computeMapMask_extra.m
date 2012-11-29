function [ mask ] = computeMapMask_extra(height, width, confs, roi, type)

if nargin < 5
    type = 1 ;
end

node_height = confs(1) ;
node_width = confs(2) ;

if type == 1

    rnd = ceil(height / node_height) ;
    cnd = ceil(width / node_width) ;

    mask = zeros(rnd, cnd) ;

    binsx = vl_binsearch(1:node_width:width+1, roi([1 3])) ;
    binsy = vl_binsearch(1:node_height:height+1, roi([2 4])) ;

    mask(binsy(1):binsy(2), binsx(1):binsx(2)) = 1;

elseif type == 2
    rnd = ceil(height / node_height) ;
    cnd = ceil(width / node_width) ;

    mask = zeros(rnd, cnd) ;
    
    linx1 = (1 - node_width/2) : node_width : width+node_width/2 + 1 ;
    linx2 = (1 - node_width/2) : node_width : width+node_width/2 + 1 ;
    liny1 = (1 - node_height/2) : node_height : height+node_height/2 + 1 ;
    liny2 = (1 - node_height/2) : node_height : height+node_height/2 + 1 ;

    binsx1 = vl_binsearch((1 - node_width/2) : node_width : width+node_width/2 + 1, roi(1)) ;
    binsx2 = vl_binsearch((1 - node_width/2) : node_width : width+node_width/2 + 1, roi(3)) ;
    binsy1 = vl_binsearch((1 - node_height/2) : node_height : height+node_height/2 + 1, roi(2)) ;
    binsy2 = vl_binsearch((1 - node_height/2) : node_height : height+node_height/2 + 1, roi(4)) ;
    
    binsx1 = binsx1 - double(mod(roi(1), linx1(binsx1)) == 0) ;
    binsx2 = binsx2 - double(mod(roi(3), linx2(binsx2)) == 0) ;
    binsy1 = binsy1 - double(mod(roi(2), liny1(binsy1)) == 0) ;
    binsy2 = binsy2 - double(mod(roi(4), liny2(binsy2)) == 0) ;
    
    binsx = [binsx1, binsx2] ;
    binsy = [binsy1, binsy2] ;

    mask(binsy(1):binsy(2), binsx(1):binsx(2)) = 1;
    
end