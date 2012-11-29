function olap = boxOverlap(allbox1, allbox2)

olap = zeros(size(allbox1, 2), size(allbox2, 2)) ;
for i = 1 : size(allbox1, 2)
    for j = 1 : size(allbox2, 2)
        box1 = vec(allbox1(:, i)) ;
        box2 = vec(allbox2(:, j)) ;

        area1 = (box1(3) - box1(1) + 1) * (box1(4) - box1(2) + 1) ;
        area2 = (box2(3) - box2(1) + 1) * (box2(4) - box2(2) + 1) ;

        box1_ = [box1(1), box1(2), box1(3) - box1(1) + 1, box1(4) - box1(2) + 1]' ;
        box2_ = [box2(1), box2(2), box2(3) - box2(1) + 1, box2(4) - box2(2) + 1]' ;

        recintsct = rectint(box1_', box2_') ;
        recunion = area1 + area2 - recintsct ;
        olap(i, j) = recintsct ./ recunion ;
    end
end
