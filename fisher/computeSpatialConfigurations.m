function ixs = computeSpatialConfigurations(frames, height, width, confs)

ixs = zeros(sum(prod(confs, 2)), size(frames, 2), 'single') ;

numSpatialX = confs(:, 2) ;
numSpatialY = confs(:, 1) ;
cnt = 1 ;

for i = 1 : length(numSpatialX)
  binsx = vl_binsearch(linspace(1, width+1, numSpatialX(i)+1), frames(1, :)) ;
  binsy = vl_binsearch(linspace(1, height+1, numSpatialY(i)+1), frames(2, :)) ;

  for p = 1 : numSpatialX(i)
      for q = 1 : numSpatialY(i)
          ixs(cnt, :) = binsx == p & binsy == q ;
          cnt = cnt + 1 ;
      end
  end
end
