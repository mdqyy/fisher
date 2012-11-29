%% Start using ERS

VOCinit
VOCprep

load /home/gavves/u015425/projects/matlab/mine/fisherNet/voc2007_bagofwords_model-nwords4000.mat
load /home/gavves/u015425/projects/matlab/mine/fisherNet/voc2007svmmodelUnnormalized-onlyBackgroundNegatives_fordetection_bag_of_words.mat

partition = 'test' ;
superpixelmethod = 'turbopixels' ;
arbelaezpath = '/home/gavves/u015425/projects/matlab/mine/fisherNet/gPb' ;
numSuperpixels = 200 ;
malikthres = 0.1 ;
saveto = '/home/gavves/u015425/projects/matlab/mine/fisherNet/voc-ers-bow-onlybgnegatives-turbopixels' ;

%% Scan category
c = 1 ;
class = VOCopts.classes{c} ;
[ids, gt] = textread(sprintf(VOCopts.clsimgsetpath, class, partition), '%s %d');
ids = ids(gt == 1) ;
ids = str2double(ids) ;
files = imfls(ids) ;

%% Get an image
for q = 1 : numel(files)
    eta(q, numel(files)) ;
    imc = imread(filepath(files(q))) ;
    img = imreadgray(filepath(files(q))) ;

    %% Get image superpixels using Turbopixels
    spix.superpixelmethod = superpixelmethod ;
    spix.numSuperpixels = numSuperpixels ;
    spix.malikthres = malikthres ;
    spix.arbelaezpath = arbelaezpath ;
    spix.impath = filepath(files(q)) ;
    if strcmp(spix.superpixelmethod, 'turbopixels')
        [phi, boundary, spmap] = superpixels(im2double(img), spix.numSuperpixels);
        spix.spmap = spmap ;
        spix.phi = phi ;
        spix.boundary = boundary ;
        spix.numSuperpixels = numel(unique(spix.spmap)) ;
    elseif strcmp(spix.superpixelmethod, 'arbelaez')
        [spmap, boundary, ucm2] = arbelaezMalikRegions(spix.impath, spix.arbelaezpath, spix.malikthres) ;
        spix.spmap = spmap ;
        spix.boundary = boundary ;
        spix.ucm2 = ucm2 ;
        spix.numSuperpixels = numel(unique(spix.spmap)) ;
    end
    
    %% Create (node, edge) graph
    graph = superpixelsToGraph(spix.boundary, spix.spmap) ;
    graph = (graph - 1)' ; % In the ERS code ONLY the graph nodes start from 0, not the superpixel labels.
    graph = [spix.numSuperpixels * ones(size(graph, 1), 1), graph] ;

    %% Construct Adjacency Matrix
    adjmat = zeros(graph(1),graph(1));
    idx = sub2ind(size(adjmat), graph(:,2)+1, graph(:,3)+1);
    adjmat(idx) = 1;

    %% Calculate Features Histogram: Bag-of-Words, etc. NOTE: The histogram should be unnormalized
    [frames, sift] = vl_phow(single(img), 'Step', model.ss, 'Sizes', model.scale) ;
    nzix = sum(sift, 1) ~= 0 ;     % OPTIONAL
    sift = single(sift(:, nzix)) ; % OPTIONAL, although if there are too many 0 vectors, performance decreases significantly
    frames = frames(:, nzix) ;
    [dummy, binsa] = bagofwords(model, sift) ;
    wordmap = zeros(size(img)) ;
    wordmap(sub2ind(size(wordmap), frames(2, :), frames(1, :))) = binsa ;
    hist = zeros(model.numWords, spix.numSuperpixels) ;
    for i = 1 : spix.numSuperpixels
        ix = find(spix.spmap == i) ;
        binsa = wordmap(ix) ;
        hist(:, i) = vl_binsum(hist(:, i), ones(size(binsa)), binsa) ;
    end

    %% Superpixel SVM score
    spix_scores = hist' * svmmodel{c}.w(:, 1) ;

    %% Solve MWCS problem
    [nodes, scores] = solveMWCS(spix_scores, adjmat) ;
    
%     try
%         [nodes, scores] = solveMWCS(spix_scores, adjmat, graph) ;
%     catch ME
%         continue ;
%     end

    %%% Visualizations
    %%%% Compute the segmentation mask result
    bcmask = zeros(size(spix.spmap));
    for n=1:length(nodes)
      idx = find(spix.spmap == nodes(n));
      bcmask(idx) = 1;
    end

    %%%% Get segmented image
    res_segim = zeros(size(imc));
    for n=1:3
      tmpim = imc(:,:,n);
      res_segim(:,:,n) = double(tmpim) .* bcmask;
    end

    save('-v7.3', [saveto '/voc_ers-c' num2str(c) '-' num2str(q) '.mat'], 'bcmask', 'res_segim', 'spix', 'nodes', 'imc') ;

end



%%

class = 1 ;
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w'} ;

% cmap = colormap(jet(30)) ;

% while true
% i = randpick(10, 1) ;
for i = 1 : 10
load(['voc-ers-bow-onlybgnegatives-turbopixels/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
% load(['voc-ers-bow/voc_ers-' num2str(i) '.mat']) ;

figure(1)
%%%% Plot results
subplot(2,2,1), imshow(imc);
freezeColors
subplot(2,2,2), imshow(spix.spmap, []), colormap('jet')
freezeColors
B = bwboundaries(bcmask) ;
% subplot(2,2,3), imshow(bcmask);
subplot(2,2,3), imshow(uint8(res_segim));
freezeColors
subplot(2,2,4)
imshow(imc);
hold on
if ~isempty(B)
    for j = 1 : size(B)
        plot(B{j}(:, 2), B{j}(:, 1), 'color', colors{j}, 'LineWidth', 5)
    end
end
% subplot(2,2,4), imshow(uint8(res_segim));
freezeColors
pause
clf
end

%%
class = 1 ;
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w'} ;

% cmap = colormap(jet(30)) ;

% while true
% i = randpick(10, 1) ;
for i = 21 : 30
% ldb = load(['voc-ers-fishernet-onlybgnegatives/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
ldb = load(['voc-ers-bow-onlybgnegatives/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
ldf = load(['voc-ers-bow-onlybgnegatives-turbopixels/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;

%%%% Plot results
figure(1)
clf
subplot(2,2,1), imshow(ldb.imc);
freezeColors
subplot(2,2,2), imshow(ldb.spix.spmap, []), colormap('jet')
freezeColors
B = bwboundaries(ldb.bcmask) ;
% subplot(2,2,3), imshow(bcmask);
subplot(2,2,3), imshow(uint8(ldb.res_segim));
freezeColors
subplot(2,2,4)
imshow(ldb.imc);
hold on
if ~isempty(B)
    for j = 1 : size(B)
        plot(B{j}(:, 2), B{j}(:, 1), 'color', colors{j}, 'LineWidth', 5)
    end
end
% subplot(2,2,4), imshow(uint8(res_segim));
freezeColors

figure(2)
clf
subplot(2,2,1), imshow(ldf.imc);
freezeColors
subplot(2,2,2), imshow(ldf.spix.spmap, []), colormap('jet')
freezeColors
B = bwboundaries(ldf.bcmask) ;
% subplot(2,2,3), imshow(bcmask);
subplot(2,2,3), imshow(uint8(ldf.res_segim));
freezeColors
subplot(2,2,4)
imshow(ldf.imc);
hold on
[br, bc] = find(ldf.spix.boundary) ;
plot(bc, br, 'y.')
if ~isempty(B)
    for j = 1 : size(B)
        plot(B{j}(:, 2), B{j}(:, 1), 'color', colors{j}, 'LineWidth', 5)
    end
end
% subplot(2,2,4), imshow(uint8(res_segim));
freezeColors
pause

end













