%% Start using ERS with FisherNet

clear

VOCinit
VOCprep

load /home/gavves/u015425/projects/matlab/mine/fisherNet/voc2007model_fordetection.mat
ldU = load('/home/gavves/u015425/projects/matlab/mine/fisherNet/voc2007svmmodelUnnormalized-onlyBackgroundNegatives_fordetection_fisherNet.mat') ;

partition = 'test' ;
superpixelmethod = 'arbelaez' ;
arbelaezpath = '/home/gavves/u015425/projects/matlab/mine/fisherNet/gPb' ;
numSuperpixels = 200 ;
malikthres = 0.1 ;
saveto = '/home/gavves/u015425/projects/matlab/mine/fisherNet/voc-ers-fishernet-onlybgnegatives-ss1' ;
model.ss = 1 ;

%% Scan category
c = 3 ;
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
    
    %% Get FisherNet on superpixels
    spix.superpixelmethod = superpixelmethod ;
    spix.numSuperpixels = numSuperpixels ;
    spix.malikthres = malikthres ;
    spix.arbelaezpath = arbelaezpath ;
    spix.impath = filepath(files(q)) ;
    [fnet, fnorm, spix] = fisher_net_superpixels( imc, spix, model, ldU.svmmodel{c}.w(:, 1) ) ;

    %% Create (node, edge) graph
    graph = superpixelsToGraph(spix.boundary, spix.spmap) ;
    graph = (graph - 1)' ; % In the ERS code ONLY the graph nodes start from 0, not the superpixel labels.
    graph = [spix.numSuperpixels * ones(size(graph, 1), 1), graph] ;

    %% Construct Adjacency Matrix
    adjmat = zeros(graph(1),graph(1));
    idx = sub2ind(size(adjmat), graph(:,2)+1, graph(:,3)+1);
    adjmat(idx) = 1;

    %% Superpixel SVM score
    spix_scores = vec(fnet) ;

    %% Solve MWCS problem
    [nodes, scores] = solveMWCS(spix_scores, adjmat) ;

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
class = 4 ;
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w'} ;

% cmap = colormap(jet(30)) ;

figure
while true
i = randpick(10, 1) ;
load(['voc-ers-bow-onlybgnegatives/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
% load(['voc-ers-bow/voc_ers-' num2str(i) '.mat']) ;

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
class = 6 ;
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w', 'r', 'g', 'b', 'm', 'c', 'y', 'w'} ;

% cmap = colormap(jet(30)) ;

% while true
% i = randpick(10, 1) ;
for i = 1 : 100
% ldb = load(['voc-ers-fishernet-onlybgnegatives/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
ldb = load(['voc-ers-bow-onlybgnegatives/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;
ldf = load(['voc-ers-fishernet-onlybgnegatives-ss1/voc_ers-c' num2str(class) '-' num2str(i) '.mat']) ;

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
if ~isempty(B)
    for j = 1 : size(B)
        plot(B{j}(:, 2), B{j}(:, 1), 'color', colors{j}, 'LineWidth', 5)
    end
end
% subplot(2,2,4), imshow(uint8(res_segim));
freezeColors
pause

end

















