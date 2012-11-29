VOCinit
imgset = 'trainval';
ids = textread(sprintf(VOCopts.imgsetpath, imgset), '%s');
%savDir = 'D:/Data/PASCAL/VOC2007/matlab/siftPhow/ss3';

N = 500000;
nFeatsPerIm = ceil(N / numel(ids)) + 1;
sift = zeros(128, nFeatsPerIm * numel(ids));

%feature configs
model.ss = 2;

cnt = 1;
for i = 1 : numel(ids)
    %eta(i, numel(ids));
    %ld = load([savDir '/' ids{i} '_allVariables.mat']);
    
    fprintf('image %d\n', i);
    imname = fullfile(VOCopts.datadir, VOCopts.dataset, 'JPEGImages', [ids{i}, '.jpg']);
    
    im = imread(imname);
    im = rgb2gray(im);
    
    [frs, descrs] = vl_phow(single(im), 'Step', model.ss);
    nzix = sum(descrs, 1) ~= 0;  % OPTIONAL
    descrs = descrs(:, nzix);    % OPTIONAL, although if there are too many 0 vectors, performance decreases significantly
    %frs = frs(:, nzix);
    
    sift(:, cnt : cnt + nFeatsPerIm - 1) = vl_colsubset(descrs, nFeatsPerIm);
    cnt = cnt + nFeatsPerIm;
end

save -v7.3 D:/Data/PASCAL/VOC2007/matlab/siftPhow/voc2007_sampledSiftVectors500K_phow_ss2.mat sift