%data configs
%VOCinit
%traintest = 'trainval';
%ids = textread(sprintf(VOCopts.imgsetpath, traintest), '%s');
imgset = 'D:/Data/PASCAL/VOC2007/JPEGImages';
imfls = dir(fullfile(imgset, '*.jpg'));
savDir = 'D:/Data/PASCAL/VOC2007/matlab/siftPhow/ss2';

%feature configs
sstep = 2;

tic
for i = 1:length(imfls)

    fprintf('image %d\n', i) ;

    imname = fullfile(imgset, imfls(i).name);
    im = imread(imname);
    im = rgb2gray(im);
    
    [frs, descrs] = vl_phow(single(im), 'Step', sstep);

    image.path = imname;
    image.frames = frs;
    image.height = size(im, 1);
    image.width = size(im, 2);
    image.step = sstep;
    image.scale = 'default-phow';
    image.sift = descrs;

    [filpat, filnam, filext] = fileparts(imfls(i).name);
    save('-v7.3', fullfile(savDir, [filnam '_allVariables.mat']), 'image'); 
end
toc


