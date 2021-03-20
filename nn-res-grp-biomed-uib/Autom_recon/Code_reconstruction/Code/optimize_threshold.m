% Script to find optimal thresholding values

%% 1. Initiation
clear variables

load im_cohendiff.mat % smoothed image stack

load imth.mat % image stack that is adaptively thresholded
load th.mat   % the adaptive threshold




%% 2. Calculate and plot result of thresholding

Nvoxels_segment_min = 20; % Minimum nr of voxels in a segment (smaller segments are discarded)
add2threshold = 0.2;     % Value added to adaptive threshold
simpleThreshold = 0.2; % Threshold for soma and apical denrite (fraction of maximum image intensity)

% simple thresholding
terskel = simpleThreshold*1000;
seg4=single(zeros(size(im_cohendiff)));
seg4(im_cohendiff>=terskel)=1;

% keep biggest object
bw = bwkeep(seg4,1,26);   % 26 denotes 3d voxel connectivity (all neighbours)
figure(343)
subplot(2,2,1)
hold off;
imagesc(max(im_cohendiff,[],3)); hold on
contour(max(bw,[],3),1,'r')
title('Simple thresholding for soma and primary branch')
 
im4=im_cohendiff./max(im_cohendiff(:));

bw2 = imth>(th+add2threshold);

c3b=bw+bw2;
c3b(c3b>=1)=1;

subplot(2,2,2)
imagesc(max(im_cohendiff,[],3)); hold on
%colormap gray
contour(max(bw2,[],3),1,'r')
title('Adaptive thresholding');

% Remove small objects (probably noise)
BW2 = bwareaopen(c3b, Nvoxels_segment_min);
subplot(2,2,3)
imagesc(max(im_cohendiff,[],3)); hold on
%colormap gray
contour(max(BW2,[],3),1,'r')
title('Combined, with small objects removed');