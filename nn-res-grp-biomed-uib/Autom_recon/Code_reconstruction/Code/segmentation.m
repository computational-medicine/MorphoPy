function segmentation(im, voxel_size, add2threshold, adaptive_filtersize, cohenc, Nvoxels_segment_min, simpleThreshold)
%% mask background (and pipette)
% and scale intensities from 0 - 1000
seg4 = mask_background(im);
clear im
disp('Neuron masked.')

%% image smoothing
im_cohendiff = smooth_image(seg4,cohenc,voxel_size);
clear seg4

%% thresholding (adaptive and simple)
BW2 = thresholding(im_cohendiff,voxel_size, add2threshold, adaptive_filtersize, Nvoxels_segment_min, simpleThreshold);

%% Crop stacks
[I2,t1,t2,t3,t4,t5,t6] = cropvolum_nevron(BW2);
save t1.mat t1; save t2.mat t2; save t3.mat t3; save t4.mat t4; save t5.mat t5; save t6.mat t6
save I2.mat I2   % cropped initially segmented image
s_crop = size(I2);
save s_crop.mat s_crop
I3 = im_cohendiff(t3:t4,t1:t2,t5:t6); % cropped, smoothed image
save I3.mat I3

end

function BW2 = thresholding(im,voxel_size, add2threshold, adaptive_filtersize, Nvoxels_segment_min, simpleThreshold)
thresh = simpleThreshold*1000;
seg4=single(zeros(size(im)));
seg4(im>=thresh)=1;

% keep biggest object
[bw] = bwkeep(seg4,1,26);   % 26 denotes 3d voxel connectivity (all neighbours)
figure('Name','Simple thresholding for soma and primary branch')
imagesc(max(im,[],3)); hold on
%colormap gray
contour(max(bw,[],3),1,'r')
save bw.mat bw;
disp('Simple threshold applied.')

if exist('th.mat','file') && exist('imth.mat','file')
    load th.mat th
    load imth.mat imth
    disp('Adaptive treshold without addition loaded')
    bw2 = (imth>(th+add2threshold));
    
else
    im4=im./max(im(:));
    save im4.mat im4
    display(add2threshold)
    disp('Adaptive filtering/thresholding...')
    tic
    [bw2,th,imth] = adaptfiltim(im4,adaptive_filtersize,add2threshold,voxel_size*10e2);
    th = th - add2threshold;
    save th.mat th;        % obtained threshold without addition
    save imth.mat imth;      % scaled image for thresholding
    save bw2x.mat bw2;       % black/white thresholded image
    disp('Image adaptive filtered/thresholded.')
    toc 
end    
c3b=bw+bw2;
c3b(c3b>=1)=1;

figure('Name','Adaptive thresholding');
imagesc(max(im,[],3)); hold on
%colormap gray
contour(max(bw2,[],3),1,'r')

% Remove small objects (probably noise)
BW2 = bwareaopen(c3b, Nvoxels_segment_min);
save BW2.mat BW2

figure('Name','Combined, with small objects removed');
imagesc(max(im,[],3)); hold on
%colormap gray
contour(max(BW2,[],3),1,'r')

end


function seg4 = mask_background(im)
% ask user to draw boundaries of background if not done before
if exist('neuronmask.mat','file')
    load('neuronmask.mat');
    disp('Mask of neuron loaded.')
else
    neuronmask = cutout(max(im,[],3),'please click boundaries of area containing the neuron; press enter to finish, backspace to correct');
    save neuronmask.mat neuronmask;
    disp('Mask for neuron created.')
end

% scale original between 0 and 1000
im = double(im);
maxorig=max(im(:));
minorig=min(im(:));

orig=((im-minorig)./(maxorig-minorig))*1000;
save orig.mat orig
seg4 =double(orig).*repmat(double(neuronmask),1,1,size(im,3));
save seg4.mat seg4

end

function im_cohendiff = smooth_image(im,cohenc,voxel_size)

if exist('im_cohendiff.mat','file')
    load im_cohendiff.mat
    disp('Smoothed image stack loaded.')
else
    tic
    disp('Smoothing image stack...')
    im_cohendiff = zeros(size(im));
    for i=1:size(im,3)
        im_cohendiff(:,:,i) = cohenhdiff(im(:,:,i),cohenc.DT,cohenc.NITER,cohenc.KAPPA,voxel_size*10e2);
    end
    save im_cohendiff.mat im_cohendiff;
    disp('Images in stack smoothed.')
    toc
end

end
