function fastmarching(NbrightPixelsSoma, fname)

load I2.mat;
load I3.mat;

%% Label segments
BW3=bwlabeln(I2);
save BW3.mat BW3

%% find soma centre
[somapkt,zsnittmax]=gensomapkt(I2,NbrightPixelsSoma);
save somapkt.mat somapkt
save zsnittmax.mat zsnittmax

%% Fast Marching from the soma
disp('Find paths to soma using fast marching...')
Forgr_c = pathToSoma(I3,BW3,somapkt,zsnittmax);
disp('Paths to soma generated.')

%% Join voxels of paths with segmentation
Join_and_save(BW3, Forgr_c, fname)
%lage_binaer_celle_2_del1_test_tynn_forgreining2(1);
disp('Segmentation done!')

% plot results
figure('Name', 'Fast Marching results')
imagesc(max(BW3,[],3))
hold on;
contour(max(Forgr_c,[],3),1,'r')
plot(somapkt(2),somapkt(1),'gx')
end

function [somapkt,zsnittmax]=gensomapkt(I2,NbrightPixelsSoma)
% Approximate soma center as centre of mass of NbrightVoxelsSoma brightest pixels
% in an additive projection of the distance transform of the segmented
% image. The distance transform highlights wide, thick structures.

dist_image=bwdist(imcomplement(I2));

[s1,s2,~]=size(dist_image);
addproj = sum(double(dist_image),3);

% Approximate soma center xy-coordinates as centre of mass of NbrightVoxelsSoma brightest pixels
[~,addproj_sort_ind]=sort(reshape(addproj,1,s1*s2),'descend');
[indtysx,indtysy]=ind2sub(size(addproj),addproj_sort_ind(1:NbrightPixelsSoma));
cmassi=mean(indtysx(:));
cmassj=mean(indtysy(:));
somapkt=[cmassi,cmassj];

% find z-coordinate with maximum distance from background at xy-coordinates
zsnitt=squeeze(dist_image(round(somapkt(1)),round(somapkt(2)),:));
[~,zsnittmax]=max(zsnitt);

end

function Join_and_save(BW3,Forgr_c,fname)
BW3b = zeros(size(BW3));
BW3b(BW3>=1) = 1;
BWwholecell = BW3b + Forgr_c;
BWwholecell(BWwholecell>=1) = 1;

% save as nifti (.nii) file
load('standardheader.mat'); % loads header
header.dim(2) = size(BWwholecell,1);
header.dim(3) = size(BWwholecell,2);
header.dim(4) = size(BWwholecell,3);
if isunix==1
    cbiWriteNifti_unix(fname,BWwholecell,header);
else
    cbiWriteNifti(fname,BWwholecell,header);
end
load('I3.mat');
if isunix==1
    cbiWriteNifti_unix('I3.nii',I3,header);
else
    cbiWriteNifti('I3.nii',I3,header);
end
end