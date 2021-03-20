% Code to run segmentation and skeletonization of the image stacks, i.e. determine which pixels belong to the
% cell volume, make a skeleton and generate a tree.
% Usage:
% - make a folder with as name the cell name (e.g. m140128_1)
%       - must contain the image stack converted to nrrd format (use FIJI).
%       - must contain a file named voxel_size.mat that contains an array
%       voxel_size = [DX, DY, DZ] corresponding to image stack voxel size in mm
% - make an empty folder named results_m140128_1
% - run the script
%
% User input:
%   - select the NRRD image stack
%   - manually denote the general area of the cell in 2D projection
%   - manually denote the soma region

%% 1. Clear previous session
clear variables
path(pathdef)       % revert to the default path dependencies

%% 2. Parameters

% % If the reconstruction is rerun, redo any previous steps? 
redo_segmentation = true;
redo_fastmarching = true;
redo_skeletonization = true;
redo_treegeneration = true;
plot_results = true;

% % Segmentation
Nvoxels_segment_min = 20; % Minimum nr of voxels in a segment (smaller segments are discarded)
add2threshold = 0.05;     % Value added to adaptive threshold
adaptivefiltersize = [2.5 2.5 5]; % [um]
%adaptivefiltersize = [10 10 20]; % [xy_voxel_sizes], used for the Diadem Olfactory Projection fibers
simpleThreshold = 0.11; % Threshold for soma and apical denrite (fraction of maximum image intensity)
NbrightPixelsSoma = 100;

% Coherence enhancing diffusion filtering parameters (smoothing)
cohenh.DT = 0.2;
cohenh.NITER = 15;
cohenh.KAPPA = 0.001;

% % Tree generation
Nsample = 1;                % every Nth point in skeleton is sampled (downsampling)
clean_tree_scaling = 0.2;   % larger values denote more aggressive cleaning
Nminpts_trees = 4;          % remove trees (branches extending from the soma) with less then N points
MST.bf = 0.1;              % balancing factor for MST (should be small if skeleton is accurate)
MST.maxdist = 10;           % [um], maximum distance between connected points in MST

% % Soma contour generation
Dth = 36;          %degrees, a larger angle cuts off more branches and makes the soma rounder
maxSomaSizeZ = 20; %[um], length of stack in z-direction to analyse: Set to double the typical soma size.
SmoothLength = 8;  %[um], should be smaller than soma diameter: Set to ~half of typical soma diameter

% % FIJI path
FIJIpath = '/Applications/Fiji.app/plugins'; % FIJIpath = 'D:\fiji-win64\Fiji.app\plugins';
FIJIscriptpath = '/Applications/Fiji.app/scripts'; % FIJIscriptpath = 'D:\fiji-win64\Fiji.app\scripts';

%% 3. Set paths

% Add directories to Matlab path
[nrrdImageStack,sourcepath] = uigetfile('*.nrrd','Select NRRD image stack');
[~, sourceFolderName] = fileparts(sourcepath(1:end-1));

resultFolderName = ['results_',sourceFolderName];
cd(sourcepath); cd('..'); mkdir(resultFolderName); cd(resultFolderName)

addpath(genpath('../../Code_reconstruction/Code'))
addpath(genpath('../../Code_reconstruction/Miji'))
addpath(genpath('../../Code_reconstruction/xtra'))
% contains imagestack of the cell (in nrrd format), manual segmentation (swc format), and the voxel sizes
addpath(genpath(['../',sourceFolderName]))
addpath(genpath(['../',resultFolderName]))
addpath('..')
% FIJI software paths
addpath(FIJIpath)
addpath(FIJIscriptpath)

%% 4. Segmentation (thresholding)
im=nrrdread(nrrdImageStack);
load voxel_size;

if (redo_segmentation || ~exist('I2.mat','file') )
    tic
    segmentation(im, voxel_size, add2threshold, adaptivefiltersize,cohenh, Nvoxels_segment_min, simpleThreshold)
    toc
else
    disp('Thresholding results found.')
end

%% 5. Segmentation (fast marching)

nii_seg_fname='segmented.nii';

if (redo_fastmarching || ~exist(nii_seg_fname,'file') )
    tic
    fastmarching(NbrightPixelsSoma, nii_seg_fname)
    toc
else
    disp('Segmentation results found.')
end

%% 6. Skeletonization

nii_seg_fname='segmented.nii';
segmented_stack=cbiReadNifti_unix(nii_seg_fname);
header=cbiReadNiftiHeader_unix(nii_seg_fname);
nii_skel_fname='skeleton.nii';

if (redo_skeletonization || ~exist(nii_skel_fname,'file'))
    % This runs Fiji from within matlab to compute the skeleton of the segmented cell
    skeleton=Nevronfijimatlab3(segmented_stack);
    cd(resultFolderName); % path is changed, so change back
    cbiWriteNifti_unix(nii_skel_fname,skeleton,header);
else
    disp('Skeletonization results found.')
end

%% 7. Tree generation and conversion to SWC-file
swcfilename = [sourceFolderName,'_auto.swc'];
swcfilename_smooth = [sourceFolderName,'_auto_smooth.swc'];

if (redo_treegeneration || ~exist(swcfilename,'file')) 
    script_mask_soma2
    tic
    binary_to_swc(Nsample,0,0,0,clean_tree_scaling,swcfilename,swcfilename_smooth,Nminpts_trees,MST,nii_seg_fname,nii_skel_fname);
    toc
    disp('Tree generated and SWC files written.')
    tree = load_tree(swcfilename);
    figure('Name','Tree after cleaning and spurious branch removal');
    plot_tree(tree,'r',[],[],[],'-3l');
else
    disp('SWC files found.')
    tree = load_tree(swcfilename);
end

%% 8. Soma contour generation
somafilename = [sourceFolderName,'_auto_soma.swc'];
somaContour = maxSomaContour(Dth, maxSomaSizeZ, SmoothLength,nii_seg_fname);
make_tree_with_soma(somaContour,tree,somafilename);

% note that converting this file to Neurolucida .ASC format using NLMorphologyConverter results in a
% proper file, containing 1 soma contour and the various trees, which is
% imported correctly into NEURON with the import tool.
% the SWC file itself should not be directly imported, since the soma contour
% is likely to be split across the file, which NEURON does not correctly import.



%% 9. Plotting
if plot_results
    figurescale = [0,1];
    im=nrrdread(nrrdImageStack);
    im = double(im);
    im = im/max(im(:));
    s = size(im);
    load voxel_size
    xim = ([1,s(1)]*voxel_size(1))*1e3;
    yim = ([-1,-s(2)]*voxel_size(2))*1e3;
    zim = ([-1,-s(3)]*voxel_size(3))*1e3;
    figure('Name','3D Tree above MIP');
    imagesc(([0,s(1)-1]*voxel_size(1))*1e3,([-0,-s(2)+1]*voxel_size(2))*1e3,max(im,[],3),figurescale)
    axis xy; hold on
    tree_shifted = tree; tree_shifted.Z = tree_shifted.Z + 100;
    h = plot_tree(tree_shifted,[1,0,0]);
    h.FaceAlpha = 0.4;
    plot3(somaContour.X,somaContour.Y,somaContour.Z+100,'g')
   
end

