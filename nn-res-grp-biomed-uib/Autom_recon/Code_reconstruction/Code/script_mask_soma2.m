
% load the relevant data
skelim = cbiReadNifti_unix(nii_skel_fname);
load('I3.mat')

%% ask user to draw soma boundaries if not done before
if exist('somamask_manual.mat','file')
    load('somamask_manual.mat');
    disp('Mask of soma loaded.')
else
    txt = 'please click boundaries of area containing the soma; press enter to finish, backspace to correct';
    somamaskXY = cutout_wskel(I3,txt,skelim);
    somamaskXZ = cutout_wskel(permute(I3,[1 3 2]),txt,permute(skelim,[1 3 2]));
    save somamask_manual.mat somamaskXY somamaskXZ
    disp('Mask for soma created.')
end

somamaskXY3d = repmat(somamaskXY,[1,1,size(I3,3)]);
somamaskXZ3d = repmat(permute(somamaskXZ,[1 3 2]),[1,size(I3,2),1]);
somamask3d = somamaskXZ3d & somamaskXY3d;
save somamask3d.mat somamask3d

%% plot results
% duplicate the skeleton and cancel out all the points inside the soma
c2 = skelim;
c2(somamask3d==1) = 0;

% % add a skeleton point in x1,x2,x3 which should be around the center of the
% % soma
% load somapkt.mat
% load zsnittmax.mat
% x1 = round(somapkt(1));
% x2 = round(somapkt(2));
% x3 = round(zsnittmax);
% c2(x1,x2,x3) = 255;

% compare the skeleton before and after
figure('Name','Skeletonization result')
subplot(1,3,1);imagesc(max(skelim,[],3));colormap(gray);axis equal; axis tight
title('Skeletonization by FIJI')
subplot(1,3,2);imagesc(max(c2,[],3));colormap(gray);axis equal; axis tight
title('Voxels inside user-denoted region (soma) removed')
subplot(1,3,3);imagesc(max(somamask3d,[],3));colormap(gray);axis equal; axis tight
title('User denoted region (soma)')
%save('soma_mask.mat','soma_mask_center','soma_mask_radius')
