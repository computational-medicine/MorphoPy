function binary_to_swc(downsamp_skel_pnts,pos1,pos2,pos3,cleanfactor,swcfilename, swcfilename_smooth, Nminpts_trees,MST,nii_seg_fname,nii_skel_fname)
%% load skeleton and segmentation
skel_shortest_path_com = cbiReadNifti_unix(nii_skel_fname);
skel_shortest_path_com = permute(skel_shortest_path_com,[2 1 3]);
BWcell = cbiReadNifti_unix(nii_seg_fname);

%% if desired, downsample points
load('somapkt.mat')
load('zsnittmax.mat')
som2(1) = somapkt(1);
som2(2) = somapkt(2);
som2(3) = zsnittmax;
dum=0;
k=find(skel_shortest_path_com==255);

k2=k(1:downsamp_skel_pnts:end);
skel_shortest_path_com2=zeros(size(skel_shortest_path_com));
skel_shortest_path_com2(k2)=255;
while 1
    nederste_z=max([round(som2(3))-dum;1]);
    test10=skel_shortest_path_com2(round(som2(2))-dum:round(som2(2))+dum,round(som2(1))-dum:round(som2(1))+dum,nederste_z:round(som2(3))+dum);
    test = find(test10==255);
    if length(test) >= 1;
        test2 = test(1);
        [xtest,ytest,ztest] = ind2sub(size(test10),test2);
        koord = [round(som2(2))-dum+xtest-1 round(som2(1))-dum+ytest-1 nederste_z+ztest-1];
        break;
    end
    dum=dum+1;
end

somind = sub2ind(size(skel_shortest_path_com),koord(1),koord(2),koord(3));

%% Find process diameters by calculating the distance transform
BWcell2=zeros(size(BWcell));
for i=1:size(BWcell,3)
    BWcell2(:,:,i)=imcomplement(BWcell(:,:,i));
end
load('voxel_size.mat');
somadist = bwdistsc(BWcell2,[voxel_size(1) voxel_size(2) voxel_size(3)]);
k3=k2;

%% get real world coordinates
load('t1.mat')
load('t2.mat')
load('t3.mat')
load('t4.mat')
load('t5.mat')
load('t6.mat')
load('seg4.mat');
fullt_volum=zeros(size(seg4));
clear seg4
skel_shortest_path_com3=zeros(size(skel_shortest_path_com));
skel_shortest_path_com3(k3)=1;

fullt_volum(t1:t2,t3:t4,t5:t6)=skel_shortest_path_com3;
k4=find(fullt_volum==1);

load('voxel_size.mat');
[s1,s2,s3] = ind2sub(size(fullt_volum),k4);
somind2 = find(k3==somind);
s1b=(s1+pos1).*voxel_size(1)*10e2;
s2b=-(s2+pos2).*voxel_size(2)*10e2;
s3b=-(s3+pos3).*voxel_size(3)*10e2;

%% solve Minimum Spanning Tree problem
% ST3 = MST_tree (somind2, s1b, s2b, s3b, 0.4, 10000, [], []);
ST3 = MST_tree (somind2, s1b, s2b, s3b, MST.bf, MST.maxdist, [], []);
figure, plot_tree(ST3,[],[],[],[],'-2l');

somadistmm=somadist.*10e2;
somadistmm2=zeros(size(somadistmm,2),size(somadistmm,1),size(somadistmm,3));
for i=1:size(somadistmm,3)
    somadistmm2(:,:,i) = flipud(rot90(somadistmm(:,:,i)));   
end

somadistmm3 = zeros(size(fullt_volum));
somadistmm3(t1:t2,t3:t4,t5:t6) = somadistmm2;

for i=1:length(ST3.X)   % used to be length(s1), which for some reason 
                        % deviates one from the nr of points in ST3 on
                        % rare occasions

    % Note that the tree toolbox stores the diameters in tree.D, while the SWC
    % file format oddly stores the radius, but denotes it with 'D'.
    xkoord = (ST3.X(i)/(voxel_size(1)*10e2))-pos1;
    ykoord = -(ST3.Y(i)/(voxel_size(2)*10e2))-pos2;
    zkoord = -(ST3.Z(i)/(voxel_size(3)*10e2))-pos3;

    ST3.D(i) = somadistmm3(round(xkoord),round(ykoord),round(zkoord))*2;
    
    % Change the flags (segment type, R) to 3=dendrite
    ST3.R(i) = 3;

end
%figure, plot_tree(ST3);
%figure, plot_tree(ST3,[],[],[],[],'-3l');
save ST3.mat ST3

% BJZ (151214): clean short branches from the tree
tree3 = clean_tree(ST3, cleanfactor, '-s');
disp('Tree is cleaned')

%% BJZ Delete the branchpoints inside the soma, except for the first point
load somamask3d.mat
somamask_large = zeros(size(fullt_volum));
somamask_large(t1:t2,t3:t4,t5:t6)=permute(somamask3d,[2,1,3]);

todelete = [];
for i=2:length(tree3.X)
    x = round((tree3.X(i)/(voxel_size(1)*10e2))-pos1);
    y = round(-(tree3.Y(i)/(voxel_size(2)*10e2))-pos2);
    z = round(-(tree3.Z(i)/(voxel_size(3)*10e2))-pos3);
    if somamask_large(x,y,z)==1;
        todelete = [todelete,i]; %#ok<AGROW>
    end
end
% remove trees with less then Nminpts_branch points
tree3 = delete_tree(tree3,todelete);     % remove the points in the soma
rootchilds = find(idpar_tree(tree3)==1); % find children of the main root, i.e. the actual start of all trees
Ndown = child_tree(tree3);               % downstream points
smalltreeroots = rootchilds(Ndown(rootchilds)<Nminpts_trees);         % find origins of trees with less then 4 points
todelete = [];                           % collect all nodes in the small subtrees
for i = smalltreeroots';
    todelete = [todelete,(find(sub_tree(tree3,i)==1))']; %#ok<AGROW>
end
tree3 = delete_tree(tree3,todelete);     % remove the subtrees
tree3 = setlabels_todendrite(tree3); % fix the labels

save tree3.mat tree3;
swc_tree(tree3, swcfilename);
disp('Branch points inside soma removed')
%% smooth tree
tree_smooth = smooth_tree(tree3);
tree_smooth = setlabels_todendrite(tree_smooth); %  Necessary to reset?

save tree_smooth.mat tree_smooth;
swc_tree(tree_smooth, swcfilename_smooth);
disp('Smoothed tree also generated')
end

function tree = setlabels_todendrite(tree)
tree.R(1) = 1;       % set starting point to be the soma
tree.D(1) = 0;       % with radius 0
for i=2:length(tree.R)
    tree.R(i) = 3;   % set description tag to 3: dendrites
end

end