function Forgr_c = pathToSoma(I3,BW3,somapkt,zsnittmax)

SpeedImage=double(I3);
SpeedImage=SpeedImage/max(SpeedImage(:));
SpeedImage(SpeedImage<=1e-4)=1e-8; % originally set to 1e-4


% find source points, i.e. the soma
% somapkt for x og y-coordinates
% zsnittmax for z koordinate.
SourcePoint=round([somapkt(1) somapkt(2) zsnittmax]);

tic
% perform fast marching
[T,~] = msfm(SpeedImage, SourcePoint', true, true);
save T.mat T
disp('Marching time from soma calculated.')
toc

%% backtracing

% don't trace from the soma area
sar = BW3(round(SourcePoint(1)),round(SourcePoint(2)),round(SourcePoint(3)));

% find shortest path from all sections to soma
antalcomp=max(BW3(:));
for j=1:antalcomp
        tic
        fprintf('Connecting component %d out of %d ...\n', j, antalcomp)
        % find the voxel in segment j that has the shortest FM arrival time
        % from the soma, and start backtracing from there
        ks= find(BW3==j);
        [minarrtime,minarrtime_ind]=min(T(ks));
        if j~=sar
            [minarrtime_x,minarrtime_y,minarrtime_z]=ind2sub(size(BW3),ks(minarrtime_ind));
            StartPoint = [minarrtime_x minarrtime_y minarrtime_z];
            ShortestPath{j} = shortestpath(T,StartPoint',SourcePoint',0.5,'rk4'); %#ok<AGROW>
        else
            ShortestPath{j} = []; %#ok<AGROW>
        end
        toc
end

Forgr_c = zeros(size(BW3));
Forgr_c_thick = zeros(size(BW3));
s_crop = size(BW3);
for j=1:antalcomp
    path = ShortestPath{j};
    if ~isempty(path)
        [stor1,stor2]=size(path);
        if stor1==3 && stor2~=3
            path=path';
        end
        if min(round(path(:,1)))>=2 && min(round(path(:,2)))>=2 && min(round(path(:,3)))>=2 && max(round(path(:,1)))<=s_crop(1)-1 && max(round(path(:,2)))<=s_crop(2)-1 && max(round(path(:,3)))<=s_crop(3)-1
            for k=1:size(path,1)
                Forgr_c(round(path(k,1)),round(path(k,2)),round(path(k,3)))=1;
                % make a thicker path for plotting purposes
                Forgr_c_thick(round(path(k,1)),round(path(k,2)),round(path(k,3)))=1;
                Forgr_c_thick(round(path(k,1))+1,round(path(k,2)),round(path(k,3)))=1;
                Forgr_c_thick(round(path(k,1))-1,round(path(k,2)),round(path(k,3)))=1;
                Forgr_c_thick(round(path(k,1)),round(path(k,2))+1,round(path(k,3)))=1;
                Forgr_c_thick(round(path(k,1)),round(path(k,2))-1,round(path(k,3)))=1;
                Forgr_c_thick(round(path(k,1)),round(path(k,2)),round(path(k,3))+1)=1;
                Forgr_c_thick(round(path(k,1)),round(path(k,2)),round(path(k,3))-1)=1;
            end
        end
    end
end

save ShortestPath.mat ShortestPath
save Forgr_c.mat Forgr_c
save Forgr_c_thick.mat Forgr_c_thick
% figure; 
% isosurface((BW3>0),0.1);
% hold on
% isosurface(Forgr_c,0.1)

end
