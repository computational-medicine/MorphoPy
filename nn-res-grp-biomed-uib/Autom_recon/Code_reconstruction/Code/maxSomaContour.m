function somaContour = maxSomaContour(Dth, maxSomaSizeZ, SmoothLength,fname)
% maxSomaContour.m
% Finds soma contours around a predetermined point in the segmented image
% stack
% BJZ: 16-04-29
% Assumes x and y voxelsize are the same.

Niter = 10;
Npntscontour = 100;
doplotslices=false; % for debugging purposes

imstack_segm=cbiReadNifti(fname);
load zsnittmax
load voxel_size
load somapkt
load t1; load t3; load t5
cntr = somapkt([2,1]);

x2 = cell(0);
y2 = cell(0);
Nslices = ceil(maxSomaSizeZ/2/(voxel_size(3)*1e3));

Nslicerange = min(max([zsnittmax-Nslices, zsnittmax+Nslices],1),size(imstack_segm,3));
slice = Nslicerange(1):Nslicerange(end);

for n = 1:length(slice);
    im = imstack_segm(:,:,slice(n));
    [x2{n},y2{n},A(n)] = findsoma2D(im, cntr,Dth,Niter,Npntscontour);
    
    
    if doplotslices
    figure(n+100)
    imagesc(im); hold on
    plot(cntr(1),cntr(2),'rx')
    plot(mean(x2{n}),mean(y2{n}),'rd')
    plot(x2{n},y2{n},'g')
    end
    
end
Asmooth = smooth(A,SmoothLength*1e-3/voxel_size(3),'sgolay');
[maxA,maxind] = max(Asmooth);
maxareaX = x2{maxind};
maxareaY = y2{maxind};
maxslice = slice(maxind);

figure; plot(slice*voxel_size(3)*1e3,A); hold on
plot(slice*voxel_size(3)*1e3,Asmooth);
plot(slice(maxind)*voxel_size(3)*1e3,maxA,'rd')
title('Soma cross-section area'); ylabel('pixels');
xlabel('z (um)')
figure; imagesc(imstack_segm(:,:,slice(maxind)))
hold on;
plot(maxareaX,maxareaY,'g')
plot(mean(maxareaX),mean(maxareaY),'rx')

% convert to real world coordinates

X = (maxareaX+t1-1)*voxel_size(1)*1e3;
Y = -(maxareaY+t3-1)*voxel_size(2)*1e3;
Z = -ones(size(X))*(maxslice+t5-1)*voxel_size(3)*1e3;

somaContour.X = X';
somaContour.Y = Y';
somaContour.Z = Z';
save somaContour.mat somaContour

end

function [x2,y2,A] = findsoma2D(im, cntr,Dth,Niter, Npntscontour)

Dth = Dth*pi/180;

% add empty border around image
imtemp = zeros(size(im)+2);
imtemp(2:end-1,2:end-1) = im;
im = imtemp;
% and move center point accordingly
cntr = cntr+1; % x+1 & y+1

im = bwselect(im,cntr(1),cntr(2),8);

c = contourc(double(im),1);
if isempty(c)
    x2=[]; y2=[];A=0; return
else
    [x,y] = C2xyz(c);
    n = 1;
    while ~inpolygon(cntr(1),cntr(2),x{n},y{n})
        n = n+1;
    end
    x = x{n}; y = y{n};
    
    c = [x;y];
    %c(1,:) = smooth(c(1,:),15);
    %c(2,:) = smooth(c(2,:),15);
    %hold on; plot(c(1,:),c(2,:),'rd-')
    
    %curv = LineCurvature2D(c');
    %figure(2);
    %scatter(c(1,:),c(2,:),[],curv); hold on
    
    
    [x,y,th,r,x2,y2] = findxy(c,cntr,Dth, Npntscontour,false);
    % repeat with new center
    for i=1:Niter-1;
        [x,y,th,r,x2,y2] = findxy(c,[mean(x2),mean(y2)],Dth, Npntscontour,false);
    end
    %[x,y,th,r,x2,y2,th2,r2] = findxy(c,[mean(x2),mean(y2)],Dth,true);
        
    A = polyarea(x2,y2);
    % correct for the 1 pixel white border that was added
    x2=x2-1;
    y2=y2-1;
end
end

function [x,y,th,r,x2,y2] = findxy(c,cntr,Dth, Npntscontour, doplot)
x = c(1,:); y = c(2,:);
r = sqrt((x-cntr(1)).^2 + (y-cntr(2)).^2);
th = mod(atan2(x-cntr(1), y-cntr(2))+pi,2*pi)-pi;

th2 = -pi:2*pi/Npntscontour:pi;
%th2a = zeros(size(th));
usepoints = false(size(r));
for i=1:length(th2)
    ind = find(abs(mod(th-th2(i)+pi,2*pi)-pi) < Dth/2);
    [~, ind2] = min(r(ind));
    usepoints(ind(ind2))=true;
    %th2a(i) = th(ind(ind2));
    
    %testind(i) = ind(ind2);
end

th2b = th(usepoints);
r2b = r(usepoints);
[th2b,sortind] = sort(th2b);
r2b = r2b(sortind);

r2 = interp1([th2b(end)-2*pi, th2b, th2b(1)+2*pi],[r2b(end) r2b r2b(1)],th2,'linear','extrap');
x2 = cntr(1)+r2.*sin(th2);
y2 = cntr(2)+r2.*cos(th2);

if doplot
hold on
figure(200)
plot(th2,r2,'r--'); hold on
plot(th,r,'.')
figure(201); hold on
plot(x2,y2,'g')
plot(cntr(1),cntr(2),'rx')
end

end



function [x,y,z] = C2xyz(C)
% C2XYZ returns the x and y coordinates of contours in a contour
% matrix and their corresponding z values. C is the contour matrix given by 
% the contour function. 
% 
%
%% Syntax
% 
%  [x,y] = C2xyz(C)
%  [x,y,z] = C2xyz(C)
% 
%% Description 
% 
% [x,y] = C2xyz(C) returns x and y coordinates of contours in a contour
% matrix C
% 
% [x,y,z] = C2xyz(C) also returns corresponding z values. 
% 
% 
%% Example
% Given a contour plot, you want to know the (x,y) coordinates of the contours, 
% as well as the z value corresponding to each contour line. 
%
% C = contour(peaks); 
% [x,y,z] = C2xyz(C);
% 
% This returns 1 x numberOfContourLines cells of x values and y values, and
% their corresponding z values are given in a 1 x numberOfContourLines
% array. If you'd like to plot a heavy black line along all of the z=0
% contours and a dotted red line along the z = -2 contours, try this: 
% 
% hold on; % Allows plotting atop the preexisting peaks plot. 
% for n = find(z==0); % only loop through the z = 0 values. 
%     plot(x{n},y{n},'k','linewidth',2)
% end
% 
% for n = find(z==-2) % now loop through the z = -2 values. 
%     plot(x{n},y{n},'r:','linewidth',2)
% end
% 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Created by Chad Greene, August 2013. 
% 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% See also contour, contourf, clabel, contour3, and C2xy.

x=cell(0);y=cell(0);z=cell(0);

m(1)=1; 
n=1;  
try
    while n<length(C)
        n=n+1;
        m(n) = m(n-1)+C(2,m(n-1))+1; 
        
    end
end

for nn = 1:n-2
     x{nn} = C(1,m(nn)+1:m(nn+1)-1); 
     y{nn} = C(2,m(nn)+1:m(nn+1)-1); 
     if nargout==3
        z(nn) = C(1,m(nn));
     end
end

end
