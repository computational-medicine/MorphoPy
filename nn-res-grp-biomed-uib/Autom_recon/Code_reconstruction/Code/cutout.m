function mask = cutout(z,txt)
    projection = max(z,[],3);
    fig = figure(198);
    imagesc(projection);
    title(txt)
    [x,y] = getline(fig);
    [X,Y] = meshgrid(1:size(z,2),1:size(z,1));
    mask = inpolygon(X,Y,x,y);
    close(fig);
    figure('Name','Masked image stack'); imagesc(max(double(z+max(z(:)/10)).*repmat(double(mask),1,1,size(z,3)),[],3))

end