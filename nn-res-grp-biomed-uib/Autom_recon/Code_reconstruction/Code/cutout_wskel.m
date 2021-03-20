function mask = cutout_wskel(im,txt,skel)
    if isempty(txt)
        txt = 'please click boundaries of area containing the neuron; press enter to finish, backspace to correct';
    end
    projection = max(im,[],3)/max(im(:));
    projskel = max(skel,[],3);
    projection(projskel>0) = 2;
    fig = figure(198);
    imagesc(projection);
    title(txt)
    [x,y] = getline(fig);
    close(fig)
    if ~isempty(x)
        [X,Y] = meshgrid(1:size(im,2),1:size(im,1));
        mask = inpolygon(X,Y,x,y);
    else
        mask = zeros(size(projection));
    end
    figure; imagesc(max(double(im+max(im(:)/10)).*repmat(double(mask),1,1,size(im,3)),[],3))
end