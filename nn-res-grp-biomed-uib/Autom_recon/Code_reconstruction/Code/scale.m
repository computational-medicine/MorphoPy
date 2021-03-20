% SCALE(IM) Scaling.
% SCALE(IM) scales the image IM between 0 and 1.
% SCALE(IM,LOW,HIGH) scales the image IM between LOW and HIGH.
%
function [im] = scale(varargin)

im = varargin{1};
if nargin == 1
    low = 0;
    high = 1;
else
    low = varargin{2};
    high = varargin{3};
end;

 if isempty(find(im))
    disp('Empty matrix, no scaling.')
    return;
end;

% must not be uint16 values!!!
im = double(im);

% scale to 0 and 1
im = im - min(im(:));

maxIm = max(im(:));
if ne(maxIm,0) | ~isempty(find(im))  
    im = im / maxIm;
end;

% scale to high and low
im = im * (high-low);
im = im + low;

