% BWSIZE Finds the size of binary segments
%
% [DIM,FASER] = BWSIZE(BW,CONN) finds the size of the disconnected binary
% segments inside BW and returns the number of pixels in DIM and the
% piecewise constant image FASER connected to the label in DIM
%
function [dim,faser] = bwsize(BW,conn)

%[faser,L] = bwlabeln(BW,conn);
[faser,L] = bwlabeln(BW);

dim = zeros(L,1);
for i = 1 : L
    % this region
    reg = eq(faser,i);
    
    % number of pixels
    dim(i) = length(find(reg));
    
end;
