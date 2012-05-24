function [regions] = kbextract(im,windows,K,v_th,mask)
%KBEXTRACT Extract Kadir-Brady salience regions from an image
tic
im=imread('cars\image_0006.jpg');
s1=size(im,1); s2=size(im,2); 
min_size=min(s1,s2);
if min_size==s1
    temp1=100/s1; temp2=ceil(s2*temp1);
    im=imresize(im, [100 temp2]);
else
    temp1=100/s2; temp2=ceil(s1*temp1);
    im=imresize(im, [temp2 100]);
end
    
%f=fspecial('gaussian', [3 3], 8);
%im=imfilter(im,f);

nr = size(im,1);
nc = size(im,2);


if nargin < 5
    mask = ones(nr,nc);
end

if nargin < 4
    v_th = 1;
end

if nargin < 3
    K = 5;
end

if nargin < 2
    windows = [10:2:25];
end

% the first step is to find local maxima in scale-space
base_regions = kbdetect(im,windows,mask);

% now it clustering must be done
regions = kbprune(base_regions, K, v_th);

kbshow(im,regions);

regions.gamma

toc
return;