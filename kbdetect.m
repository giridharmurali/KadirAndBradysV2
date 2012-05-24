function [regions] = kbdetect(in_im, windows, mask)
%KBDETECT Compute Kadir-Brady salience for all pixels

if size(in_im,3) ~= 1
    im = rgb2gray(in_im);
else
    im = in_im;
end

nr = size(im,1);
nc = size(im,2);

% find pixels that we are going to examine
[r,c] = find(mask);
nPix = length(r);

% get how many scales we are doing
nScales = length(windows);

USE_INTENSITY = 1;

if USE_INTENSITY % setting this to 1 will use intensity as the salience feature
    QUANTIZATION = 16;
    edges = [0:QUANTIZATION:256];
else % here we use 
    QUANTIZATION = 2*pi / 8;
    edges = [-pi:QUANTIZATION:pi];

    % compute derivatives and then the orientation image
    dx = imfilter(im, [-1 0 1], 'replicate');
    dy = imfilter(im, [-1 0 1]', 'replicate');

    oriIm = atan2(double(dy),double(dx));
end

last_h = zeros(length(edges)-1,nPix);
ss_x = zeros(nScales, nPix);
entropy = zeros(nScales, nPix);

out_scales = nScales;

% now iterate
for s_count=1:nScales
    win_size = windows(s_count);
    this_win = floor( (win_size)/2 );

    if (this_win+1 > nr/2 || this_win+1 > nc/2)
        out_scales = s_count-1;
        break;
    end

    for i=1:nPix
        % ignore points too close to the edge of the image
        if (s_count == 2 && r(i) == 93 && c(i) == 93)
            dummy = 1;
        end

        min_r = r(i)-this_win;
        min_c = c(i)-this_win;
        max_r = min_r + win_size-1;
        max_c = min_c + win_size-1;

        if (min_r < 1)
            min_r = 1;
        end

        if (max_r > nr)
            max_r = nr;
        end

        if (min_c < 1)
            min_c = 1;
        end

        if (max_c > nc)
            max_c = nc;
        end

%         if (min_r < 1 || max_r > nr || min_c < 1 || max_c > nc)
%             continue;
%         end

        % compute the histogram of intensity values in this region
        if USE_INTENSITY
            patch = im(min_r:max_r, min_c:max_c);
        else
            patch = oriIm(min_r:max_r, min_c:max_c);
        end

        h = histc(patch(:), edges);
        h = h(1:end-1);
        h = h/sum(h);

        % save it
%         scale_space(s_count,i,:) = h;
        idx = find(h > 0);
        entropy(s_count,i) = -sum( h(idx).*log(h(idx)) );

        if (s_count >= 2)
            dif = abs(h-last_h(:,i));
            factor = windows(s_count)^2/(2*windows(s_count)-1);
            factor1 = (windows(s_count)+1)^2/(2*windows(s_count)+1);
            ss_x(s_count,i) = factor * sum(dif);
            ss_x_new(s_count,i)=factor1*sum(dif);

            if (s_count == 2)
                ss_x(s_count-1,i) = ss_x(s_count,i);
                ss_x_new(s_count-1,i) = ss_x_new(s_count,i); % New modification - it can also be zero as it is not going to matter
            end
        end

        last_h(:,i) = h;
    end
end

% now find local maxima in scale space by looking at the second derivative
% being less than zero. calculate the weights and smooth since the first
% derivative calculation will be noisy
% fx = [-1 1 0]';
fxx = [1 -2 1]';

% ss_x = imfilter(scale_space(1:out_scales,:,:), fx, 'replicate');
% ss_x = sum(abs(ss_x),3);
% factor_vec = (windows(1:out_scales).^2)./(2*windows(1:out_scales)-1);
% factor = repmat(factor_vec',1,nPix);
% ss_x = ss_x.*factor;
% weight = imfilter(ss_x(1:out_scales,:), [1/3,1/3,1/3], 'replicate');
weight = ss_x(1:out_scales,:);
weight1= ss_x_new(1:out_scales,:);
% New modification - change the weight 
temp=weight1(:,2:end);
temp(:,end+1)=zeros(size(temp,1),1);
new_weight = (weight+temp)/2;
new_weight(:,end)=new_weight(:,end-1);

% idx = find(scale_space == 0);
% scale_space(idx) = 1;
% 
% entropy = -sum(scale_space(1:out_scales,:,:) .* ...
%     log(scale_space(1:out_scales,:,:)),3);
ss_xx = imfilter(entropy(1:out_scales,:), fxx, 'replicate');
ss_xx(1,:) = 0;
[int_pts] = find(ss_xx < 0);

% weight these points by the weighting function, which is the scale window
% size times the first derivative of the scale-space function

regions.gamma = entropy(int_pts) .* weight(int_pts); % New modification - change 'weights' to 'new_weights'
[scales,locs] = ind2sub([out_scales nPix], int_pts);

regions.scale = windows(scales)';
regions.r = r(locs);
regions.c = c(locs);

% for i=1:length(locs)
%     regions.gamma(i,1) = entropy(scales(i),locs(i)) * weight(scales(i),locs(i));
%     regions.scale(i,1) = windows(scales(i));
%     regions.r(i,1) = r(locs(i));
%     regions.c(i,1) = c(locs(i));
% end

return;