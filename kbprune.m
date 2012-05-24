function [regions] = kbprune(candidates, K, v_th)
%KBPRUNE Remove weak regions from Kadir-Brady candidates

% Altered 8/2/2008:
%   Now we sort the gamma values in temp before we get the distances,
%   instead of the other way around

% apply a global threshold to the gamma vals
thresh_val = .6 * max(candidates.gamma);
indices = find(candidates.gamma > thresh_val);

temp.gamma = candidates.gamma(indices);
temp.scale = candidates.scale(indices);
temp.r = candidates.r(indices);
temp.c = candidates.c(indices);

% sort the gamma values, and order everything by that
[sgamma, sidx] = sort(temp.gamma, 'descend');
temp.gamma = sgamma;
temp.scale = temp.scale(sidx);
temp.r = temp.r(sidx);
temp.c = temp.c(sidx);

% create a distance matrix
n = length(temp.gamma);

if K+1 > n
    regions.gamma = temp.gamma(1);
    regions.r = temp.r(1);
    regions.c = temp.c(1);
    regions.scale = temp.scale(1);
    return;
end

D = zeros(n,n);

pts = [temp.c, temp.r, temp.scale];

% fill it with distances
for i=1:n
    pt = pts(i,:);

    dists = sqrt(sum( (pts-repmat(pt,size(pts,1),1)).^2,2 ) );

    D(i,:) = dists';
    D(:,i) = dists;

%     for j=i+1:n
%         D(i,j) = (temp.scale(i)-temp.scale(j))^2 + ...
%                  (temp.r(i)-temp.r(j))^2 + ...
%                  (temp.c(i)-temp.c(j))^2;
%         D(j,i) = D(i,j);
%     end
end

% D = sqrt(D);

nReg = 0;
regions = [];
pos = zeros(3,K+1);

% now do the pruning process
% [Mx,Midx] = sort(temp.gamma, 'descend');
for i=1:n
    % index = Midx(i);
    index = i;

    pos(1,1) = temp.c(index);
    pos(2,1) = temp.r(index);
    pos(3,1) = temp.scale(index);

    [sD,sidx] = sort(D(index,:));

    for j=1:K
        pos(1,j+1) = temp.c(sidx(j+1));
        pos(2,j+1) = temp.r(sidx(j+1));
        pos(3,j+1) = temp.scale(sidx(j+1));
    end

    cent = mean(pos,2);

    v = var(sqrt(sum((pos-repmat(cent,1,K+1)).^2,1)));

    if (v > v_th)
        continue;
    end

    % now that we know the regions is "suffiently clustered", make sure the
    % region is "far enough" from already clustered regions
    if (nReg > 0)
        d = sqrt(sum(([regions.c, regions.r, regions.scale] - ...
                 repmat(cent', nReg, 1)).^2,2));
        if (sum(cent(3) >= d) == 0)
            nReg = nReg+1;
            regions.c(nReg,1) = cent(1);
            regions.r(nReg,1) = cent(2);
            regions.scale(nReg,1) = cent(3);
            regions.gamma(nReg,1) = temp.gamma(index);
        end
    else
        nReg = nReg+1;
        regions.c(nReg,1) = cent(1);
        regions.r(nReg,1) = cent(2);
        regions.scale(nReg,1) = cent(3);
        regions.gamma(nReg,1) = temp.gamma(index);
    end
end

return;