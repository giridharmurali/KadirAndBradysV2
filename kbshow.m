function kbshow(im,regions)
%KBSHOW Display regions found by the Kadir-Brady algorithm

figure(1)
subplot('position', [0 0 1 1]);
imagesc(im);
colormap gray;
hold on;

theta = [-.03:.01:2*pi];

for i=1:length(regions.gamma)
    r = regions.scale(i)/2;

    x = r*cos(theta);
    y = r*sin(theta);

    X = x+regions.c(i);
    Y = y+regions.r(i);
    
    j=num2str(i);

    plot(X,Y,'r-');
    text(regions.c(i),regions.r(i),j,'Fontsize',18,'Color',[1 1 0]);
end


% Plot six regions with highest saliency
figure(2)
subplot('position', [0 0 1 1]);
imagesc(im);
colormap gray;
hold on;

theta = [-.03:.01:2*pi];

for i=1:6
    r = regions.scale(i)/2;

    x = r*cos(theta);
    y = r*sin(theta);

    X = x+regions.c(i);
    Y = y+regions.r(i);
    
    j=num2str(i);

    plot(X,Y,'r-');
    text(regions.c(i),regions.r(i),j,'Fontsize',18,'Color',[1 1 0]);
end


index=1;

figure(3)
for i=1:6
    xmin=regions.c(i)-(regions.scale(i)/2);
    ymin=regions.r(i)-(regions.scale(i)/2);
    height=regions.scale(i);
    width=regions.scale(i);
    cropped_image=imcrop(im,[xmin ymin width height]);
    resized_image=imresize(cropped_image, [11 11]);
    subplot(1,6,i), imshow(resized_image);
    index=index+1;
end


return;