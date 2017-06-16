%% This function converts combines 2D FLAIR images to single T1rho image

function out_img = T1rho(in_img,slice_num)
[sx sy sz] = size(in_img(1,1).img);
time = 0:20:100;
time(1) = 1;

%% create image mask for ,threshold set at 25
mask = zeros(sx, sy);
mask(in_img(1,1).img(:,:,slice_num) >= 25) = 1; %thresholding
mask_tmp = mask;
%find biggest connected component
CC = bwconncomp(mask);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
mask_tmp(CC.PixelIdxList{idx}) = 0;
mask = mask - mask_tmp;
clear mask_tmp
%edge smoothing
mask = imfill(mask);
E = edge(mask,'canny');
%Dilate the edges
Ed = imdilate(E,strel('disk',2));
%Filtered image
Ifilt = imfilter(mask,fspecial('gaussian'));
%Use Ed as logical index into I to and replace with Ifilt
mask(Ed) = Ifilt(Ed);


%% T1rho fitting
fprintf(['T1rho mapping in progress... (' num2str(slice_num) '/' num2str(sz) ') \n']);
out_img = zeros(sx,sy);
h1 = figure('name',['Layer ' num2str(slice_num)]);
set(h1,'Position',[664 51 560 865]);
%hwait = waitbar(0,['T1rho mapping in progress... (' num2str(slice_num) '/' num2str(sz) ')']);
%set(hwait,'Position',[467.25 676.5 270 56.25]);
%count = 0;
pix_num = sx*sy;

for i = 1:sx
    for j = 1:sy
        %waitbar(count/pix_num,hwait);
        if mask(i,j) == 1
            parfor k = 1:6
                intensity(k) = double(in_img(1,k).img(i,j,slice_num));
            end
            if mean(intensity) > 0
                y = log(intensity);
                y(isinf(y)) = 0;
                coef = polyfitB(time(2:6),y(2:6),1,y(1));
                fit_val = polyval(coef,time);
                out_img(i,j) = abs(1/coef(1));
                    figure(h1);
                    subplot(2,1,1);
                    hold off
                    scatter(time,y,'b*');
                    xlim([0 100]);
                    ylim([0 10]);
                    xlabel('Spin lock time (ms)');                     
                    ylabel('Log of magnetization intensity');
                    hold on
                    plot(time,fit_val,'k-.','LineWidth',1);
                    title(['(' num2str(i) ',' num2str(j) '), T1rho = ' num2str(round(out_img(i,j))) 'ms']);
                clear A coef
             clear intensity
            end
        end
        %count = count + 1;
    end
    subplot(2,1,2);
    T1_disp = imrotate(out_img,-90);
    T1_disp = flipdim(T1_disp ,2);  
    imshow(T1_disp/400);
end
%delete(hwait);
close all;
end