tStart = tic;
clc;
clear all;
close all;

path = 'D:\T1Rho\NIFIT\';
time = 0:20:100;
for i = 1:6
    nii(i) = load_untouch_nii([path num2str(time(i)) '.nii']);    
end
time(1) = 1;

[sx sy sz] = size(nii(1,1).img);
test_num = 33; %slice number
mask = zeros(sx, sy);
mask(nii(1,1).img(:,:,test_num) >= 25) = 1;
mask_tmp = mask;
CC = bwconncomp(mask);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
mask_tmp(CC.PixelIdxList{idx}) = 0;
mask = mask - mask_tmp;
%mask = imerode(mask,strel('disk',2));
mask = imfill(mask);
%mask = imdilate(mask,strel('disk',2));
clear mask_tmp

T1 = zeros(sx,sy);
h1 = figure;
set(h1,'Position',[44 61 560 865]);
%hwait = waitbar(0,'T1rho mapping in progress...');
%set(hwait,'Position',[467.25 676.5 270 56.25]);
%count = 0;
pix_num = sx*sy;
for i = 1:sx
    for j = 1:sy
        %waitbar(count/pix_num,hwait);
        if mask(i,j) == 1
            parfor k = 1:6
                intensity(k) = double(nii(1,k).img(i,j,test_num));
            end
            if mean(intensity) > 0 %&& mean(intensity) <= 1500 
                y = log(intensity);
                y(isinf(y)) = 0;
                coef = polyfitB(time(2:6),y(2:6),1,y(1));
                fit_val = polyval(coef,time);
                T1(i,j) = abs(1/coef(1));
                %if coef(1) < 0
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
                    title(['(' num2str(i) ',' num2str(j) '), T1rho = ' num2str(round(T1(i,j))) 'ms']);
                %end
                %if coef(1) < 0
                %    T1(i,j) = -1/coef(1);
                %end
                %{
                if abs(coef(1)) <= 0.05 %&& coef(3) <= 100
                    coef
                    A = nlinfit(time,intensity,F,x0); 
                    if A(2) < 500 && A(2) > 0;
                        figure(h)
                        hold off
                        scatter(time,intensity,'b*');
                        xlim([0 100]);
                        ylim([0 100]);
                        xlabel('Spin lock time (ms)');   
                        ylabel('Magnetization intensity');
                        hold on;
                        fit_val = A(1)*(exp(-time./A(2)));
                        plot(time,fit_val,'k-.','LineWidth',1);
                        T1(i,j) = A(2);              
                        title(['(' num2str(i) ',' num2str(j) '), T1rho = ' num2str(round(T1(i,j))) 'ms, M0 = ' num2str(round(A(1)))]);    
                        pause(0.005);
                    end
                %}
                %end
                clear A coef
             clear intensity
            end
        end
        %count = count + 1;
    end
    subplot(2,1,2);
    T1_disp = imrotate(T1,-90);
    T1_disp = flipdim(T1_disp ,2);  
    imshow(T1_disp/400);
end
delete(hwait);
h3 = figure;
imshow(T1_disp/400);
saveas(h3,'D:\T1Rho\T1rho\T1r_new.jpg','jpg');
close all;
tEnd = toc;
fprintf('Total time for mapping = %d minutes and %f seconds \n',floor(tEnd/60),rem(tEnd,60));