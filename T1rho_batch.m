clc;
clear all;
close all;

path = 'D:\T1Rho\NIFIT\';
time = 0:20:100;
parfor i = 1:6
    nii(i) = load_untouch_nii([path num2str(time(i)) '.nii']);    
end

[sx sy sz] = size(nii(1,1).img);
out_img = nii(1,1);
out_img.img = zeros(sx,sy,sz);
parfor i = 0:sz-1
    i
    img_tmp(:,:,sz-i) = T1rho(nii,sz-i);
end
for i = 1:sz
    out_img.img(:,:,i) = img_tmp(:,:,i);
end
clear img_tmp
save([path 'T1rho.mat']);
save_untouch_nii(out_img,[path 'T1rho.nii']);
%system('shutdown -s');