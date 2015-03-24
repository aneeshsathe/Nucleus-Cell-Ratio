function [out_img, thresh_fac]= an_gray(I)
%AN_GRAY returns graythresholded img
%   Detailed explanation goes here
thresh_fac=graythresh(I);
out_img=I>thresh_fac;

end

