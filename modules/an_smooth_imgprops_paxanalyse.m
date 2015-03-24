function [bw_re_prop,bw_sum_cell ]= an_smooth_imgprops_paxanalyse(cell_stack)
%AN_SMOOTH_IMGPROPS_PAXANALYSE Summary of this function goes here
%   Detailed explanation goes here

%%
smooth_cell=smooth3(smooth3(cell_stack,'gaussian'),'gaussian');



sum_smooth_cell=sum(smooth_cell,3);

h = fspecial('gaussian', [3 3],0.5);
sum_smooth_cell=filter2(h, sum_smooth_cell);


%% crop cell

% bs_sum_smooth_cell=sum_smooth_cell-imopen(sum_smooth_cell, strel('disk',250)); %% DEFAULT

bs_sum_smooth_cell=imopen(sum_smooth_cell, strel('disk',5));
%bs_sum_smooth_cell=bs_sum_smooth_cell-imopen(bs_sum_smooth_cell, strel('disk',250));
% bs_sum_smooth_cell=sum_smooth_cell-imopen(sum_smooth_cell, strel('ball',250,3));
%bs_sum_smooth_cell=bs_sum_smooth_cell-imopen(bs_sum_smooth_cell, strel('ball',15,5));
bw_sum_cell=bs_sum_smooth_cell>graythresh(bs_sum_smooth_cell);
bw_sum_cell=imfill(bw_sum_cell,'holes');

%         bw_sum_cell=bs_sum_smooth_cell>graythresh(sum_cell);
%         bw_sum_cell=sum_smooth_cell>graythresh(sum_cell);
%         bw_sum_cell=imclearborder(bw_sum_cell);



% se2=strel('disk',2);
% bw_sum_cell=imclose(bw_sum_cell,se2);
% bw_sum_cell=imdilate(bw_sum_cell,se2);
% imshow(bw_sum_cell)
bw_re_prop=regionprops(bw_sum_cell,'Area','BoundingBox','PixelIdxList','PixelList');

end

