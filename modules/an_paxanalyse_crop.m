function crop_img = an_paxanalyse_crop( in_img,pix_list,bb)
%AN_PAXANALYSE_CROP Performs 3d cropping based on input pixel list and
%bounding box
%   Detailed explanation goes here
%%

% bb
cmc=zeros(size(in_img,1),size(in_img,2),class(in_img));
% crop_img=zeros(bb(4)+1,bb(3)+1,size(in_img,3));
% whos
% for count=1:size(in_img,3)
%     tmp=in_img(:,:,count);
%     cmc(pix_list)=tmp(pix_list);
%     crop_img(:,:,count)=imcrop(cmc,bb);
% end
%     

% for count=1:size(in_img,3)
%     tmp=in_img(:,:,count);
%     cmc(pix_list)=tmp(pix_list);
%      pro_crop_img=imcrop(in_img(:,:,count),bb);
%         
%         crop_img(:,:,count)=pro_crop_img-imopen(pro_crop_img,strel('disk',25));
% 
% end
%     



% whos
% for count=1:size(in_img,3)
%         pro_crop_img=imcrop(in_img(:,:,count),bb);
%         
%         crop_img(:,:,count)=pro_crop_img-imopen(pro_crop_img,strel('disk',25));
% end


for count=1:size(in_img,3)
        crop_img(:,:,count)=imcrop(in_img(:,:,count),bb);        
        
end

end

