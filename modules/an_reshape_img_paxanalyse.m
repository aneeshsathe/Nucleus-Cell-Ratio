function out_img= an_reshape_img_paxanalyse( in_img )
%AN_RESHAPE_IMG_PAXANALYSE Summary of this function goes here
%   Detailed explanation goes here

out_img=reshape(in_img,[],size(in_img,2)*size(in_img,3));

end

