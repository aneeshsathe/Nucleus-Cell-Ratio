function [ BS_orig_cell_stack] = an_noise_rem_paxanalyse( bw_sum_cell,orig_cell_stack )
%AN_NOISE_REM_PAXANALYSE Summary of this function goes here
%   Detailed explanation goes here
%%

sum_cell=bw_sum_cell;

thresh_cell=imcomplement(sum_cell);

sum_orig_cell=sum(double(orig_cell_stack),3);

noise_img=thresh_cell.*sum_orig_cell;

 noise_val=mean(mean(noise_img))/size(orig_cell_stack,3);

BS_orig_cell_stack=orig_cell_stack-noise_val;






%%

end

