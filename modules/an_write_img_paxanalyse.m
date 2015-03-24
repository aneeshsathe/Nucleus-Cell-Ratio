function []=an_write_img_paxanalyse(out_img,max_img,proj_cell_nuc,img_write_path,series_count,obj_count,t_count)
%AN_WRITE_IMG_PAXANALYSE Summary of this function goes here
%   Detailed explanation goes here
% 
% 
% out_img=[reshape(crop_cell,[],size(cell_bw,2)*size(cell_bw,3));
%     reshape(cell_bw,[],size(cell_bw,2)*size(cell_bw,3));
%     reshape(nuc_bw,[],size(cell_bw,2)*size(cell_bw,3));
%     reshape(crop_nuc,[],size(cell_bw,2)*size(cell_bw,3))%;
%     %reshape(phos_stack(cell_pix_list),[],size(cell_bw,2)*size(cell_bw,3));
%     ];
% max_img=max(orig_cell_stack,[],3);
% %sum_img=sum(orig_cell_stack,3,'double');
% %         imtool(out_img,[])




if t_count==1
    
    if isa(out_img,'gpuArray')
        imwrite(double(gather(out_img)),[img_write_path, 'Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
        'tiff');%,'Compression', 'jpeg', 'RowsPerStrip', 16);
    
    imwrite(double(gather(max_img)),[img_write_path, 'Max_Proj_Series',num2str(series_count),'_cell_', num2str(obj_count), '.tiff'],...
        'tiff');
    imwrite(double(gather(proj_cell_nuc)),[img_write_path, 'Proj_BW',num2str(series_count),'_cell_', num2str(obj_count), '.tiff'],...
        'tiff');
    else
    imwrite(out_img,[img_write_path, 'Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
        'tiff');%,'Compression', 'jpeg', 'RowsPerStrip', 16);
    
    imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_cell_', num2str(obj_count), '.tiff'],...
        'tiff');
    
    imwrite(proj_cell_nuc,[img_write_path, 'Proj_BW',num2str(series_count),'_cell_', num2str(obj_count), '.tiff'],...
        'tiff');
    % imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
    % 'tiff');
    end
elseif t_count>1
    imwrite(out_img,[img_write_path, 'Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
        'tiff','WriteMode','append');%, 'Compression', 'jpeg', 'RowsPerStrip', 16);
    
    imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_', '_cell_', num2str(obj_count), '.tiff'],...
        'tiff','WriteMode','append');
    
    imwrite(proj_cell_nuc,[img_write_path, 'Proj_BW',num2str(series_count),'_', '_cell_', num2str(obj_count), '.tiff'],...
        'tiff','WriteMode','append');
    
    %imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
    % 'tiff','WriteMode','append');
    
end


end

