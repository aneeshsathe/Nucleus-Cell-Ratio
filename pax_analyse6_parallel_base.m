function [ xls_label ] = pax_ana_6_par_test(root,img_write_path,file_list, file_count,series_to_read,preset_vals )
%PAX_ANA_6_PAR_TEST Summary of this function goes here
%   Detailed explanation goes here
[file_ext,cell_chanl,nuc_chanl,phos_chanl,cell_2d_area_thresh,cell_3d_area_thresh,nuc_3d_area_thresh]=deal(preset_vals{1:7});
file_path=[root, file_list(file_count).name]; %modification
        %% initialize file reader and get number of series
%         reader = an_gen_reader( file_path,1 );
        [reader,series_num,pro_series_count ]= an_gen_reader(file_path,series_to_read);
        
        
        %%
        for series_count=pro_series_count
            %for series_count=1:series_num%use this for all series
            %for series_count=[2,6,10 ]%use this to specify series in [series1,series2, series18, etc]
            reader.setSeries(series_count-1)
            C_num=reader.getSizeC;
            Z_num=reader.getSizeZ;
            T_num=reader.getSizeT;
            
            for t_count=1:T_num
                %% get cell and nuc channels
                %get cell
                counter=0;
                %[orig_cell_stack, orig_nuc_stack]=deal(zeros(512,512,Z_num));
                
                z_range=1:Z_num;
                
                orig_cell_stack=bfGetStack(reader,cell_chanl,z_range,t_count-1);
                orig_nuc_stack=bfGetStack(reader,nuc_chanl,z_range,t_count-1);
                orig_phos_stack=bfGetStack(reader,phos_chanl,z_range,t_count-1);
                
                cell_stack=mat2gray(orig_cell_stack);
                nuc_stack=mat2gray(orig_nuc_stack);
                phos_stack=mat2gray(orig_phos_stack);
                %             smooth_cell=smooth3(cell_stack);
                %             smooth_nuc=smooth3(nuc_stack);
                smooth_cell=smooth3(smooth3(cell_stack,'gaussian'),'gaussian');
                smooth_nuc=smooth3(smooth3(nuc_stack,'gaussian'),'gaussian');
                smooth_phos=smooth3(smooth3(phos_stack,'gaussian'),'gaussian');
                
                sum_cell=sum(cell_stack,3)+sum(phos_stack,3);
               
                sum_smooth_cell=sum(smooth_cell,3)+sum(smooth_phos,3);
                re_smooth_cell=reshape(smooth_cell,[],size(smooth_cell,2)*size(smooth_cell,3));
                
                %             figure, imshow(sum_cell,[])
                h = fspecial('gaussian', [3 3],0.5);
                sum_smooth_cell=filter2(h, sum_smooth_cell);
                
                
                %% crop cell
                bs_sum_smooth_cell=sum_smooth_cell-imopen(sum_smooth_cell, strel('disk',250));
                %bs_sum_smooth_cell=bs_sum_smooth_cell-imopen(bs_sum_smooth_cell, strel('disk',250));
                % bs_sum_smooth_cell=sum_smooth_cell-imopen(sum_smooth_cell, strel('ball',250,3));
                %bs_sum_smooth_cell=bs_sum_smooth_cell-imopen(bs_sum_smooth_cell, strel('ball',15,5));
                bw_sum_cell=bs_sum_smooth_cell>graythresh(bs_sum_smooth_cell);
                
                
                %         bw_sum_cell=bs_sum_smooth_cell>graythresh(sum_cell);
                %         bw_sum_cell=sum_smooth_cell>graythresh(sum_cell);
                %         bw_sum_cell=imclearborder(bw_sum_cell);
                
                se2=strel('disk',1);
                bw_sum_cell=imclose(bw_sum_cell,se2);
                bw_sum_cell=imdilate(bw_sum_cell,se2);
                
                bw_re_prop=regionprops(bw_sum_cell,'Area','BoundingBox','PixelIdxList','PixelList');
                bw_re_prop=regionprops(bw_sum_cell,'Area','BoundingBox','PixelIdxList','PixelList');
                
                if ~isempty(bw_re_prop)
                    pix_list=bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).PixelIdxList;
                    crop_bw=imcrop(bw_sum_cell,bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).BoundingBox);
                    bb=bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).BoundingBox;
                    
                    obj_idx=find([bw_re_prop.Area]>cell_2d_area_thresh);
                    
                    for obj_count=1:size(obj_idx,2)
                        
                        bb=bw_re_prop(obj_idx(obj_count)).BoundingBox;
                        pix_list=bw_re_prop(obj_idx(obj_count)).PixelIdxList;
                        [tmp, cmc, cc, cn, cp, cop]=deal(zeros(size(cell_stack,1),size(cell_stack,2)));
                        [coc, noc]=deal(zeros(size(cell_stack,1),size(cell_stack,2),'uint16'));
                        
                        %crop images
                        for count=1:size(smooth_cell,3)
                            tmp=smooth_cell(:,:,count);
                            %                         cmc=tmp;
                            cmc(pix_list)=tmp(pix_list);
                            crop_smo_cell(:,:,count)=imcrop(cmc,bb);
                            
                            tmp=cell_stack(:,:,count)+phos_stack(:,:,count);
                            %                         cc=tmp;
                            cc(pix_list)=tmp(pix_list);
                            crop_cell(:,:,count)=imcrop(cc,bb);
                            
                            tmp=nuc_stack(:,:,count);
                            cn=tmp;
                            %                         cn(pix_list)=tmp(pix_list);
                            crop_nuc(:,:,count)=imcrop(cn,bb);
                            
                            %             tmp=phos_stack(:,:,count);
                            %             cp(pix_list)=tmp(pix_list);
                            %             crop_phos(:,:,count)=imcrop(cc,bb);
                            
                            tmp=orig_cell_stack(:,:,count);
                            %                         coc=tmp;
                            coc(pix_list)=tmp(pix_list);
                            crop_orig_cell(:,:,count)=imcrop(coc,bb);
                            
                            tmp=orig_nuc_stack(:,:,count);
                            %                         coc=tmp;
                            noc(pix_list)=tmp(pix_list);
                            crop_orig_nuc(:,:,count)=imcrop(noc,bb);
                            
                            %             tmp=orig_phos_stack(:,:,count);
                            %             cop(pix_list)=tmp(pix_list);
                            %             crop_orig_phos(:,:,count)=imcrop(cc,bb);
                            
                            %             imshow(crop_smo_cell(:,:,count),[])
                        end
                        
                        
                        %% get pix list
                        se3=strel('ball', 55, 55);
                        crop_cell=imopen(crop_cell,se3);
                        crop_nuc=imopen(crop_nuc,se3);
                        cell_thresh_fac=graythresh(crop_cell);
                        nuc_thresh_fac=graythresh(crop_nuc);
                        
                        %apply thresh
                        %             unre_crop_smo_cell=reshape(re_crop_smo_cell,[],size(crop_smo_cell,2),size(crop_smo_cell,3));
                        %             cell_bw=unre_crop_smo_cell>cell_thresh_fac;
                        cell_bw=crop_cell>cell_thresh_fac;
                        nuc_bw=crop_nuc>nuc_thresh_fac;
                        
                        cell_bw=cell_bw|nuc_bw;
                        %fill holes
                        cell_bw=bwareaopen(imfill(imdilate(cell_bw,se2),'holes'),cell_3d_area_thresh);
                        nuc_bw=bwareaopen(imfill(nuc_bw,'holes'),nuc_3d_area_thresh);
                        % make new cell_bw which includes both cell and nuclear
                        
                        %get pix list
                        prop_cell_bw=regionprops(cell_bw,'Area', 'PixelIdxList','FilledImage');
                        prop_nuc_bw=regionprops(nuc_bw,'Area', 'PixelIdxList','BoundingBox','FilledImage');
                        
                        if (sum([prop_cell_bw.Area]>cell_3d_area_thresh)&&sum([prop_nuc_bw.Area]>nuc_3d_area_thresh))
                            %                         disp('passed')
                            
                            cell_pix_list=prop_cell_bw([prop_cell_bw.Area]>cell_3d_area_thresh).PixelIdxList;
                            
                            %to select only those nuclear pixels that are most
                            %in the cell pixels
                            rat_member=zeros(size(prop_nuc_bw,1));
                            for count=1:size(prop_nuc_bw,1)
                                nuc_pix_list=prop_nuc_bw(count).PixelIdxList;
                                nuc_s=size(nuc_pix_list,1);
                                cell_s=size(cell_pix_list,1);
                                member_pix=sum(ismember(nuc_pix_list,cell_pix_list));
                                % all(sum(ismember(nuc_pix_list,cell_pix_list)))
                                rat_member(count)=member_pix/size(nuc_pix_list,1);
                            end
                            
                            [~,in_cell_nuc_pix]=max(rat_member);
                            nuc_pix_list=prop_nuc_bw(in_cell_nuc_pix).PixelIdxList;
                            
                            %get data
                            cell_area=max([prop_cell_bw.Area]);
                            nuc_area=max([prop_nuc_bw.Area]);
                            
                            cell_int=sum(crop_orig_cell(cell_pix_list));
                            nuc_int=sum(crop_orig_cell(nuc_pix_list));
                            nuc_nuc_int=sum(crop_orig_nuc(nuc_pix_list));
                            nuc_height=prop_nuc_bw(in_cell_nuc_pix(1)).BoundingBox(1,end);
                            nuc_proj_area=regionprops(sum(prop_nuc_bw(in_cell_nuc_pix(1)).FilledImage,3)>0,'Area');
                            cell_proj_area=regionprops(sum(prop_cell_bw.FilledImage,3)>0,'Area');
                            %% cell and nucleus projected areas
%                             pro_cell_proj=sum(cell_bw,3);
%                             pro_nuc_proj=sum(nuc_bw,3);                            
%                             pro_cell_proj_prop=regionprops(pro_cell_proj);
%                             pro_nuc_proj_prop=regionprops(pro_nuc_proj);
                            
%                             imshow(pro_nuc_proj)
%                             hold on
%                             scatter(pro_nuc_proj_prop(end).BoundingBox(1),pro_nuc_proj_prop(end).BoundingBox(2))
%                             
                            
                            %%
                            
                            %             phos_cell_int=sum(crop_orig_phos(cell_pix_list));
                            %             phos_nuc_int=sum(crop_orig_phos(nuc_pix_list));
                            
                            %%
                            
                            %% gather data
                            [pro_xls_label{1,1:18,obj_count}]=deal(...
                                [['Series_',num2str(series_count)],['Cell', num2str(obj_count)],['_', num2str(file_count),],file_list(file_count).name(1:end-4)],...
                                t_count, ... %1-3
                                nuc_int, cell_int, nuc_int/cell_int, nuc_area,cell_area,...%4-8
                                nuc_nuc_int,nuc_int/nuc_nuc_int,nuc_int/nuc_area, cell_int/cell_area,...%9-12
                                (nuc_int/nuc_area)/(cell_int/cell_area),cell_thresh_fac,nuc_thresh_fac,nuc_int/(cell_int-nuc_int),...%13-15
                                nuc_height,cell_proj_area.Area,nuc_proj_area.Area...%16-18
                                );
                            
                            
                            %% make image to write
                            
                            out_img=[mat2gray(imadjust(reshape(crop_orig_cell,[],size(cell_bw,2)*size(cell_bw,3))));
                                reshape(cell_bw,[],size(cell_bw,2)*size(cell_bw,3));
                                reshape(nuc_bw,[],size(cell_bw,2)*size(cell_bw,3));
                                mat2gray(imadjust(reshape(crop_orig_nuc,[],size(cell_bw,2)*size(cell_bw,3))))%;
                                %reshape(phos_stack(cell_pix_list),[],size(cell_bw,2)*size(cell_bw,3));
                                ];
                            max_img=max(orig_cell_stack,[],3);
                            proj_cell_nuc=[sum(cell_bw,3) sum(nuc_bw,3)]>1;
                            %sum_img=sum(orig_cell_stack,3,'double');
                            %         imtool(out_img,[])
                            %% write image
                            if t_count==1
                                imwrite(out_img,[img_write_path, 'Series',num2str(series_count), '_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff');%,'Compression', 'jpeg', 'RowsPerStrip', 16);
                                
                                imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff');
                                
                                % imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                % 'tiff');
                                
                                imwrite(proj_cell_nuc,[img_write_path, 'Proj_BW',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff');
                            elseif t_count>1
                                imwrite(out_img,[img_write_path, 'Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff','WriteMode','append');%, 'Compression', 'jpeg', 'RowsPerStrip', 16);
                                
                                imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff','WriteMode','append');
                                
                                %imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                % 'tiff','WriteMode','append');
                                
                                imwrite(proj_cell_nuc,[img_write_path, 'Proj_BW',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count),file_list(file_count).name(1:end-4), '.tiff'],...
                                    'tiff','WriteMode','append');
                            end
                            clear out_img
                            
                            disp(['Series ' num2str(file_count) ' of ' num2str(size(file_list,1)) ' time ' num2str(t_count) ' of ' num2str(T_num)])
                        end % for 3d area thresh
                        %% clear var
                        clear rat_member crop_smo_cell crop_cell crop_nuc crop_orig_cell crop_orig_nuc crop_phos crop_orig_phos
                    end % for object
%                     xls_label=vertcat(xls_label,reshape(permute(pro_xls_label, [1 3 2]),[],18));
                    xls_label(1,1:18)=pro_xls_label(:);
                    
                    clear data
                    clear cell_stack phos_stack orig_phos_stack nuc_stack orig_cell_stack cell_stack nuc_stack smooth_cell smooth_nuc
                    clear sum_cell sum_smooth_cell re_smooth_cell nuc_bw orig_nuc_stack pix_list pro_xls_label
                    clear cumu_nuc_int
                end % if bw_prop is empty
                
            end %for T count
            
        end %for series count
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

