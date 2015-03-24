%% Pax_analyse6_nd2 
% based on Pax_analyse_5 has been modularised and background subtraction
% has been added. 

clear
clc
%% THINGS TO CHANGE IN THIS SECTION
% >>>>>>>CHANGE FILE PATH BELOW!<<<<<<<<
%root_root='/mnt/mbi/images/naotaka/'; %this is base path in laksa/satay
root_root='C:\Users\Nyx\Desktop\test\';

%>>>>>>>>IF MIXED DATA SET SERIES BELOW. IF NOT MIXED SET TO 0
series_to_read=0;%0; %[23,24,25,26,27]%use this to specify series. when not in use set to zero

%>>>>>> CHANGE EXTENSION TO .nd2 BELOW<<<<<<<<<<
file_ext='*.nd2'; %e.g. '*.mvd2' '*.nd2'

cell_chanl=2;
nuc_chanl=1;
% phos_chanl=3;
cell_2d_area_thresh=100000;%nikon A1R 1024x1024
% cell_2d_area_thresh=10000; %spinning disk normal cells
%  cell_2d_area_thresh=1000; %spinning disk fak-/- or vin-/- (small cells)
cell_3d_area_thresh= 300000; %Nikon A1R for 1024x1024 cells
%cell_3d_area_thresh= 30000;%spinning disk normal cells
%  cell_3d_area_thresh= 3000;%spinning disk fak-/- or vin-/- (small cells)
nuc_3d_area_thresh=7000;
%  nuc_3d_area_thresh=700;


%% DO NOT MODIFY PROGRAM AFTER THIS POINT
add_path1=fullfile(fileparts(mfilename('fullpath')), 'image_tools', 'bfmatlab', filesep);
add_path2=fullfile(fileparts(mfilename('fullpath')), 'dlmcell', filesep);
addpath(add_path1)
addpath(add_path2)

%% read image

%root_root='C:\Users\Aneesh\Documents\Data\spinning disk\';

fold_list=dir(root_root); %% check contents of folder
for fold_count=3:size(fold_list,1)
    %% set root specify output path and make output directory 
    [root,img_write_path,file_list] = an_set_paths(root_root, fold_list(fold_count).name,file_ext);
    
    
    
    for file_count=1:size(file_list,1)
        %file_path=[root, file_list.name]; %use this normally
        file_path=[root, file_list(file_count).name]; %modification
        %% initialize file reader and get number of series
        reader = an_gen_reader( file_path );
        %% initialize storage cell
        [xls_label{1,1:14}]=deal(...
            'Cell','Time', ... %1-2
            'Nuc Intensity', 'Cell Intensity','Ratio', 'Nuc Area','Cell Area',...%3-7
            'Nuc Label Int ','nuc_coloc_ratio','Avg Nuc Int', 'Avg Cell Int',...%8-11
            'Avg Int Ratio','Cell thresh val','Nuc thresh val'...%12-14
            );
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
                
                cell_stack=mat2gray(orig_cell_stack);
                nuc_stack=mat2gray(orig_nuc_stack);
                %             smooth_cell=smooth3(cell_stack);
                %             smooth_nuc=smooth3(nuc_stack);
                smooth_cell=smooth3(smooth3(cell_stack,'gaussian'),'gaussian');
                smooth_nuc=smooth3(smooth3(nuc_stack,'gaussian'),'gaussian');
                
                sum_cell=sum(cell_stack,3);
                sum_smooth_cell=sum(smooth_cell,3);
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
                            
                            tmp=cell_stack(:,:,count);
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
                        cell_thresh_fac=graythresh(crop_cell);
                        nuc_thresh_fac=graythresh(crop_nuc);
                        
                        %apply thresh
                        %             unre_crop_smo_cell=reshape(re_crop_smo_cell,[],size(crop_smo_cell,2),size(crop_smo_cell,3));
                        %             cell_bw=unre_crop_smo_cell>cell_thresh_fac;
                        cell_bw=crop_cell>cell_thresh_fac;
                        nuc_bw=crop_nuc>nuc_thresh_fac;
                        %fill holes
                        cell_bw=bwareaopen(imfill(imdilate(cell_bw,se2),'holes'),cell_3d_area_thresh);
                        nuc_bw=bwareaopen(imfill(nuc_bw,'holes'),nuc_3d_area_thresh);
                        % make new cell_bw which includes both cell and nuclear
                        cell_bw=cell_bw|nuc_bw;
                        %get pix list
                        prop_cell_bw=regionprops(cell_bw,'Area', 'PixelIdxList');
                        prop_nuc_bw=regionprops(nuc_bw,'Area', 'PixelIdxList');
                        
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
                            %             phos_cell_int=sum(crop_orig_phos(cell_pix_list));
                            %             phos_nuc_int=sum(crop_orig_phos(nuc_pix_list));
                            
                            %%
                            
                            %% gather data
                            [pro_xls_label{1,1:14,obj_count}]=deal(...
                                [['Series_',num2str(series_count)],['Cell', num2str(obj_count)],['_', num2str(file_count),]],t_count, ... %1-3
                                nuc_int, cell_int, nuc_int/cell_int, nuc_area,cell_area,...%4-8
                                nuc_nuc_int,nuc_int/nuc_nuc_int,nuc_int/nuc_area, cell_int/cell_area,...%9-12
                                (nuc_int/nuc_area)/(cell_int/cell_area),cell_thresh_fac,nuc_thresh_fac...%13-15
                                );
                            
                            
                            %% make image to write
                            
                            out_img=[reshape(crop_cell,[],size(cell_bw,2)*size(cell_bw,3));
                                reshape(cell_bw,[],size(cell_bw,2)*size(cell_bw,3));
                                reshape(nuc_bw,[],size(cell_bw,2)*size(cell_bw,3));
                                reshape(crop_nuc,[],size(cell_bw,2)*size(cell_bw,3))%;
                                %reshape(phos_stack(cell_pix_list),[],size(cell_bw,2)*size(cell_bw,3));
                                ];
                            max_img=max(orig_cell_stack,[],3);
                            %sum_img=sum(orig_cell_stack,3,'double');
                            %         imtool(out_img,[])
                            %% write image
                            if t_count==1
                                imwrite(out_img,[img_write_path, 'Series',num2str(series_count), '_', num2str(file_count), '_cell_', num2str(obj_count), '.tiff'],...
                                    'tiff');%,'Compression', 'jpeg', 'RowsPerStrip', 16);
                                
                                imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count), '.tiff'],...
                                    'tiff');
                                % imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
                                % 'tiff');
                                
                            elseif t_count>1
                                imwrite(out_img,[img_write_path, 'Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count), '.tiff'],...
                                    'tiff','WriteMode','append');%, 'Compression', 'jpeg', 'RowsPerStrip', 16);
                                
                                imwrite(max_img,[img_write_path, 'Max_Proj_Series',num2str(series_count),'_', num2str(file_count), '_cell_', num2str(obj_count), '.tiff'],...
                                    'tiff','WriteMode','append');
                                
                                %imwrite(uint16(sum_img),[img_write_path, 'Sum_Proj_Series',num2str(series_count), '_cell_', num2str(obj_count), '.tiff'],...
                                % 'tiff','WriteMode','append');
                                
                            end
                            clear out_img
                            disp(['Series ' num2str(series_count) ' of ' num2str(series_num) ' time ' num2str(t_count) ' of ' num2str(T_num)])
                        end % for 3d area thresh
                        %% clear var
                        clear rat_member crop_smo_cell crop_cell crop_nuc crop_orig_cell crop_orig_nuc crop_phos crop_orig_phos
                    end % for object
                    xls_label=vertcat(xls_label,reshape(permute(pro_xls_label, [1 3 2]),[],14));
                    clear data
                    clear cell_stack phos_stack orig_phos_stack nuc_stack orig_cell_stack cell_stack nuc_stack smooth_cell smooth_nuc
                    clear sum_cell sum_smooth_cell re_smooth_cell nuc_bw orig_nuc_stack pix_list pro_xls_label
                    clear cumu_nuc_int
                end % if bw_prop is empty
                
            end %for T count
            
        end %for series count
    end %for file list count
    %% write data
    
    xls_file_name=[img_write_path, fold_list(fold_count).name '.txt'];
    
    dlmcell(xls_file_name,xls_label,',')
    
    %write excel file if OS is windows
    if ispc
        xls_file_name2=[img_write_path, fold_list(fold_count).name '.xls'];
        xlswrite(xls_file_name2,xls_label,'Sheet1' );
    end
    
    
    disp('done writing excel file')
    clear xls_label final_data_out out_* xls_file_*
    clear nuc_avg_out nuc_all_out nuc_cumu_x out_cumu_nuc_int
    clear cell_avg_out cell_all_out cell_cumu_x out_cumu_cell_int
    
end %for fold count
