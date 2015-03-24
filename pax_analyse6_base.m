%% Pax_analyse6
% based on Pax_analyse_5 has been modularised and background subtraction
% has been added.

clear
clc
tic
%% THINGS TO CHANGE IN THIS SECTION
% >>>>>>>CHANGE FILE PATH BELOW!<<<<<<<<
%root_root='/mnt/mbi/images/naotaka/'; %this is base path in laksa/satay
root_root='D:\Aneesh\test\';

%>>>>>>>>IF MIXED DATA SET SERIES BELOW. IF NOT MIXED SET TO 0
series_to_read=0; %[23,24,25,26,27]%use this to specify series. when not in use set to zero

%>>>>>> CHANGE EXTENSION TO .nd2 BELOW<<<<<<<<<<
file_ext='mvd2'; %e.g. '*.mvd2' '*.nd2'

cell_chanl=3;
nuc_chanl=1;
meas_chan=2;
% phos_chanl=3;
cell_2d_area_thresh=1000; %spinning disk normal cells
%  cell_2d_area_thresh=1000; %spinning disk fak-/- or vin-/- (small cells)
cell_3d_area_thresh= 50000;%spinning disk normal cells
%  cell_3d_area_thresh= 3000;%spinning disk fak-/- or vin-/- (small cells)
nuc_3d_area_thresh=70000;
%  nuc_3d_area_thresh=700;

%% DO NOT MODIFY PROGRAM AFTER THIS POINT
% mfilename()
% an_add_req_paths(fileparts(mfilename('fullpath')))

add_path1=fullfile(fileparts(mfilename('fullpath')), 'image_tools', 'bfmatlab', filesep);
add_path2=fullfile(fileparts(mfilename('fullpath')), 'dlmcell', filesep);
add_path3=fullfile(fileparts(mfilename('fullpath')), 'modules', filesep);
addpath(add_path1)
addpath(add_path2)
addpath(add_path3)

%% read image

[Outfiles, outDir]= an_recdir( root_root,file_ext);
%root_root='C:\Users\Aneesh\Documents\Data\spinning disk\';

% fold_list=dir(root_root); %% check contents of folder
for fold_count=1:numel(Outfiles)
    root=outDir{fold_count};
    file_path=Outfiles{fold_count};
    img_write_path=fullfile(root, ['results_', date],filesep);
    mkdir(img_write_path);
%     [root,img_write_path,file_list] = an_set_paths(root_root, fold_list(fold_count).name,file_ext);
    
    %file_path=[root, file_list.name]; %use this normally
%     file_path=[root, file_list(end).name]; %modification
    %% initialize file reader and get number of series
    [reader,series_num,pro_series_count ]= an_gen_reader(file_path,series_to_read);
    
    %    reader=bfGetReader(file_path);
    %         if series_to_read==0
    %             series_num=reader.getSeriesCount;
    %             pro_series_count=1:series_num;
    %         else
    %             pro_series_count=series_to_read;
    %         end
    
    
    
    %% initialize storage cell
    [xls_label{1,1:16}]=deal(...
        'Cell','Time', ... %1-2
        'Nuc Intensity', 'Cell Intensity','Ratio', 'Nuc Area','Cell Area',...%3-7
        'Nuc Label Int ','nuc_coloc_ratio','Avg Nuc Int', 'Avg Cell Int',...%8-11
        'Avg Int Ratio','Cell thresh val','Nuc thresh val','Std in Nuc', 'ProjCellArea'...%12-15
        );
    %%
    for series_count=1:pro_series_count
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
            
            z_range=1:Z_num;
            
                        orig_cell_stack=bfGetStack(reader,meas_chan,z_range,t_count-1);
                        orig_nuc_stack=bfGetStack(reader,nuc_chanl,z_range,t_count-1);
            
                        temp=bfGetStack(reader,cell_chanl,z_range,t_count-1);
%             orig_cell_stack=bfGetStack(reader,meas_chan,z_range,t_count-1);
%             orig_nuc_stack=bfGetStack(reader,nuc_chanl,z_range,t_count-1);
%             
%             temp=bfGetStack(reader,cell_chanl,z_range,t_count-1);
            %             temp2=reshape(temp,[],size(temp,1)*size(temp,2));
            %             temp3=reshape(adapthisteq(temp2),size(temp));
            %             cell_stack=mat2gray(temp3);
            cell_stack=mat2gray(temp);
            nuc_stack=mat2gray(orig_nuc_stack);
            pro_cell_stack=imfill(smooth3(mat2gray(cell_stack+nuc_stack)));
            %             cell_stack=reshape(imfill(reshape(pro_cell_stack,[],size(pro_cell_stack,1)*size(pro_cell_stack,2))),size(pro_cell_stack));
            %             smooth_cell=smooth3(cell_stack);
            %             smooth_nuc=smooth3(nuc_stack);
            
            [bw_re_prop,bw_sum_cell ] = an_smooth_imgprops_paxanalyse(cell_stack);
            
            
            %             BS_orig_cell_stack= orig_cell_stack;
            %             BS_orig_nuc_stack = nuc_stack;
            
            BS_orig_cell_stack= an_noise_rem_paxanalyse(bw_sum_cell,orig_cell_stack);
            BS_orig_nuc_stack = an_noise_rem_paxanalyse(bw_sum_cell,nuc_stack);
            
            [pro_xls_label{1,1:16,1}]=deal(zeros);
            if ~isempty(bw_re_prop)
                %                 pix_list=bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).PixelIdxList;
                %                 crop_bw=imcrop(bw_sum_cell,bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).BoundingBox);
                %                 bb=bw_re_prop([bw_re_prop.Area]==max([bw_re_prop.Area])).BoundingBox;
                
                obj_idx=find([bw_re_prop.Area]>cell_2d_area_thresh);
                
                for obj_count=1:size(obj_idx,2)
                    
                    bb=bw_re_prop(obj_idx(obj_count)).BoundingBox;
                    pix_list=bw_re_prop(obj_idx(obj_count)).PixelIdxList;
                    %                     [tmp, cmc, cc, cn, cp, cop]=deal(zeros(size(cell_stack,1),size(cell_stack,2)));
                    %                     [coc, noc]=deal(zeros(size(cell_stack,1),size(cell_stack,2),'uint16'));
                    
                    %crop images
                    %                     crop_smo_cell=an_paxanalyse_crop( smooth_cell,pix_list,bb );
                    crop_cell=an_paxanalyse_crop( cell_stack,pix_list,bb );
                    crop_nuc=an_paxanalyse_crop( nuc_stack,pix_list,bb );
                    crop_orig_cell=an_paxanalyse_crop(BS_orig_cell_stack ,pix_list,bb );
                    crop_orig_nuc=an_paxanalyse_crop( BS_orig_nuc_stack,pix_list,bb );
                    
                    
                    
                    %% get pix list
                    %                     se2=strel('disk',15);
                    
                    %                     temp2=reshape(crop_cell,[],size(crop_cell,1)*size(crop_cell,2));
                    %                     temp3=reshape(adapthisteq(temp2),size(crop_cell));
                    %                     [cell_bw,cell_thresh_fac]=an_gray(smooth3(temp3));
                    
                    [cell_bw,cell_thresh_fac]=an_gray(imfill(smooth3(crop_cell)));
                    [nuc_bw,nuc_thresh_fac]=an_gray(smooth3(crop_nuc));
                    cell_bw=cell_bw|nuc_bw;
                    %fill holes
                    se2=strel('disk',5);
                    cell_bw=bwareaopen(imfill(imdilate(cell_bw,se2),6,'holes'),cell_3d_area_thresh);
                    nuc_bw=bwareaopen(imfill(nuc_bw,'holes'),nuc_3d_area_thresh);
                    
                    
                    %se2=strel('ball',15,5);
                    %                     cell_bw=bwareaopen(imfill(imdilate(cell_bw,se2),'holes'),cell_3d_area_thresh);
                    %                     nuc_bw=bwareaopen(imfill(nuc_bw,'holes'),nuc_3d_area_thresh);
                    
                    % make new cell_bw which includes both cell and nuclear
                    
                    
                    %                     imshow([an_reshape_img_paxanalyse(cell_bw);an_reshape_img_paxanalyse(nuc_bw)])
                    %%
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
                        nuc_std=std(double(crop_orig_cell(nuc_pix_list)));
                        %             phos_cell_int=sum(crop_orig_phos(cell_pix_list));
                        %             phos_nuc_int=sum(crop_orig_phos(nuc_pix_list));
                        
                        %%
                        
                        %% gather data
                        [pro_xls_label{1,1:16,obj_count}]=deal(...
                            [['Series_',num2str(series_count)],['Cell', num2str(obj_count)]],t_count, ... %1-3
                            nuc_int, cell_int, nuc_int/cell_int, nuc_area,cell_area,...%4-8
                            nuc_nuc_int,nuc_int/nuc_nuc_int,nuc_int/nuc_area, cell_int/cell_area,...%9-12
                            (nuc_int/nuc_area)/(cell_int/cell_area),cell_thresh_fac,nuc_thresh_fac,...%13-15
                            nuc_std,bw_re_prop(obj_idx(obj_count)).Area...%16
                            );
                        
                        
                        %% make image to write
                        
                        out_img=[an_reshape_img_paxanalyse(mat2gray(crop_orig_cell));
                            an_reshape_img_paxanalyse(crop_cell);
                            an_reshape_img_paxanalyse(cell_bw);
                            an_reshape_img_paxanalyse(nuc_bw);
                            an_reshape_img_paxanalyse(crop_nuc)%;
                            %reshape(phos_stack(cell_pix_list),[],size(cell_bw,2)*size(cell_bw,3));
                            ];
                        max_img=max(BS_orig_cell_stack,[],3);
                        %sum_img=sum(orig_cell_stack,3,'double');
                        %         imtool(out_img,[])
                        %% write image
                        an_write_img_paxanalyse(out_img,max_img,img_write_path,series_count,obj_count,t_count)
                        
                        clear out_img
                        disp(['Series ' num2str(series_count) ' of ' num2str(series_num) ' time ' num2str(t_count) ' of ' num2str(T_num)])
                    end % for 3d area thresh
                    %% clear var
                    clear rat_member crop_smo_cell crop_cell crop_nuc crop_orig_cell crop_orig_nuc crop_phos crop_orig_phos
                end % for object
                xls_label=vertcat(xls_label,reshape(permute(pro_xls_label, [1 3 2]),[],16));
                clear data
                clear cell_stack phos_stack orig_phos_stack nuc_stack orig_cell_stack cell_stack nuc_stack smooth_cell smooth_nuc
                clear sum_cell sum_smooth_cell re_smooth_cell nuc_bw orig_nuc_stack pix_list pro_xls_label
                clear cumu_nuc_int
            end % if bw_prop is empty
            
        end %for T count
        
    end %for series count
    
    %% write data
    [~,b,~]=fileparts(file_path);
    xls_file_name=[img_write_path, b, '.txt'];
    
    dlmcell(xls_file_name,xls_label,',')
    
    %write excel file if OS is windows
    if ispc
        xls_file_name2=[img_write_path, b '.xls'];
        xlswrite(xls_file_name2,xls_label,'Sheet1' );
    end
    
    
    disp('done writing excel file')
    clear xls_label final_data_out out_* xls_file_* b
    clear nuc_avg_out nuc_all_out nuc_cumu_x out_cumu_nuc_int
    clear cell_avg_out cell_all_out cell_cumu_x out_cumu_cell_int
    
end %for fold count
toc