function [ root,img_write_path,file_list ] = an_set_paths( root_root,fold_name,file_ext )
%AN_SET_PATHS Summary of this function goes here
%   Detailed explanation goes here
root=fullfile(root_root,fold_name,filesep); %change path as required
    img_write_path=fullfile(root, ['results_', date],filesep);
    mkdir(img_write_path);

    file_list=dir(fullfile(root, file_ext));
    
end

