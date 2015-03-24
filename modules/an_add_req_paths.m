function [ ] = an_add_req_paths(in_name)
%ADD_REQ_PATHS Summary of this function goes here
%   Detailed explanation goes here
% add_path1=fullfile(fileparts(mfilename('fullpath')), 'image_tools', 'bfmatlab', filesep)
% add_path2=fullfile(fileparts(mfilename('fullpath')), 'dlmcell', filesep)
% addpath(add_path1)
% addpath(add_path2)

add_path1=fullfile(in_name, 'image_tools', 'bfmatlab', filesep);
add_path2=fullfile(in_name, 'dlmcell', filesep);
addpath(add_path1)
addpath(add_path2)
end

