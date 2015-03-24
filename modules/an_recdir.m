function [Outfiles, outDir]= an_recdir( baseDir,searchExpression )
% function [folders, path]= an_recdir( root,wldcard )
%AN_RECDIR Summary of this function goes here
%   Detailed explanation goes here
%%
dstr = dir(baseDir);%search current directory and put results in structure
Outfiles = {};
outDir = {};
for II = 1:length(dstr)
    
    if ~dstr(II).isdir && ~strcmp(dstr(II).name(1),'.') && ~isempty(regexp(dstr(II).name(end-4:end),searchExpression,'match'))
        %         dstr(II).name(end-4:end)
        %look for a match that isn't a directory
        Outfiles{length(Outfiles)+1} =fullfile(baseDir, dstr(II).name);
        outDir{length(outDir)+1} = fullfile(baseDir,filesep);
    elseif dstr(II).isdir && ~strcmp(dstr(II).name,'.') && ~strcmp(dstr(II).name,'..')
        %if it is a directory(and not current or up a level), search in that
        pname = fullfile(baseDir,dstr(II).name);
        %         OutfilesTemp=an_recdir(pname,searchExpression);
        [OutfilesTemp, outDirTemp]=an_recdir(pname,searchExpression);
        if ~isempty(OutfilesTemp)
            %if recursive search is fruitful, add it to the current list
            Outfiles((length(Outfiles)+1):(length(Outfiles)+length(OutfilesTemp))) = OutfilesTemp;
            outDir((length(outDir)+1):(length(outDir)+length(outDirTemp))) = outDirTemp;
        end
    end
end


end

