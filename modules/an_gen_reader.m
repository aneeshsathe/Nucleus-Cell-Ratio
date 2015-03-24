function [reader,series_num,pro_series_count ]= an_gen_reader( file_path,series_to_read )
%AN_GEN_READER Takes file path and generates reader for use by bfopen
%   Detailed explanation goes here

reader=bfGetReader(file_path);
        if series_to_read==0
            series_num=reader.getSeriesCount;
            pro_series_count=1:series_num;
        else
            pro_series_count=series_to_read;
        end


end

