function img_stack=bfGetStack(reader,chanl,z_range,T_count)
%% to return desired z-stack at desired time point in a 5-D file

%get image parameters
dim_order=reader.getDimensionOrder;
C_num=reader.getSizeC;
Z_num=reader.getSizeZ;
T_num=reader.getSizeT;

%use appropriate order depending on dimension order
if strcmp(dim_order,'XYCZT')
    %determine relevant planes
    planes=((chanl:C_num:C_num*Z_num)+C_num*Z_num*T_count);
elseif strcmp(dim_order,'XYZCT')
    %planes=0;
end

%read planes
counter=0;
for z_count=planes(z_range)
    counter=counter+1;
    img_stack(:,:,counter)=bfGetPlane(reader,z_count);
end

end