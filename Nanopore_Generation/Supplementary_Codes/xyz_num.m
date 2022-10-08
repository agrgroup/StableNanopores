function [xyz_in] = xyz_num(polys,n,num,type)
%Extracts data from sheets 
%input type 1 for shape factor, 2 for minor axis and 3 for major axis
str=string(n);
if type==1
filename3="shape_factor"+str+".xlsx";
elseif type==2 
 filename3="minor_axis"+str+".xlsx";   
else
    filename3="major_axis"+str+".xlsx";
end
data=xlsread(filename3);
polys_index=data(:,1);

if num>1
num=polys_index(end);
else
  num=polys_index(1);
end

A=polys{num};
[xyz_in]=xyz_A(n,num,type,A);
end