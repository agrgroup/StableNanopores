% clc;
% clear all
function [count]=search(llsf,ulsf,llma,ulma,llmi,ulmi,n)
%Input
% n=8;
%llsf=lower limit of shape factor range
%ulsf=upper limit of shape factor range
%llma=lower limit of major axis range
%ulma=upper limit of major axis range
%llmi=lower limit of minor axis range
%ulmi=upper limit of minor axis range
% n=6;
% llsf=0.7;%lower limit of shape factor range
% ulsf=0.9;%upper limit of shape factor range
% llma=7;%lower limit of major axis range
% ulma=9.5;%upper limit of major axis range
% llmi=6;%lower limit of minor axis range
% ulmi=9;%upper limit of minor axis range

str=string(n);
load("stable_nanopores"+str+".mat");
%Extraction of data from sheets
filename1="shape_factor"+str+".xlsx";
datasf=xlsread(filename1);
datasf2=datasf;
sfIndexes = (datasf2(:,3) >= llsf) & (datasf2(:,3) <=ulsf);
 if nnz(sfIndexes)>0
  datasf2=datasf2(sfIndexes,:);
  maIndexes = (datasf2(:,4) >=llma) & (datasf2(:,4) <=ulma);
  if nnz(maIndexes)>0
       datasf2=datasf2(maIndexes,:);
    Indexes = (datasf2(:,5) >=llmi) & (datasf2(:,5) <=ulmi);
    if nnz(Indexes)>0
       datasf2=datasf2(Indexes,:);
       count=nnz(Indexes)
       for i=1:count
           A=polys_ind{datasf2(Indexes(i),1)};
           xyz_A(n,A);
       end 
    else
        count=0;
        return;
    end
  else
    count=0;
    return;
  end
else
   count=0; 
   return;
end

end

