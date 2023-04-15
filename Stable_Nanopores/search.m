function [count]=search(llsf,ulsf,llma,ulma,llmi,ulmi,n)
%Input:
%n=Number of triangles in the polyiamond (Number of atoms missing in the nanopore)
%llsf=lower limit of shape factor range
%ulsf=upper limit of shape factor range
%llma=lower limit of major axis range
%ulma=upper limit of major axis range
%llmi=lower limit of minor axis range
%ulmi=upper limit of minor axis range

%Output:
%count:Returns the count of nanopores that has the specified range of properties
%given as input
%An xyz file named "inrange_search"+"n".xyz (n=Number of triangles in the 
%polyiamond or number of atoms missing in the nanopore) containing the coordinates of
%nanopores that fall in the input range is also generated. 

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

