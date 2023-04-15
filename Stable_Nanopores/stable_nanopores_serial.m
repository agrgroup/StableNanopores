clear all
clc
tic
%https://oeis.org/A001420
fixed_oeis=[2, 3, 6, 14, 36, 94, 250, 675, 1838, 5053, 14016, 39169, 110194, 311751, 886160, 2529260, 7244862, 20818498, 59994514, 173338962];
%https://oeis.org/A070765
free_oeis=[1, 1, 1, 3, 4, 12, 24, 66, 159, 444, 1161, 3226, 8785, 24453, 67716, 189309, 528922, 1484738, 4172185, 11756354, 33174451, 93795220, 265565628];

%Inputs:
%n=Number of triangles in the polyiamond (Number of atoms missing in the nanopore)
n=14

remove_dangling=1;%Set remove_dangling=1 or 0, remove_dangling=1 means that dangling moeities and dangling bonds will be omitted and stable nanores will be saved, remove_dangling=0 means that dangling moeities and dangling bonds will not omitted and free polyiamonds without holes will be saved. This can be used to verify the code with the OEIS numbers, 
print=1;%print=1 or 0, if print=1, saves the coordinates in xyz file and viceversa

%Outputs:
%If remove_dangling=1,
%Saves .mat file named stable_nanoporesx.mat (x=n) which contains polys_ind
%and ind. Polys_ind is a cell array that contains all representative
%polyiamond coordinates. ind is the indices number with respect to the
%initially generated set of fixed polyiamonds with inverted triangle at
%origin.
%If remove_dangling=0;
%Saves .mat file named free_polyiamondsx.mat (x=n) which contains polys_ind
%and ind. Polys_ind is a cell array that contains all representative
%polyiamond coordinates. ind is the indices number with respect to the
%initially generated set of fixed polyiamonds with inverted triangle at
%origin. All free polyiamonds without holes will be saved here.

if n==1 || n==2
    xyz(n,[]);
    ind_stable=1;
    if n==1
    polys_ind={[0,0]};
    str=string(n);
    save("stable_nanopores"+str+".mat",'ind_stable','polys_ind')
    else 
     polys_ind={[0,0],[0,1]};
    str=string(n);
    save("stable_nanopores"+str+".mat",'ind_stable','polys_ind')
    end
else
    
    [numpolys1,polys1] = fixed_inverted(n); %generates set of fixed polyiamonds with inverted triangle at origin

%Uncomment the followig lines to check the accuracy with respect to OEIS numbers of fixed polyiamonds

%[numpolysu,polysu] = fixed_upright(n);%generates set of fixed polyiamonds with upright triangle at origin
%fixed_polyiamonds=numpolysu+numpolys1
%oeis_fixed_polyiamonds=fixed_oeis(n)

     fixed_polyiamonds_with_inverted_triangle_at_origin=numpolys1
     
    data=zeros(numpolys1,3);
    %data_uid1={};
     for i=1:numpolys1 
       A=polys1{i};%A=representative coordinates of the polyiamond
       [ho,vrv,uid,sett]=holes_db_elimination(A,remove_dangling);%Checks for nanopores with dangling bonds and non-bonded atoms or polyiamond with holes and vertices shared by 5 triangles
%ho indicates to preence of holes and vertices shared by 5 triangles. vrv is the vertex repetition vector. uid is the unique id of polyiamonds.     
data_uid{i}=[];       
if ho==1 %Saves if no holes and vertices shared by 5 triangles are present          
         %data(i,:)= [vrv uid i];
         data(i,:)= [vrv uid i];
         data_uid{i}=sett;
       end
     end

%finds and deletes rows containg zeroes
In=find(data(:,3));
data=data(In,:);
data_uid=data_uid(In);
%Use sortrows to organize data according to vertex repetition vector
[data,In]=sortrows(data);
data_uid=data_uid(In);
fdata=flip(data(:,1));
[C,ia,ic] = unique(fdata);
si=[0;-(ia-length(data)-1)];%Extracts where a particular vertex repetition vector index ends
ind1=[];
ind2=[];

%Eliminating symmetrical isomers by comparison amongst polyiamonds having same vertex repetition vector
for i=2:size(si,1)
 strong=data_uid(:,si(i-1)+1:si(i));
       [indi1,indi2]=free_symmetry_operations(strong,data(si(i-1)+1:si(i),3));
      ind1=[ind1,indi1'];  
      ind2=[ind2,indi2'];  
end
if remove_dangling==0%Saves free polyiamonds
    for w=1:size(ind2,2)
       polys_ind{w}=polys1{ind2(w)};
    end
    free_polyiamonds=numel(ind2)
    OEIS_number=free_oeis(n)
    save("free_polyiamonds"+str+".mat",'ind2','polys_ind')
toc
else     
%Finding properties, eliminating nanopores with dangling moeties and saving
%data and xyz files
    ind_stable=[];
    order=[];
    minors=[];
    majors=[];
    shape_g=[];
    
    for w=1:size(ind2,2)
       A=polys1{ind2(w)};%A=representative coordinates of the polyiamond
%The function calculates and returns the values of properties (major axis, minor axis and shapefactor), checks for presence of dangling moeities and prints the coordinates in xyz file if the calculated values fall in the range given as input).
       [sf,ma,mi,moeityr]= properties_moeity_xyz(A,n);
       if moeityr==0
%Saving properies
          ind_stable=[ind_stable;ind2(w)];
          majors=[majors; ma];
          minors=[minors; mi];
          shape_g=[shape_g;sf];
       end
    end
    for i=1:numel(ind_stable)
       polys_ind{i}=polys1{ind_stable(i)};%Saving Reperesntative coordinates of stable nanopores 
    end   
% Sorting based on Shape factor
    [shape_sorted,inde]=sort(shape_g);
%Sorting based on Minor axis
    [minors_sorted,indemi]=sort(minors);
%Sorting based on Major axis
    [majors_sorted,indema]=sort(majors);
    format long g
    tab=[inde shape_sorted];
    str=string(n);

%Saving properties after sorting them based on their sorted shape factor values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted shape factor values
%4th and 5th column denotes to their corresponding major axis and minor axis values respectively
    writematrix([inde ind_stable(inde) shape_sorted majors(inde) minors(inde)],"shape_factor"+str+".xlsx")
  
 %Saving properties after sorting them based on their minor axis values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted minor axis values
%4th and 5th column denotes to their corresponding major axis values and shape factor values respectively
    writematrix([indemi ind_stable(indemi) minors_sorted majors(indemi) shape_g(indemi)],"minor_axis"+str+".xlsx")
  
%Saving properties after sorting them based on their major axis values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted major axis values
%4th and 5th column denotes to their corresponding minor axis values and shape factor values respectively
    writematrix([indema ind_stable(indema) majors_sorted minors(indema) shape_g(indema)],"major_axis"+str+".xlsx")

 %Saving of Data
    str=string(n);
    number_of_stable_nanopores=numel(ind_stable)
    save("stable_nanopores"+str+".mat",'ind_stable','polys_ind')
end
end
% if print==1
% xyz_from_data(n);
% end
% %Searching for nanopores with the below range of properties. Comment the
% %following to skip the searching process.
% llsf=0;%lower limit of shape factor range
% ulsf=1;%upper limit of shape factor range
% llma=0;%lower limit of major axis range
% ulma=inf;%upper limit of major axis range
% llmi=0;%lower limit of minor axis range
% ulmi=inf;%upper limit of minor axis range
% count=search(llsf,ulsf,llma,ulma,llmi,ulmi,n);%Number of nanopores that fall in the above given range
toc
