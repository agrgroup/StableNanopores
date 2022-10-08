clc
tic
%https://oeis.org/A001420
fixed_oeis=[2, 3, 6, 14, 36, 94, 250, 675, 1838, 5053, 14016, 39169, 110194, 311751, 886160, 2529260, 7244862, 20818498, 59994514, 173338962];
%https://oeis.org/A070765
free_oeis=[1, 1, 1, 3, 4, 12, 24, 66, 159, 444, 1161, 3226, 8785, 24453, 67716, 189309, 528922, 1484738, 4172185, 11756354, 33174451, 93795220, 265565628];

%Inputs:

%n=Number of triangles in the polyiamond (Number of atoms missing in the nanopore)
n=8
remove_dangling=1;%Set remove_dangling=1 or 0, remove_dangling=1 means that dangling moeities and dangling bonds will be omitted and stable nanores will be saved, remove_dangling=0 means that dangling moeities and dangling bonds will not omitted and free polyiamonds without holes will be saved. This can be used to verify the code with the OEIS numbers, 
print=1;%print=1 or 0, if print=1, saves the coordinates in xyz file and viceversa
%Input Ranges of Properties (set to default values to include all
%nanopores)
%rlsf=lower limit of shape factor range
%rusf=upper limit of shape factor range
%rlma=lower limit of major axis range
%ruma=upper limit of major axis range
%rlmi=lower limit of minor axis range
%rumi=upper limit of minor axis range
rlsf=0;
rusf=inf;
rlma=0;
ruma=inf;
rlmi=0;
rumi=inf;

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

if n==1
    xyz(n,[])
else
    
[numpolys1,polys1] = fixed_inverted(n); %generates set of fixed polyiamonds with inverted triangle at origin

%Uncomment the followig lines to check the accuracy with respect to OEIS numbers of fixed polyiamonds

%[numpolysu,polysu] = fixed_upright(n);%generates set of fixed polyiamonds with upright triangle at origin
%fixed_polyiamonds=numpolysu+numpolys1
%oeis_fixed_polyiamonds=fixed_oeis(n)

fixed_polyiamonds_with_inverted_triangle_at_origin=numpolys1

data=zeros(numpolys1,3);

parfor i=1:numpolys1    
 numb=i;
 A=polys1{i};%A=representative coordinates of the polyiamond
 [ho,vrv,uid]=holes_db_elimination(A,remove_dangling);%Checks for nanopores with dangling bonds and non-bonded atoms or polyiamond with holes and vertices shared by 5 triangles
%ho indicates to preence of holes and vertices shared by 5 triangles. vrv is the vertex repetition vector. uid is the unique id of polyiamonds.     
 if ho==1 %Saves if no holes and vertices shared by 5 triangles are present          
  data(i,:)= [vrv uid i];
 end
end
end
%finds and deletes rows containg zeroes
In=find(data(:,3));
data=data(In,:);
%Use sortrows to organize data according to vertex repetition vector
data=sortrows(data);
fdata=flip(data(:,1));
[C,ia,ic] = unique(fdata);
si=[0;-(ia-length(data)-1)];%Extracts where a particular vertex repetition vector index ends
ind=[];
%Eliminating symmetrical isomers by comparison amongst polyiamonds having same vertex repetition vector
parfor i=2:size(si,1)
  indi=free_symmetry_op(data(si(i-1)+1:si(i),:));
  ind=[ind,indi];     
end
if remove_dangling==0%Saves free polyiamonds
    parfor w=1:size(ind,2)
     polys_ind{w}=polys1{ind(w)};
    end
    free_polyiamonds=numel(ind)
    OEIS_number=free_oeis(n)
    save("free_polyiamonds"+str+".mat",'ind','polys_ind')
    toc
else
     
%Finding properties, eliminating nanopores with dangling moeties and saving data and xyz files
ind_stable=[];
order=[];
minors=[];
majors=[];
shape_g=[];

parfor w=1:size(ind,2)
 A=polys1{ind(w)};%A=representative coordinates of the polyiamond
%This function calculates and returns the values of properties (major axis, minor axis and shapefactor), checks for presence of dangling moeities and prints the coordinates in xyz file if the calculated values fall in the range given as input).
   [sf,ma,mi,moeityr,inrange]= properties_moeity_xyz(A,rlsf,rusf,rlma,ruma,rlmi,rumi,n,print);
 if moeityr==0 && inrange==1
     %Saving properies
         ind_stable=[ind_stable;ind(w)];
    majors=[majors; ma];
    minors=[minors; mi];
    shape_g=[shape_g;sf];
 end
end

parfor i=1:numel(ind_stable)
   polys_ind{i}=polys1{ind_stable(i)};%Saving Reperesntative coordinates of stable nanopores 
end
    
 %Sorting tale based on Shape factor
 [shape_sorted,inde]=sort(shape_g);
 %Sorting tale based on Minor axis
 [minors_sorted,indemi]=sort(minors);
 %Sorting tale based on Major axis
 [majors_sorted,indema]=sort(majors);
format long g
tab=[inde shape_sorted];
str=string(n);

%Saving properties after sorting them based on their sorted shape factor values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted shape factor values
%4th and 5th column denotes to their corresponding major axis and minor axis values respectively
 writematrix([inde ind_stable(inde) shape_sorted majors(inde) minors(inde)],"shape_factor"+str+".xlsx",'WriteMode','append')
 
 %Saving properties after sorting them based on their minor axis values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted minor axis values
%4th and 5th column denotes to their corresponding major axis values and shape factor values respectively
 writematrix([indemi ind_stable(indemi) minors_sorted majors(indemi) shape_g(indemi)],"minor_axis"+str+".xlsx",'WriteMode','append')
  
%Saving properties after sorting them based on their major axis values.
%1st column,inde denotes to their indices in the ind 
%2nd column denotes to their indices in the initially generated set of fixed polyiamonds with inverted triangle at origin
%3rd column denotes the sorted major axis values
%4th and 5th column denotes to their corresponding minor axis values and shape factor values respectively
 writematrix([indema ind_stable(indema) majors_sorted minors(indema) shape_g(indema)],"major_axis"+str+".xlsx",'WriteMode','append')

 %Saving of Data
str=string(n);
number_of_stable_nanopores=numel(ind_stable)
save("stable_nanopores"+str+".mat",'ind_stable','polys_ind')
end
toc
