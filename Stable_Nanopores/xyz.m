function [Coordinates]=xyz(n,cent)
%This function prints the coordinates of nanopore to an xyz file
%Input:
%cent=cent of traingles in polyiamond or coordinates of missing atoms in
%the nanopore
%n=Number of triangles in the polyiamond (Number of atoms missing in the nanopore)
if n<=5
    a=8;
else 
    a=4;

end
if mod(n,2)~=0
   N=n+a-1;
else
   N=n+a;    
end
num=2*(N-1)*N+(2*N);

%The 2-D graphene layer coordinates can be coded by dividing it into to
%types of patterns; one starting from x=0(pattern1) and other from
%x=1.232(pattern2)
o=-1:2:num;
q=0.5:2:num; 
 
x1(1,1)=0;
x1(2,1)=0;
x2(1,1)=1.232;
x2(2,1)=1.232;
for i=3:2:2*N %X axis coordinates
    x1(i,1)=x1((i-1),1)+2.464;
    x1(i+1,1)=x1((i-1),1)+2.464;
    x2(i,1)=x2((i-1),1)+2.464;
    x2(i+1,1)=x2((i-1),1)+2.464;
 end

for j=1:N %Y axis Coordinates
    x3(j,1)=(j+o(1,j))*1.4226;
    x3(j,2)=(j+o(1,j)+1)*1.4226;
    x4(j,1)=(j+q(1,j))*1.4226;
    x4(j,2)=(j+q(1,j)+1)*1.4226;
end

X=zeros(num,3);

x=[x1' x2'];
y=x'; %X Axis Coordinates

for i=1:4*N:num %X axis coordinates are repeating every 4*N times.
   X(i:i-1+(4*N),1)=y(:,1) ; 
end
%X matrix now contains all X coordinates and Y and Z coordinates as zero.

for i=1:N %repmat function repeats the Y coodinates obtained earlier so that the pattern followed by Y coordinates can be made.
  s1(i,:)=repmat(x3(i,:),1,N); %Pattern1 Y Coordinates which shows recurrence after every 2*N eg.1-36,73-108.. when N=18
  s2(i,:)=repmat(x4(i,:),1,N); %Pattern2 Y Coordinates which shows recurrence after every 2*N eg.37-72,109-144.. when N=18
end

k=0;
for i=1:4*N:num %odd iterations of 2*N for pattern1
    k=k+1;
    X(i:i+2*N-1,2)=s1(k,:)'; %Filling Solution matrix with pattern 1 Y Coordinates.
end
l=0;
for i=2*N+1:4*N:num %even iterations of 2*N for pattern2
    l=l+1;
    X(i:i+2*N-1,2)=s2(l,:)'; %Filling Solution matrix with pattern 1 Y Coordinates.
end
Xr=round(X*100);

if n==1
    X(42, :) = [];
elseif n==2
    X(68, :)=[];
X(67, :) = [];

else    
    mi=min(cent(:,1));
    if mi<0     
       cent(:,1)=cent(:,1)-(mi*2);
    end
    cent(:,1)=cent(:,1)+(2.464*4);
    cent(:,2)=cent(:,2)+(1.4226*5);
    mi=min(cent(:,1));
%Finding and deleting atoms that is to be removed frm the sheet to get the deired nanopore.
    B=round(cent*100);
    [ia, ib] = ismember(Xr, B, 'rows');
    X(ia, :) = [];
end
Coordinates=X; %Solution matrix;
sC=size(Coordinates);
error=0;
if sC~=num-n
   error=1
   return;
end

if error~=1
% filename=sprintf('nanopore_coordinates_%d.xyz',n) ;
  filename=sprintf('stable_nanopores_%d.xyz',n) ;
  FileID=fopen(filename,'a+');
  fprintf(FileID,'%d\n',num-n);
  fprintf(FileID,'\n');
  fprintf(FileID,'C %f %f %f\n',Coordinates');
  fclose(FileID);
else 
  sC(1,1);   
end
end
