function [Coordinates]=xyz_A(n,number,type,A)
%This function uses the same algorithm of xyz.m. This is just a different version of that to generate xyz files using representative coordinates (input A) and saved data.
s=24.64;
cent=[];
sh=s/2;
c=2;

for i=1:size(A)    
        pgons=nsidedpoly(3,'Center',[0 (A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)],'SideLength',s);
        if mod(A(i,1)+A(i,2),2)~=0     
           cent=[cent;[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5),0]];
        else
           cent=[cent;[(A(i,1)*sh),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5),0]];
        end
end
cent=cent*0.1;
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
else
    
 mi=min(cent(:,1));
 if mi<0
     
   cent(:,1)=cent(:,1)-(mi*2);
 end
 
 cent(:,1)=cent(:,1)+(2.464*4);
 cent(:,2)=cent(:,2)+(1.4226*5);
%Removes missing atoms from the sheet
 B=round(cent*100);
 [ia, ib] = ismember(Xr, B, 'rows');

 X(ia, :) = [];
end
Coordinates=X; %Solution matrix;
sC=size(Coordinates);
error=0;
if sC~=num-n
  error=1;
end
str=string(n);
str1=string(number);
if type==1
  filename="sf_order_"+str+"_"+str1+".xyz";
elseif type==2 
  filename="mi_order_"+str+"_"+str1+".xyz";
else
  filename="ma_order_"+str+"_"+str1+".xyz";
end
if error~=1
  FileID=fopen(filename,'a+');
  fprintf(FileID,'%d\n',num-n);
  fprintf(FileID,'\n');
  fprintf(FileID,'C %f %f %f\n',Coordinates');
  fclose(FileID);
  filename2=sprintf('min_max_%d.xyz',n) ;
  FileID=fopen(filename2,'a+');
  fprintf(FileID,'%d\n',num-n);
  fprintf(FileID,'\n');
  fprintf(FileID,'C %f %f %f\n',Coordinates');
  fclose(FileID);
else 
  sC(1,1);   
end
end
