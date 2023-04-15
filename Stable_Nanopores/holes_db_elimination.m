 function [ho,nu,sett_s,sett]=holes_db_elimination(A,remove_dangling)
 %This function identifies nanopores with dangling bonds and non-bonded atoms or polyiamond with holes and vertices shared by 5 triangles

 %Input
 %If remove_dangling =1, nanopores with dangling bonds needs to be identified
 %A=representative coordinates of the polyiamond
 
 %Output 
 %ho indicates to prsence of holes and vertices shared by 5 triangles. nu is the vertex repetition vector. sett_s is the unique id of polyiamonds.     
%ho=1 indicates no holes and no dangling bonds
n=size(A,1);%n=Number of triangles in the polyiamond (Number of atoms missing in the nanopore)
vr=n+2;
s=24.64;%side of triangle, multiplied by 10
sh=s/2;
c=2;
edge_list=[];
bl_in=[];
sett_s=[];
sett=[];
bl=zeros(n*3,2);%Initiating boundary of triangles in the polyiamond
for i=1:size(A)
       pgons=nsidedpoly(3,'Center',[0 (A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)],'SideLength',s);
       if mod(A(i,1)+A(i,2),2)~=0%If sum of representative coordnates are odd, it denotes to an upright triangle.         
         pgone=nsidedpoly(3,'Center',[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5)],'SideLength',s);
         [bx,by]=boundary (pgone,1);
         b=unique(round([bx,by],1),'rows');
         bl(3*i-2:3*i,:)=b;
         bl_un=unique(bl,'stable','rows');
         
         e_i=find(ismember(bl_un,b,'rows')==1)';%Finding new vertices to be added
         bl_in=[bl_in,e_i];%Adding new vertices to the list
         e_now=[e_i(1,1),e_i(1,2);e_i(1,2),e_i(1,3);e_i(1,1),e_i(1,3)];%Creating new edges
         edge_list=[edge_list;e_now];%Adding new edges
       else
         pgonr=rotate(pgons,180,[(A(i,1)*sh/2),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)]);%If sum of representative coordnates are even, it denotes to an inverted triangle, hence inverted triangle at origin 
         [bx,by]=boundary (pgonr,1);
         b=unique(round([bx,by],1),'rows');
         bl(3*i-2:3*i,:)=b;
         bl_un=unique(bl,'stable','rows');    
         
         e_i=find(ismember(bl_un,b,'rows')==1)';%Finding new vertices to be added
         bl_in=[bl_in,e_i];%Adding new vertices to the list
         e_now=[e_i(1,1),e_i(1,2);e_i(1,2),e_i(1,3);e_i(1,1),e_i(1,3)];%Creating new edges
         edge_list=[edge_list;e_now];%Adding new edges
       end
end
M=bl_in';
[~,~,ic] = unique(M,'stable')    ;       % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);  % Count Occurrences
maph = h(ic);% Map Occurrences To ‘ic’ Values
Result = unique([M, maph],'rows','stable');
n1=numel(find(Result(:,2)==1)); %Number of vertices shared by 1 triangle
n2=numel(find(Result(:,2)==2));%Number of vertices shared by 2 triangles
n3=numel(find(Result(:,2)==3));%Number of vertices shared by 3 triangles
n4=numel(find(Result(:,2)==4));%Number of vertices shared by 4 triangles
n5=numel(find(Result(:,2)==5));%Number of vertices shared by 5 triangles
n55=n5;
n6=numel(find(Result(:,2)==6));%Number of vertices shared by 6 triangles
if remove_dangling==0 %does not remove those nanopores with dangling bonds
    n55=0;%For cases with and without dangling bonds the value to be verified is set to 0, such that all cases are considered. 
end
if n6+n1+n2+n3+n4+n5+n6==vr && n55==0 %Checking for holes and dangling bonds
    ho=1;
    nu=[n1,n2,n3,n4,n5,n6];
    M2=edge_list;
    [~,~,ic] = unique(M2, 'rows','stable');% Unique Values By Row, Retaining Original Order
    h = accumarray(ic, 1); % Count Occurrences
    maph2 = h(ic);        % Map Occurrences To ‘ic’ Values
    M2(find(maph2==2),:)=[]; 
    edge_boundary=M2;
    ver_b=M2(1,:);
    edge_boundary(1,:)=[];
%Ordering of rim vertices and creation of unique id
    for j=3:length(M2)+1
        if  ver_b(1,1)==ver_b(end)
           break;
        end
        loc=find(edge_boundary(:,1)==ver_b(j-1));
    
        if isempty(loc)==1
          loc=find(edge_boundary(:,2)==ver_b(j-1));
          ver_b=[ver_b,edge_boundary(loc,1)];
          edge_boundary(loc,:)=[];
        else
          ver_b=[ver_b,edge_boundary(loc,2)];
          edge_boundary(loc,:)=[];
        end
     end
     ver_b_open=ver_b;%Ordered indices of boundary of polyiamond
     ver_b_open(end)=[];
     vrs=Result(ver_b_open,2);%Unique id
     sett=vrs';
     sett_ss=strrep(strjoin(string(sett)),' ','');
     sett_s=str2double(sett_ss);
     nus=strrep(strjoin(string(nu)),' ','');
     nu=nus;
else%Presence of dangling bonds or holes
     ho=0;
     nu=0;
end
end
