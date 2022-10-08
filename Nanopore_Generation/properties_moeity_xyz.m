function [sf,ma,mi,moeity,inrange]= properties_moeity_xyz(A,rlsf,rusf,rlma,ruma,rlmi,rumi,n,print)
%The function calculates and returns the values of properties (major axis, minor axis and
%shapefactor), checks for presence of dangling moeities and prints the
%coordinates in xyz file if the calculated values fall in the range given
%as input).

%Inputs:
%A=representative coordinates of the polyiamond
%rlsf=lower limit of shape factor range
%rusf=upper limit of shape factor range
%rlma=lower limit of major axis range
%ruma=upper limit of major axis range
%rlmi=lower limit of minor axis range
%rumi=upper limit of minor axis range
%n=Number of triangles in the polyiamond (Number of atoms missing in the
%nanopore)
%print=1 or 0, if print=1, saves the coordinates in xyz file and vice
%versa

%Outputs:
%sf=shape factor
%ma=major axis
%mi=minor axis
%moeity=1 or 0, Moeity =1 denotes that dangling moeities are present moeity
%=0 indicates presence of no moeities
%inrange=1 or 0, Inrange =1 means all properties fall in the given range and vice versa. 

s=24.64;%side of triangle, multiplied by 10
moeity=0;  
inrange=0;
bl=zeros(n*3,2);%Initiating boundary of triangles in the polyiamond
cent=[];%Initiating centres of triangles in the polyiamond
sh=s/2;
c=2;
%Polyiamond generation
for i=1:size(A)

    pgons=nsidedpoly(3,'Center',[0 (A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)],'SideLength',s);
    if mod(A(i,1)+A(i,2),2)~=0%If sum of representative coordnates are odd, it denotes to an upright triangle.              
        pgone=nsidedpoly(3,'Center',[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5)],'SideLength',s);
        [bx,by]=boundary (pgone,1);
        b=unique([bx,by],'rows');
        bl(3*i-2:3*i,:)=b;
        cent=[cent;[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5),0]];
    else
        pgonr=rotate(pgons,180,[(A(i,1)*sh/2),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)]);%If sum of representative coordnates are even, it denotes to an inverted triangle, hence inverted triangle at origin 
        [bx,by]=boundary (pgonr,1);
        b=unique([bx,by],'rows');
        bl(3*i-2:3*i,:)=b;
        cent=[cent;[(A(i,1)*sh),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5),0]];
     end
end
b1=bl;
bli=unique(round(b1,4),'rows');
si= size(bli,1);
s=1.4226*10;%scaling C-C bondlength in Angstroms 
bl=zeros(si*6,2);
edge_list=[];
sf=[];
ma=[];
mi=[];
%Polyhex generation
for i=1:si   
      pgonh=nsidedpoly(6,'Center',bli(i,:),'SideLength',s);%hexagon generation
      pgonhr=rotate(pgonh,90,[bli(i,:)]);
      [bx,by]=boundary (pgonhr,1);
      b=unique(round([bx,by],1),'stable','rows');
      bl(6*i-5:6*i,:)=b;
      bl_un=unique(bl,'stable','rows');
      pos=[];
%Adding new vertices
      for j=1:6
         pos_i=find(ismember(bl_un,b(j,:),'rows')==1)';
         pos=[pos,pos_i];
      end
      e_i=pos;
%Adding Edges
      e_now=[e_i(1,1),e_i(1,2);e_i(1,2),e_i(1,3);e_i(1,3),e_i(1,4);e_i(1,4),e_i(1,5);e_i(1,5),e_i(1,6);e_i(1,6),e_i(1,1)];
      edge_list=[edge_list;e_now];
end
M=bl;
[~,~,ic] = unique(M, 'rows','stable')    ;       % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);  % Count Occurrences
maph = h(ic); % Map Occurrences To ‘ic’ Values
Result = unique([M, maph],'rows','stable');
%Finding inner points in polyhex (points shared by 3 hexagons)
n3=find(Result(:,3)==3);
centn=cent;
centn(:,3)=[];
cent=cent*0.1; 

%To make sure only n atoms are removed
indB = find(ismember(round(Result(n3,1:2),1),round(centn,1),'rows'));
n3=n3(indB);
%Removing edges containg inner points, to get the rim
for q=1:numel(n3)
    f1=find(edge_list(:,1)==n3(q));
    edge_list(f1,:)=[];
    f2=find(edge_list(:,2)==n3(q));
    edge_list(f2,:)=[];
end
le=length(edge_list);%Number of edges in the rim

%Arrangling edges such that each row in the edge list is sorted.
for i=1:le
    temp=sort(edge_list(i));
    edge_lists(i,:)=temp;%sorting
end
lm=length(unique(edge_lists,'rows'));%Number of edges in the rim excluding repeating oe self-intersecting edges like (v1,v2) and (v2,v1).
%If le and lm are not same, this indicates self-intersecting or doubly traversed edges showing presence of dangling moeities   
if le~=lm
    %Presence of moeity and returning
    moeity=1;
    sf=[];
    ma=[];
    mi=[];
    return;
else 
     moeity=0;%No moeities found
end

%In case no moeities are present, using the rim of the nanopore, properties
%are calculated

%Ordering of edges to get the rim
M2=edge_list;
edge_boundary=M2;
ver_b=M2(1,:);
edge_boundary(1,:)=[];
nume=numel(unique(M2));
for j=3:length(M2)+1
  if nume==numel(unique(ver_b))
     break;
  end
loc=find(edge_boundary(:,1)==ver_b(end));
if numel(loc)>1
    edges_loc=edge_boundary(loc,:);
    for r=1:numel(loc)
       fl=fliplr(edges_loc(r,:));
       if sum(ismember(edge_boundary,fl,'rows'))>0
           locn=loc(r);
           break;
       else
           locn=loc(r);
       end
    end
    if edge_boundary(locn,2)==ver_b(end-1)
        locn=loc(numel(loc)-1);
    end
    ver_b=[ver_b,edge_boundary(locn,2)];
    edge_boundary(locn,:)=[];
elseif numel(loc)==1
    ver_b=[ver_b,edge_boundary(loc,2)];
    edge_boundary(loc,:)=[];
end
 if isempty(loc)==1
    loc=find(edge_boundary(:,2)==ver_b(end));
 if numel(loc)>1
    edges_loc=edge_boundary(loc,:);
    for r=1:numel(loc)
       fl=flipr(edges_loc(r,:));
       if sum(ismember(edge_boundary,fl,'rows'))>0
           locn=loc(r); 
           break;
       else
           locn=loc(r);
       end
    end
    if edge_boundary(locn,1)==ver_b(end-1)
        locn=loc(numel(loc)-1);
    end
     ver_b=[ver_b,edge_boundary(locn,1)];
     edge_boundary(locn,:)=[];
  else
    ver_b=[ver_b,edge_boundary(loc,1)];
    edge_boundary(loc,:)=[];
  end
 end
end

ver_b_open=[ver_b,ver_b(1,1)];
P=Result(ver_b_open,1:2)';%Coordinates of rim
ver_b_open(end)=[];%Rim indices 

%Calculation of Major Axis
x=P(1,:);
y=P(2,:);
D = sqrt((x-x.').^2 + (y-y.').^2) ;% distance to/from all points
D(D<1E-6) = NaN; % remove self-distances
D=round(D,2);
Du=unique(D,'sort');
Du(isnan(Du))=[];
Du=sort(Du,'descend');
exc=1;
for i=1:numel(Du)
    if exc==0
        break;
    end
    Dmax=Du(i);
    [r,c]=find(D==Du(i));
    if i==1
        mxidxm=[r,c];
    end
    mxidx = [r,c];
 for l=1:length(mxidx)
   DVMx = x(mxidx(l,:));
   DVMy = y(mxidx(l,:));
   mxidxn=mxidx(l,:);
   lx1=linspace(DVMx(1,1),DVMx(1,2),100);
   ly1= linspace(DVMy(1,1),DVMy(1,2),100);
   [in1,on1]= inpolygon(lx1,ly1,x,y);%Checking if the major axis is going outside the boundary and to check if major axis is intersecting more than 2 points in the rim. 
   [intersx,~]=polyxpoly(DVMx,DVMy,x,y);%From mapping toolbox
   Dmaxx=Dmax;slopex1=x-DVMx(1,1);
   major_axis=Dmaxx*0.1;
   if  sum(on1)<=2  && all(in1)==1 && sum(in1)<=numel(in1) && length(intersx)<=2%&&sum(ismember(slope,[1,1],'rows'))==0 
      major_axis=Dmaxx*0.1;%major axis value
      exc=0;
      break;
   end
 end
end

if exc==1
    DVMx = x(mxidxm(1,:));
    DVMy = y(mxidxm(1,:));
    Dmaxx=Dmax(1);
end

%Minor Axis Calculation
ang = atan2d(diff(DVMy),diff(DVMx));
% each candidate perpendicular should have at least one endpoint on a vertex
% find distance from vertices V to opposite edge of polygon
perplen = zeros(numel(x),1);
pix = zeros(numel(x),1);
piy = zeros(numel(x),1);
for v = 1:numel(x)
    % skip if V belongs to major axis points
    if ismember(v,mxidxn); continue; end    
    % project a perpendicular from V
    px = x(v) + [1 -1]*Dmaxx*cosd(ang+90);
    py = y(v) + [1 -1]*Dmaxx*sind(ang+90);    
    % find intersections between perpendicular and polygon
    [xi,yi] = polyxpoly(x,y,px,py);    
    % find distance to points of intersection from V, maximize
    pid = [xi.'-x(v); yi.'-y(v)];
    if length(pid)==0
       continue;
    end
    [perplen(v) mpididx] = max(sqrt(pid(1,:).^2 + pid(2,:).^2));    
    % this is the other end of the candidate perpendicular
    pix(v) = xi(mpididx);
    piy(v) = yi(mpididx);
 end

% find maximal perpendicular vector
[~,idx] =sort(perplen,'descend');


for l=1:length(idx)
  DVPx = [x(idx(l,:)) pix(idx(l,:))];
  DVPy = [y(idx(l,:)) piy(idx(l,:))];
  lx=linspace(x(idx(l,:)),pix(idx(l,:)),100);
  ly= linspace(y(idx(l,:)),piy(idx(l,:)),100);
  [in,on] = inpolygon(lx,ly,x,y);
  [intersx,~]=polyxpoly(DVPx,DVPy,x,y);
  Dmin=hypot(DVPx(1,1)-DVPx(1,2),DVPy(1,1)-DVPy(1,2));
  minor_axis=Dmin*0.1;
  if all(in)==true && sum(on)>=2 && sum(in)<=numel(in)&& length(intersx)<=2
          minor_axis=Dmin*0.1; %Minor Axis value  
          break;
  end

end

A_hex=3^(1.5)*0.5*(1.4226)^2*si;%Area of polyhex
sf=sqrt(4*A_hex/pi)/major_axis;%Equivalent circular diameter divided by major axis

%If calculated properties fall in the input range, return the values of
%properties and save the coordinates if needed (if print =1)
if (rlsf<=sf) && (sf<=rusf) && (rlma<=major_axis) && (major_axis<=ruma) && (rlmi<=minor_axis) && (minor_axis<=rumi) 
    inrange=1;   %Indicates that properties fall in the given range
    ma=major_axis;
    mi=minor_axis;
    if print==1
       xyz(n,cent);%Saving Coordinates in xyyz file
    end
else
    inrange=0;
    sf=[];
    ma=[];
    mi=[];
         
end
end
