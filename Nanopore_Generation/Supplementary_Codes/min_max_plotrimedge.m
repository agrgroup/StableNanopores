function [sf,ma,mi]= min_max_plotrimedge(A,n,k);
%This function uses the same algorithm used for generating all the properties and saving them. This is just a different version of that to generate and plot them separately
s=24.64;
bl=zeros(n*3,2);
cent=[];
sh=s/2;
c=2;
%Polyiamond generation
for i=1:size(A)    
       pgons=nsidedpoly(3,'Center',[0 (A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)],'SideLength',s);
       if mod(A(i,1)+A(i,2),2)~=0         
          pgone=nsidedpoly(3,'Center',[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5)],'SideLength',s);
          [bx,by]=boundary (pgone,1);
          b=unique([bx,by],'rows');
          bl(3*i-2:3*i,:)=b;
          cent=[cent;[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5),0]];
       else
          pgonr=rotate(pgons,180,[(A(i,1)*sh/2),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)]);       
          [bx,by]=boundary (pgonr,1);
          b=unique([bx,by],'rows');       
          bl(3*i-2:3*i,:)=b;
          cent=[cent;[(A(i,1)*sh),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5),0]];
       end
end
b1=bl;
bli=unique(round(b1,4),'rows');
[si,f]= size(bli);
s=1.4226*10;
fc=[1 0 0];
bl=zeros(si*6,2);
edge_list=[];
sf=[];
ma=[];
mi=[];
%polyhex generation
for i=1:si
      pgonh=nsidedpoly(6,'Center',bli(i,:),'SideLength',s);
      pgonhr=rotate(pgonh,90,[bli(i,:)]);
      [bx,by]=boundary (pgonhr,1);
      b=unique(round([bx,by],1),'stable','rows');
      bl(6*i-5:6*i,:)=b;
      bl_un=unique(bl,'stable','rows');
      pos=[];
      for j=1:6
         pos_i=find(ismember(bl_un,b(j,:),'rows')==1)';
         pos=[pos,pos_i];
      end
      e_i=pos;
      e_now=[e_i(1,1),e_i(1,2);e_i(1,2),e_i(1,3);e_i(1,3),e_i(1,4);e_i(1,4),e_i(1,5);e_i(1,5),e_i(1,6);e_i(1,6),e_i(1,1)];
      edge_list=[edge_list;e_now];
end
%Finding the rim of polyiamond and finding its connectivity
M=bl;
[Mu,ia,ic] = unique(M, 'rows','stable'); % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1);  % Count Occurrences
maph = h(ic); % Map Occurrences To ‘ic’ Values
Result = unique([M, maph],'rows','stable');
n3=find(Result(:,3)==3);
cent(:,3)=[];

indB = find(ismember(round(Result(n3,1:2),1),round(cent,1),'rows'));
n3=n3(indB);
for q=1:numel(n3)
    f1=find(edge_list(:,1)==n3(q));
    edge_list(f1,:)=[];
    f2=find(edge_list(:,2)==n3(q));
    edge_list(f2,:)=[];
end
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
points=Result(ver_b_open,1:2);
ver_b_open(end)=[];
P=points;
P=P';
Pu=unique(P','rows')';
x=P(1,:);
y=P(2,:);
cent=cent*0.1;  

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
    l=1;
    for l=1:length(mxidx)
        DVMx = x(mxidx(l,:));
        DVMy = y(mxidx(l,:));
        mxidxn=mxidx(l,:);
        lx1=linspace(DVMx(1,1),DVMx(1,2),100);
        ly1= linspace(DVMy(1,1),DVMy(1,2),100);
        [in1,on1]= inpolygon(lx1,ly1,x,y);
        [intersx,intersy]=polyxpoly(DVMx,DVMy,x,y);
        Dmaxx=Dmax;
        major_axis=Dmaxx*0.1;
        if  sum(on1)<=2  && all(in1)==1 && sum(in1)<=numel(in1) && length(intersx)<=2
           major_axis_v=Dmaxx*0.1-(1.7*2);
           major_axis=Dmaxx*0.1;
           exc=0;
           break;
        end
    end
end
% angle of DVM
if exc==1
    DVMx = x(mxidxm(1,:));
    DVMy = y(mxidxm(1,:));
    Dmaxx=Dmax(1);
end
ang = atan2d(diff(DVMy),diff(DVMx));

% each candidate perpendicular should have at least one endpoint on a vertex
% so consider all vertices V except those belonging to DVM
% find distance from V to opposite edge of polygon
perplen = zeros(numel(x),1);
pix = zeros(numel(x),1);
piy = zeros(numel(x),1);
for v = 1:numel(x)
% skip if V belongs to DVM
  
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
    [intersx,intersy]=polyxpoly(DVPx,DVPy,x,y);
    Dmin=hypot(DVPx(1,1)-DVPx(1,2),DVPy(1,1)-DVPy(1,2));
    minor_axis=Dmin*0.1;
    if all(in)==true && sum(on)>=2 && sum(in)<=numel(in)&& length(intersx)<=2 
          minor_axis=Dmin*0.1;   
          break;
    end
end

%Plotting

A_hex=3^(1.5)*0.5*(1.4226)^2*si;
sf=sqrt(4*A_hex/pi)/major_axis;
ma=major_axis;
mi=minor_axis;
st=string(n);
sfs=string(sf);
mas=string(major_axis);
mis=string(minor_axis);
figure;
plot(x,y,'b','Linewidth',3); hold on;
plot(DVMx,DVMy,'r','Linewidth',2);
hold on;
plot(DVPx,DVPy,'k','Linewidth',3); 
axis equal;
axis off;
figure(8);
subplot(2,3,k)
plot(x,y,'b','Linewidth',3); hold on;
plot(DVMx,DVMy,'r','Linewidth',2);
hold on;
plot(DVPx,DVPy,'k','Linewidth',3); 
axis equal;
axis off;
end
