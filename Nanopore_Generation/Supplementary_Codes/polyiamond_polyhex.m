function polyiamond_polyhex(n,number)
%Input n(number of triangles in polyiamond) and the index of that in the saved cell array polys_ind.
%Plots polyiamonds and polyhexes.
close all;
str=string(n);
 load("stable_nanopores"+str+".mat");
 A=polys_ind{number};
    s=2.464*10;
sh=s/2;
c=2;
    fc=[0 0 1];
    bl=[];
    %Polyiamond plotting
figure;
for i=1:size(A)
    
        pgons=nsidedpoly(3,'Center',[0 (A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)],'SideLength',s);
       if mod(A(i,1)+A(i,2),2)~=0
         
        pgone=nsidedpoly(3,'Center',[(A(i,1)*sh), (A(i,2)*3^0.5*sh)+(((c-1)*sh)/3^0.5)],'SideLength',s);
         plot(pgone,'FaceColor','r','LineWidth',1);
[bx,by]=boundary (pgone,1);
b=unique([bx,by],'rows');
bl=[bl;b];

       else
         pgonr=rotate(pgons,180,[(A(i,1)*sh/2),(A(i,2)*3^0.5*sh)+((c*sh)/3^0.5)]);
         hold on;
                  
       plot(pgonr,'FaceColor','r','LineWidth',1);
       [bx,by]=boundary (pgonr,1);
b=unique([bx,by],'rows');
bl=[bl;b];

        
       end
end
axis equal;
     axis off;
    
    %Polyhex plotting
    % figure %Uncomment the line to make plots in separate figures
     bl=unique(round(bl,4),'rows');
si= size(bl,1);
    s=1.4226*10;
  sh=s/2;
    fc=[0 0 1];
   for i=1:si
   
        pgonh=nsidedpoly(6,'Center',bl(i,:),'SideLength',s);
      pgonhr=rotate(pgonh,90,[bl(i,:)]);
      hold on;
        plot(pgonhr,'FaceColor',fc,'LineWidth',1);
        
    end
axis equal;
     axis off;
     
