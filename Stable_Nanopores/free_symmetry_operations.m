function [ind_noholes,ind_org]=free_symmetry_operations(data_uid,ki)
%This function eliminates symmetrical isomers by comparison amongst polyiamonds having same vertex repetition vector
strong=data_uid;
vrs_snow=strong;%unique ids
ki1=1:numel(ki);%Indices

ind_noholes=[];
ind_org=[];
vrs_comp=strong{1};
now1=strrep(strjoin(string(vrs_comp)),' ','');
set_compare1=strcat(now1,now1);%appending for rotation
nowi1=reverse(now1);
set_compare2=strcat(nowi1,nowi1);%appending for reflection
ind_org=ki(1);
ind=ki1(1);

for j=2:numel(ki1)  
        vrs_comp=strrep(strjoin(string(vrs_snow{j})),' ','');
        log1=0;
        log2=0;       
        for u=1:length(set_compare2)
            log1= contains(set_compare2(u),vrs_comp);%finding whether reflection  isomer already is counted
            if log1~=1
               log2= contains(set_compare1(u),vrs_comp);%finding whether rotation isomer already is counted
            end
            if log2==1 || log1==1%isomer aready counted once
                  break;
            else%isomer not counted, adding it to the set
                 if u==length(set_compare1)
                     now=vrs_comp;                  
                     nows=strcat(now,now); %appending for rotation
                     nowi=reverse(now);
                     nowsi=strcat(nowi,nowi);
                     ind=[ind,ki1(j)];
                     ind_org=[ind_org,ki(j)];
                     
                     set_compare1=[set_compare1,nows];% adding the rotation string to the set as any of its isomers are already there
                     set_compare2=[set_compare2,nowsi];% adding the reflection string to the set as any of its isomers are already there
                     break;
                 end
            end
         end
end
ind_noholes=ind';
ind_org=ind_org';
end
