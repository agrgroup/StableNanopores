function new_ind=free_symmetry_op(data)
%This function eliminates symmetrical isomers by comparison amongst polyiamonds having same vertex repetition vector
vrs_snow=string(data(:,2));%unique ids

ki=data(:,3);%Indices
new_ind=[];
original_ind_output=[];
vrs_comp=vrs_snow(1,:);
now1=strrep(strjoin(string(vrs_comp)),' ','');
set_compare1=strcat(now1,now1);%appending for rotation
nowi1=reverse(now1);
set_compare2=strcat(nowi1,nowi1);%appending for reflection
ind=ki(1);
     for j=2:numel(ki)  
        vrs_comp=strrep(strjoin(string(vrs_snow(j,:))),' ','');
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
    ind=[ind,ki(j)];
  set_compare1=[set_compare1,nows];% adding the rotation string to the set as any of its isomers are already there
  set_compare2=[set_compare2,nowsi];% adding the reflection string to the set as any of its isomers are already there

  break;
                 end
             end
         end
        
    end
  new_ind=[new_ind;ind];

end