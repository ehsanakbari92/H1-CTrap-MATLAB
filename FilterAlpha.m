function [DATAFlt, DATAFlt2, DATAFlt3, DATAFlt4, DATAFlt5] = FilterAlpha(AnlzDATA,DATASorted,ImbFlt)
cnt=1;
cnt1=1;
cnt2=1;
cnt3=1;
cnt4=1;
DATAFlt=[];
DATAFlt2=[];
DATAFlt3=[];
DATAFlt4=[];
DATAFlt5=[];
for i=1:length(DATASorted)
    if 10<=length(DATASorted{i}(:,1))
        DATASortedFlt{cnt}=DATASorted{i};
        AnlzDATAFlt(cnt,:)=AnlzDATA(i,:);
        DATA=[(cnt-1)*ones(length(DATASorted{i}(:,1)),1) DATASorted{i}(:,2:4)];
        DATAFlt=[DATAFlt; DATA];
        cnt=cnt+1;
    end
end
for i=1:length(AnlzDATAFlt(:,1))
    if ImbFlt(2)<=AnlzDATAFlt(i,6) 
        DATA=[(cnt1-1)*ones(length(DATASortedFlt{i}(:,1)),1) DATASortedFlt{i}(:,2:4)];
        DATAFlt4=[DATAFlt4; DATA];
        cnt1=cnt1+1;
    elseif AnlzDATAFlt(i,9)<ImbFlt(1) 
        DATA=[(cnt2-1)*ones(length(DATASortedFlt{i}(:,1)),1) DATASortedFlt{i}(:,2:4)];
        DATAFlt3=[DATAFlt3; DATA];
        cnt2=cnt2+1;
    else
        DATA=[(cnt3-1)*ones(length(DATASortedFlt{i}(:,1)),1) DATASortedFlt{i}(:,2:4)];
        DATAFlt2=[DATAFlt2; DATA];
        cnt3=cnt3+1;
    end
end
for i=1:length(AnlzDATAFlt(:,1))
    if ImbFlt(2)<=AnlzDATAFlt(i,6) ||  ImbFlt(1)<=AnlzDATAFlt(i,9)
        DATA=[(cnt4-1)*ones(length(DATASortedFlt{i}(:,1)),1) DATASortedFlt{i}(:,2:4)];
        DATAFlt5=[DATAFlt5; DATA];
        cnt4=cnt4+1;
    end
end
end