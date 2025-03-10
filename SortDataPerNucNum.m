function [DATAGreenSorted, DATARedSorted] = SortDataPerNucNum(DATAGreen,DATARed)
%%
ExpGreen=DATAGreen(:,1);
CounterGreen=DATAGreen(:,2);
TimeGreen=DATAGreen(:,3);
PositionGreen=DATAGreen(:,4);
PixelGreen=DATAGreen(:,5);
ExpRed=DATARed(:,1);
CounterRed=DATARed(:,2);
TimeRed=DATARed(:,3);
PositionRed=DATARed(:,4);
PixelRed=DATARed(:,5);

%%
%for i=1:length(Counter)
 %   b=regexp(Pixel(i),'\d+(\.)?(\d+)?','match')
  %  PixelNum(i,1)=str2double(strjoin([b{:}],''));
%end
%%
ExpIndxGreen=[];
cnt=1;
ExpCntGreen=ones(length(ExpGreen),1);
for j=1:length(ExpGreen)-1
    
    if ExpGreen(j+1)-ExpGreen(j)~=0
        
        
        ExpIndxGreen(cnt,1)=j+1;
        cnt=cnt+1;
        ExpCntGreen(j,1)=cnt;
    end
end
ExpIndxGreen=[1; ExpIndxGreen; length(ExpGreen)+1]; 
for i=1:length(ExpIndxGreen)-1
    DATAGreenSorted{i}=DATAGreen(ExpIndxGreen(i):ExpIndxGreen(i+1)-1,:);
end
%%
ExpIndxRed=[];
cnt=1;
ExpCntRed=ones(length(ExpRed),1);
for j=1:length(ExpRed)-1
    
    if ExpRed(j+1)-ExpRed(j)~=0
        
        
        ExpIndxRed(cnt,1)=j+1;
        cnt=cnt+1;
        ExpCntRed(j,1)=cnt;
    end
end
ExpIndxRed=[1; ExpIndxRed; length(ExpRed)+1]; 
for i=1:length(ExpIndxRed)-1
    DATARedSorted{i}=DATARed(ExpIndxRed(i):ExpIndxRed(i+1)-1,:);
end
end