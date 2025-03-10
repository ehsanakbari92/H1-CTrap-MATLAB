function [DATAGreenSorted] = SortDataPerDNAOnlyTest(DATAGreen)
%%
ExpGreen=DATAGreen(:,1);
CounterGreen=DATAGreen(:,2);
TimeGreen=DATAGreen(:,3);
PositionGreen=DATAGreen(:,4);
PixelGreen=DATAGreen(:,5);

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
end