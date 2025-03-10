function [Anlz, MSDTable, MSDAve, DataSorted2, DataSortedPerTest] = MSDAnls(DATA)
%%
Exp=DATA(:,1);
Counter=DATA(:,2);
Time=DATA(:,3);
Position=DATA(:,4);
Pixel=DATA(:,5);
%NucN=DATA(:,6);

%%
%for i=1:length(Counter)
 %   b=regexp(Pixel(i),'\d+(\.)?(\d+)?','match')
  %  PixelNum(i,1)=str2double(strjoin([b{:}],''));
%end
%%
ExpIndx=[];
cnt=1;
for j=1:length(Exp)-1
    
    if Exp(j+1)-Exp(j)~=0
        
        
        ExpIndx(cnt,1)=j+1;
        cnt=cnt+1;
    end
end
ExpIndx=[ExpIndx; length(Exp)]; 
if length(ExpIndx)~=1
    for i=1:length(ExpIndx)-1
        if i~=length(ExpIndx)-1
        Counter(ExpIndx(i):ExpIndx(i+1)-1,1)=Counter(ExpIndx(i):ExpIndx(i+1)-1,1)+Counter(ExpIndx(i)-1,1)+1;
        else
            Counter(ExpIndx(i):ExpIndx(i+1),1)=Counter(ExpIndx(i):ExpIndx(i+1),1)+Counter(ExpIndx(i)-1,1)+1;
        end
    end
end
%% Seperating the data based on each test and unique nuc number. 
% ExpIndx2=[1 ExpIndx; length(Exp)]; 
% if length(ExpIndx)~=2
%     for i=1:length(ExpIndx)-1
%         DataSortedPerTest{i}=DATA(ExpIndx(i):ExpIndx(i+1),:);
%     end
% else
%     DataSortedPerTest=DATA;
% end

%%
for j=0:Counter(end)

    cnt=1;
    for i=1:length(Counter)
        if Counter(i)==j
            DataSorted{j+1}(cnt,:)=[Counter(i) Time(i) Position(i) Pixel(i)];
            cnt=cnt+1;
        end
    end
end
%% Grouping the traces based on the nuc numbers;

% GN1=50;
% GN1=100;
% cntG1=1;
% cntG2=1;
% cntG3=1;
% for j=1:cnt-1
%     if DataSorted{cnt}(1,4)<=GN1
%         DataSortedG1{cntG1}=DataSorted{cnt};
%         cntG1=cntG1+1;
%     elseif DataSorted{cnt}(1,4)<=GN2
%         DataSortedG1{cntG2}=DataSorted{cnt};
%         cntG2=cntG2+1;
%     else
%         DataSortedG1{cntG2}=DataSorted{cnt};
%         cntG3=cntG3+1;
%     end
% end
%% MSD - Averaged across trajectories

for i=1:Counter(end)+1
    L=size(DataSorted{i});
    LengthInd(i,:)=[i;L(1)];
    i;
    for j=1:L(1)
        j;
        DataSorted{i}(j,5)=(DataSorted{i}(j,2)-DataSorted{i}(1,2));
        DataSorted{i}(j,6)=(DataSorted{i}(j,3)-DataSorted{i}(1,3))^2;
    end
end
SortLengthInd=sortrows(LengthInd,2);
for i=1:Counter(end)+1
    DataSorted2{i}=DataSorted{SortLengthInd(i,1)};
    %plot(DataSorted2{i}(:,2),DataSorted2{i}(:,3)-DataSorted2{i}(1,3),'r-')
    %hold on
end
for k=1:Counter(end)+1
    if k==1
        for i=1:SortLengthInd(k,2)
            TimeSum=0;
            MSDSum=0;
            MSDSTDMtx=[];
            for j=k:Counter(end)+1
                DataSorted2{j}(i,5);
                TimeSum=TimeSum+DataSorted2{j}(i,5);
                MSDSum=MSDSum+DataSorted2{j}(i,6);
                MSDSTDMtx(j)=DataSorted2{j}(i,6);
            end
             TimeAve(i)=TimeSum/(Counter(end)+1);
             MSD(i)=MSDSum/(Counter(end)+1);
             MSDSTD(i)=std(MSDSTDMtx)/sqrt(length(MSDSTDMtx));
        end
    else
        for i=SortLengthInd(k-1,2)+1:SortLengthInd(k,2)
            TimeSum=0;
            MSDSum=0;
            MSDSTDMtx=[];
            for j=k:Counter(end)+1
                TimeSum=TimeSum+DataSorted2{j}(i,5);
                MSDSum=MSDSum+DataSorted2{j}(i,6);
                MSDSTDMtx(j)=DataSorted2{j}(i,6);
            end
             TimeAve(i)=TimeSum/(Counter(end)+2-k);
             MSD(i)=MSDSum/(Counter(end)+2-k);
             MSDSTD(i)=std(MSDSTDMtx)/sqrt(length(MSDSTDMtx));
        end
    end
end
%Ind=1:Counter(end)+1;
%SortLengthInd=[Ind' SortLengthInd];
%%
%% MSD - Averaged for each trajectory

for k=1:Counter(end)+1
    for n=1:SortLengthInd(k,2)-1
        MSDSum=0;
        for i=1:SortLengthInd(k,2)-n
            MSDSum=MSDSum+(DataSorted2{k}(i+n,3)-DataSorted2{k}(i,3))^2;
        end
        if (SortLengthInd(k,2)-n+1)<(SortLengthInd(k,2)/10)
            break
        else
            MSDTable{k}(n+1,1)=DataSorted2{k}(n+1,2)-DataSorted2{k}(1,2);
            MSDTable{k}(n+1,2)=MSDSum/(SortLengthInd(k,2)-n);
        end
    end
end
%%
for i=1:Counter(end)+1
    SortLengthInd(i,3)=[length(MSDTable{i})];
end
%%
for k=1:Counter(end)+1
    if k==1
        for i=1:SortLengthInd(k,3)
            TimeSum=0;
            MSDSum=0;
            MSDSTDMtx2=[];
            MSDMax=MSDTable{1}(i,2);
            MSDMin=MSDTable{1}(i,2);
            for j=k:Counter(end)+1
                TimeSum=TimeSum+MSDTable{j}(i,1);
                MSDSum=MSDSum+MSDTable{j}(i,2);
                MSDSTDMtx2(j)=MSDTable{j}(i,2);
                if MSDTable{j}(i,2)<MSDMin
                    MSDMin=MSDTable{j}(i,2);
                end
                if MSDTable{j}(i,2)>MSDMax
                    MSDMax=MSDTable{j}(i,2);
                end
            end
             MSDAve(i,1)=TimeSum/(Counter(end)+1);
             MSDAve(i,2)=MSDSum/(Counter(end)+1);
             MSDAve(i,3)=std(MSDSTDMtx2)/sqrt(length(MSDSTDMtx2));
             MSDAve(i,4)=MSDMin;
             MSDAve(i,5)=MSDMax;
        end
    else
        for i=SortLengthInd(k-1,3)+1:SortLengthInd(k,3)
            TimeSum=0;
            MSDSum=0;
            MSDSTDMtx2=[];
            MSDMax=MSDTable{k}(i,2);
            MSDMin=MSDTable{k}(i,2);
            for j=k:Counter(end)+1
                TimeSum=TimeSum+MSDTable{j}(i,1);
                MSDSum=MSDSum+MSDTable{j}(i,2);
                MSDSTDMtx2(j)=MSDTable{j}(i,2);
                if MSDTable{j}(i,2)<MSDMin
                    MSDMin=MSDTable{j}(i,2);
                end
                if MSDTable{j}(i,2)>MSDMax
                    MSDMax=MSDTable{j}(i,2);
                end
            end
             MSDAve(i,1)=TimeSum/(Counter(end)+2-k);
             MSDAve(i,2)=MSDSum/(Counter(end)+2-k);
             MSDAve(i,3)=std(MSDSTDMtx2)/sqrt(length(MSDSTDMtx2));
             MSDAve(i,4)=MSDMin;
             MSDAve(i,5)=MSDMax;
        end
    end
end
%%
%LinFitLength=min(SortLengthInd(:,3));
%LinFitLength=3;
N=ceil(sqrt((Counter(end)+1)));
%hold on
for j=1:Counter(end)+1
    if length(MSDTable{j}(:,1))<10
        LinFitLength=length(MSDTable{j}(:,1));
    else
        LinFitLength=10;
    end
    %Anlz(j,:)=[DataSorted{j}(end,2)-DataSorted{j}(1,2) mean(DataSorted{j}(:,4))];
    %plot(MSDTable{j}(1:LinFitLength,1),MSDTable{j}(1:LinFitLength,2))
    F = fittype ( @(a,x) (a*x),'independent','x');
    F2 = fittype ( @(a,b,x) (a*x+b),'independent','x');
    [fitted_curve,gof] = fit(MSDTable{j}(1:LinFitLength,1), MSDTable{j}(1:LinFitLength,2), F, 'StartPoint', [1]);
    %subplot(N,N,j)
    %plot(MSDTable{j}(1:LinFitLength,1), MSDTable{j}(1:LinFitLength,2))
    MSDSlopeD=(fitted_curve.a)/2;
    Rsq=gof.rsquare;
    j;
    TF=isinf(log(MSDTable{j}(2:LinFitLength,2)));
    if sum(TF)==0 && 2<LinFitLength
        [fitted_curve2,gof] = fit(log(MSDTable{j}(2:LinFitLength,1)), log(MSDTable{j}(2:LinFitLength,2)), F2,  'StartPoint', [1, 0.1]);
        MSDAlpha=fitted_curve2.a;
        MSDAlphaD=exp((fitted_curve2.b))/2;
        RsqAlpha=gof.rsquare;
    else
        MSDAlpha=0;
        MSDAlphaD=0;
        RsqAlpha=0;
    end
        %P=polyfit(MSDTable{j}(1:LinFitLength,1),MSDTable{j}(1:LinFitLength,2),1);
    %Anlz(j,:)=[(j-1) DataSorted2{j}(end,2)-DataSorted2{j}(1,2) mean(DataSorted2{j}(:,4)) P(1)];
    Anlz(j,:)=[(j-1) DataSorted2{j}(end,2)-DataSorted2{j}(1,2) mean(DataSorted2{j}(:,4))  MSDSlopeD Rsq MSDAlpha MSDAlphaD RsqAlpha MSDTable{j}(LinFitLength,2)];
    end
%AnlzComp=[AnlzComp;Anlz];

end