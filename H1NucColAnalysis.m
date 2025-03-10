clearvars -except  AnlzComp 
clc
cdoriginal=cd;
%%
[f_name1, path_name, filter_index] = uigetfile('.csv','Select first image of sequence');
cd(path_name);
DATAF1pN=readmatrix(f_name1,'Sheet','F=1pN');
DATAF1pNRed=readmatrix(f_name1,'Sheet','F=1pN_NucPos');
DATAF1pNFRET=readmatrix(f_name1,'Sheet','F=1pN_FRET');
DATAF1pNCy5Time=readmatrix(f_name1,'Sheet','F=1pN_Cy5PosLastTime');
cd(cdoriginal);
%%
[DATAGreenF1pN, DATARedF1pN] = SortDataPerNucNum(DATAF1pN,DATAF1pNRed);
[DATAFretF1pN] = SortDataPerNucNum(DATAF1pNFRET,DATAF1pNRed);
%%
close all
NucTimeSortedComp=[];
NF1pNComp=[];
for i=1:length(DATAGreenF1pN)
    [AnlzF1pN, MSDTableF1pN, MSDAveF1pN,DATASortedF1pN] = MSDAnls(DATAGreenF1pN{i});
    [AnlzF1pNFret, MSDTableF1pNFret, MSDAveF1pNFret,DATASortedF1pNFret] = MSDAnls(DATAFretF1pN{i});    
    [NucPos] = OctPosAnlz(DATARedF1pN{i})
    [NucTime] = OctLastTimeAnlz(DATANucTimeF1pN{i})
    NucPosSorted=sort(NucPos(:,2))
    NucTimeSorted=sortrows(NucTime,2)
    NucTimeSortedComp=[NucTimeSortedComp; NucTimeSorted(3)];
    L=size(NucPosSorted);
    NF1pNComp=[NF1pNComp;L(1)];
    Width=0.2;
    figure(i)
    hold on
    for j=1:L(1)
        plot([0 NucTimeSorted(j,3)],[NucPos(j,2) NucPos(j,2)],'r-');
        plot([0 NucTimeSorted(j,3)],[NucPos(j,2)+Width NucPos(j,2)+Width],'r--');
        plot([0 NucTimeSorted(j,3)],[NucPos(j,2)-Width NucPos(j,2)-Width],'r--');
    end
    for j=1:length(DATASortedF1pN)
        yy=[];
        yy = smooth(DATASortedF1pN{j}(:,3),'rlowess');
        yy2 = smooth(yy,'rlowess');
        yy2(:,2)=zeros(length(yy),1);
        yy2(:,3)=zeros(length(yy),1);
        yy2(:,4)=zeros(length(yy),1);
        for k=1:length(yy2)
            for m=1:L(1)
                if (NucPosSorted(m)+Width)<yy2(k)
                    yy2(k,2)=m;
                elseif abs(NucPosSorted(m)-yy2(k))<=Width
                    yy2(k,2)=L(1)+m;
                    yy2(k,4)=NucTimeSorted(m,3);
                end
            end
            yy2(k,3)=L(1);
        end
        figure(i)
        hold on
        yyaxis left
        plot(DATASortedF1pN{j}(:,2),yy2(:,1),'g.-') 
        DATASortedF1pN{j}(:,6:9)=yy2;
        NucInd=find(yy2(:,2)>L(1))
    end
    DATASortedF1pNComp{i}=DATASortedF1pN;
    DATASortedF1pNFretComp{i}=DATASortedF1pNFret;
    pause (1)
    imagename = strcat('ImageF1pN',num2str(i),'.tif');
    cd(path_name);
    saveas(gcf,imagename);
    cd(cdoriginal);
end
%%
cnt=1;
TawF1pN=[];
for i=1:length(DATASortedF1pNComp)
    
   for j=1:length(DATASortedF1pNComp{i})
        
%   for i=1
%       for j=3

        NucChange=DATASortedF1pNComp{i}{j}(:,7)-DATASortedF1pNComp{i}{j}(:,8);
        NucInd=find(NucChange>0);
        NucIndChange=[];
        for k=1:length(NucInd)-1
            NucIndChange(k,1)=NucInd(k+1)-NucInd(k);
        end
        NucChangeInd=find(NucIndChange>1);

        k=1; 
        while k<length(NucInd)
            k
            if NucIndChange(k)~=1
                if NucInd(k)~=1 
                    TawF1pN(cnt,:)=[i j NucInd(k) NucInd(k) 0 DATASortedF1pNComp{i}{j}(NucInd(k)+1,7)-DATASortedF1pNComp{i}{j}(NucInd(k)-1,7)];
                    k=k+1;
                    cnt=cnt+1
                    % figure(i)
                    % plot(DATASortedF1pNComp{i}{j}(NucInd(k),2),DATASortedF1pNComp{i}{j}(NucInd(k),6),'c.-') 
                    % pause(0.1)
                elseif NucInd(k)==1
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(k) 0 Inf];
                            cnt=cnt+1
                            k=k+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k),2),DATASortedF1pNComp{i}{j}(NucInd(k),6),'b.-') 
                            % pause(0.1)
                end
            else
                for m=k:length(NucIndChange)
                    if NucIndChange(m)~=1
                        if NucInd(k)~=1 && NucInd(m)~=length(NucChange)
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m) DATASortedF1pNComp{i}{j}(NucInd(m),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) DATASortedF1pNComp{i}{j}(NucInd(m)+1,7)-DATASortedF1pNComp{i}{j}(NucInd(k)-1,7)];
                            cnt=cnt+1
                            pause(0.1)
                            k=m+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'k.-') 
                            % pause(0.1)
                            break
                        elseif NucInd(k)==1
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m) DATASortedF1pNComp{i}{j}(NucInd(m),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) Inf];
                            cnt=cnt+1
                            pause(0.1)
                            k=m+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'b.-') 
                            % pause(0.1)
                            break
                        elseif NucInd(m)==length(NucChange)
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m) DATASortedF1pNComp{i}{j}(NucInd(m),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) Inf];
                            cnt=cnt+1
                            pause(0.1)
                            k=m+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'b.-') 
                            % pause(0.1)
                            break
                        end
                    elseif m==length(NucIndChange)
                        if NucInd(m+1)==length(NucChange)
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m+1) DATASortedF1pNComp{i}{j}(NucInd(m+1),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) Inf];
                            cnt=cnt+1
                            pause(0.1)
                            k=length(NucIndChange)+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'b.-') 
                            % pause(0.1)
                            break
                        elseif NucInd(k)==1
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m+1) DATASortedF1pNComp{i}{j}(NucInd(m+1),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) Inf];
                            cnt=cnt+1
                            pause(0.1)
                            k=length(NucIndChange)+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'b.-') 
                            % pause(0.1)
                            break
                        else
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m+1) DATASortedF1pNComp{i}{j}(NucInd(m+1),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) DATASortedF1pNComp{i}{j}(NucInd(m+1)+1,7)-DATASortedF1pNComp{i}{j}(NucInd(k)-1,7)];
                            cnt=cnt+1
                            pause(0.1)
                            k=length(NucIndChange)+1;
                            % figure(i)
                            % plot(DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),2),DATASortedF1pNComp{i}{j}(NucInd(k):NucInd(m),6),'k.-') 
                            % pause(0.1)
                            break
                        end
                    end
                end
            end
        end
    end
end
TawName = strcat('TawF1pN', num2str(Width*1000),'nmWindow','.xlsx');
cd(path_name);
xlswrite(TawName,TawF1pN);
cd(cdoriginal);
%%
TawF1pNPooled=TawF1pN;
TawF1pNPooled=TawF1pN;
cnt=1;
for i=1:length(TawF1pNPooled)-1
    if 0.4<=TawF1pNPooled(i,5)
        if 0.4<=TawF1pNPooled(i+1,5) && TawF1pNPooled(i,1)==TawF1pNPooled(i+1,1) && TawF1pNPooled(i,2)==TawF1pNPooled(i+1,2) && (TawF1pNPooled(i+1,3)-TawF1pNPooled(i,4))<=10
            TawF1pNPooled(i,:)=[TawF1pNPooled(i,1) TawF1pNPooled(i,2) TawF1pNPooled(i,3) TawF1pNPooled(i+1,4) TawF1pNPooled(i,5)+TawF1pNPooled(i+1,5) TawF1pNPooled(i,6)];
            TawF1pNPooled(i+1,:)=[];
            cnt=cnt+1
        end
    end
end
TawName = strcat('TawF1pN', num2str(Width*1000),'nmWindowPooled','.xlsx');
cd(path_name);
xlswrite(TawName,TawF1pNPooled);
cd(cdoriginal);

%% CPD plot of Colocalization times 
figure ()
[F1pNf,F1pNx]=ecdf(TawF1pNPooled(:,5));
plot(F1pNx,F1pNf,'r-')

%% Analyzing Nuc Bypass Probability

TawF1pNRand=TawF1pN(randperm(length(TawF1pN)),:);
LAnlzF1pN=size(TawF1pNRand(:,1));
TawF1pNRandG1=TawF1pNRand(1:LAnlzF1pN(1)/3,:);
TawF1pNRandG2=TawF1pNRand(LAnlzF1pN(1)/3:2*LAnlzF1pN(1)/3,:);
TawF1pNRandG3=TawF1pNRand(2*LAnlzF1pN(1)/3:end,:);

TawF1pNRandG1UnBound=TawF1pNRandG1(TawF1pNRandG1(:,5)<0.4,:);
TawF1pNRandG1Flt=TawF1pNRandG1UnBound( ~any( isinf( TawF1pNRandG1UnBound ),2),:);
TawF1pNRandG1FltCrossed = TawF1pNRandG1Flt( TawF1pNRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF1pNRandG1FltCrossed)/length(TawF1pNRandG1);
BoundRatioG1=1-(length(TawF1pNRandG1UnBound)/length(TawF1pNRandG1));
RefRatioG1=1-CrossRatioG1-BoundRatioG1;


TawF1pNRandG2UnBound=TawF1pNRandG2(TawF1pNRandG2(:,5)<0.4,:);
TawF1pNRandG2Flt=TawF1pNRandG2UnBound( ~any( isinf( TawF1pNRandG2UnBound ),2),:);
TawF1pNRandG2FltCrossed = TawF1pNRandG2Flt( TawF1pNRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF1pNRandG2FltCrossed)/length(TawF1pNRandG2);
BoundRatioG2=1-(length(TawF1pNRandG2UnBound)/length(TawF1pNRandG2));
RefRatioG2=1-CrossRatioG2-BoundRatioG2;


TawF1pNRandG3UnBound=TawF1pNRandG3(TawF1pNRandG3(:,5)<0.4,:);
TawF1pNRandG3Flt=TawF1pNRandG3UnBound( ~any( isinf( TawF1pNRandG3UnBound ),2),:);
TawF1pNRandG3FltCrossed = TawF1pNRandG3Flt( TawF1pNRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF1pNRandG3FltCrossed)/length(TawF1pNRandG3);
BoundRatioG3=1-(length(TawF1pNRandG3UnBound)/length(TawF1pNRandG3));
RefRatioG3=1-CrossRatioG3-BoundRatioG3;


CrossRatioF1pN=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF1pN=[CrossRatioF1pN; mean(CrossRatioF1pN); std(CrossRatioF1pN)/sqrt(3)]

BoundRatioF1pN=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF1pN=[BoundRatioF1pN; mean(BoundRatioF1pN); std(BoundRatioF1pN)/sqrt(3)]

RefRatioF1pN=[RefRatioG1; RefRatioG2; RefRatioG3];
RefRatioF1pN=[RefRatioF1pN; mean(RefRatioF1pN); std(RefRatioF1pN)/sqrt(3)]
%% Analyzing FRET

TawF1pNBound=TawF1pNPooled(0.4<=TawF1pNPooled(:,5),:);
%TawF1pNBoundCont=TawF1pNCont(0.4<=TawF1pNCont(:,5),:);
EFretF1pN=[];
AveEFretF1pN=[];
EFretCont=[];
cnt=1;
for c=1:length(TawF1pNBound)
    i=TawF1pNBound(c,1);
    j=TawF1pNBound(c,2);
    k=TawF1pNBound(c,3);
    m=TawF1pNBound(c,4);
    TimeCheck=[DATASortedF1pNComp{i}{j}(k:m,2)-DATASortedF1pNComp{i}{j}(k:m,9)];
    TimeCheckInd=find(TimeCheck<0)
    if ~isempty(TimeCheckInd)
        m=k+TimeCheckInd(end)-1;
        SignalGreen=DATASortedF1pNComp{i}{j}(k:m,4);
        SignalFret=DATASortedF1pNFretComp{i}{j}(k:m,4);
        SignalFretCorrected=SignalFret-0.1*SignalGreen--0.1*SignalFret;
        FretE=SignalFretCorrected./(SignalFretCorrected+SignalGreen);
        EFretF1pN{cnt}=[DATASortedF1pNComp{i}{j}(k:m,2) FretE];
        AveEFretF1pN(cnt,1)=mean(EFretF1pN{cnt}(:,2));
        AveEFretF1pN(cnt,2)=TawF1pNBound(c,5);
        AveEFretF1pN(cnt,3)=k;
        cnt=cnt+1;
    end
end

for i=1:length(AveEFretF1pN)
    if isnan(AveEFretF1pN(i,1))
        AveEFretF1pN(i,:)=[];
    end
end