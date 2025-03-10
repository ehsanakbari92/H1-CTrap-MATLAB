clearvars -except  AnlzComp 
clc
cdoriginal=cd;
%%
[f_name1, path_name, filter_index] = uigetfile('.csv','Select first image of sequence');
cd(path_name);
DATAF1pN=readmatrix(f_name1,'Sheet','F=1pN');
cd(cdoriginal);
%%
[DATAGreenF1pN] = SortDataPerDNAOnlyTest(DATAF1pN);
%%

for i=1:length(DATAGreenF1pN)
    i
    [AnlzF1pN, MSDTableF1pN, MSDAveF1pN,DATASortedF1pN] = MSDAnls(DATAGreenF1pN{i});
    [DATAF1pNFlt, DATAF1pNFlt2, DATAF1pNImmobile{i}, DATAF1pNMobile{i}]=FilterAlpha(AnlzF1pN,DATASortedF1pN,ImbFlt);
    %DATAF1pNFltImbMbl{i}=[DATAF1pNFlt3; DATAF1pNFlt4];
end

%%
NucPosComp=[];
for i=1:length(DATARedF1pN)
    [NucPos] = rand(4,1)*15;
    NucInd=0:3;
    NucPosSorted{i}=[NucInd' sort(NucPos(:,1))]
    NucPosComp=[NucPosComp; NucPosSorted{i}]
    NucNum(i,1)=length(NucPosSorted{i})
end
NucIndRand=randperm(length(DATARedF1pN))
%%
close all
Cnt=1;
for i=1:length(DATAF1pNMobile)
    if ~isempty(DATAF1pNMobile{i})
    [AnlzF1pN, MSDTableF1pN, MSDAveF1pN,DATASortedF1pN] = MSDAnls2(DATAF1pNMobile{i});
    L=size(NucPosSorted{NucIndRand(i)});
    Width=0.2;
    figure(i)
    for j=1:L(1)
        yline(NucPosSorted{NucIndRand(i)}(j,2),'r-');
        yline(NucPosSorted{NucIndRand(i)}(j,2)+Width,'r--');
        yline(NucPosSorted{NucIndRand(i)}(j,2)-Width,'r--');
    end
    for j=1:length(DATASortedF1pN)
        yy=[];
        yy = smooth(DATASortedF1pN{j}(:,3),'rlowess');
        yy2 = smooth(yy,'rlowess');
        yy2(:,2)=zeros(length(yy),1);
        for k=1:length(yy2)
            for m=1:L(1)
                if (NucPosSorted{NucIndRand(i)}(m)+Width)<yy2(k)
                    yy2(k,2)=m;
                elseif abs(NucPosSorted{NucIndRand(i)}(m)-yy2(k))<=Width
                    yy2(k,2)=L(1)+m;
                end
            end
            yy2(k,3)=L(1);
        end
        figure(i)
        hold on
        yyaxis left
        %plot(DATASortedF1pN{j}(:,2),DATASortedF1pN{j}(:,3),'g.-')
        %plot(DATASortedF1pN{j}(:,2),yy,'k.-')
        plot(DATASortedF1pN{j}(:,2),DATASortedF1pN{j}(:,3),'k.-') 
        plot(DATASortedF1pN{j}(:,2),yy(:,1),'r.-') 
        plot(DATASortedF1pN{j}(:,2),yy2(:,1),'g.-') 
        DATASortedF1pN{j}(:,6:8)=yy2;
        NucInd=find(yy2(:,2)>L(1))
    end
    DATASortedF1pNComp{i,:}=DATASortedF1pN;
    %Cnt=Cnt+1;
    %DATASortedF1pNFretComp{i,:}=DATASortedF1pNFret;
    pause (1)
    imagename = strcat('ImageF1pN',num2str(i),'.tif');
    cd(path_name);
    saveas(gcf,imagename);
    cd(cdoriginal);
    end
end

%%
cnt=1;
TawF1pN=[];
for i=1:length(DATASortedF1pNComp)
% 
    for j=1:length(DATASortedF1pNComp{i})
        
   % for i=1
   %     for j=5

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
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(k) 0 -Inf];
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
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m) DATASortedF1pNComp{i}{j}(NucInd(m),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) -Inf];
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
                            TawF1pN(cnt,:)=[i j NucInd(k) NucInd(m+1) DATASortedF1pNComp{i}{j}(NucInd(m+1),2)-DATASortedF1pNComp{i}{j}(NucInd(k),2) -Inf];
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
%%

%%
figure (1)
[F1pNf,F1pNx]=ecdf(TawF1pN(:,5));
plot(F1pNx,F1pNf,'k-')

%% Analyzing Nuc Bypass Probability

TawF1pNRand=TawF1pN(randperm(length(TawF1pN)),:);
LAnlzF1pN=size(TawF1pNRand(:,1));
TawF1pNRandG1=TawF1pNRand(1:LAnlzF1pN(1)/3,:);
TawF1pNRandG2=TawF1pNRand(LAnlzF1pN(1)/3:2*LAnlzF1pN(1)/3,:);
TawF1pNRandG3=TawF1pNRand(2*LAnlzF1pN(1)/3:end,:);
TawF1pNRandG1UnBound=TawF1pNRandG1(TawF1pNRandG1(:,5)<0.4,:);
TawF1pNRandG1Nuc=TawF1pNRandG1UnBound( any( isinf( TawF1pNRandG1UnBound ),2),:);
TawF1pNRandG1OffNuc=TawF1pNRandG1Nuc(find( 0<TawF1pNRandG1Nuc (:,6)),:);
TawF1pNRandG1OnNuc=TawF1pNRandG1Nuc(find( TawF1pNRandG1Nuc (:,6)<0),:);
TawF1pNRandG1Flt=TawF1pNRandG1UnBound( ~any( isinf( TawF1pNRandG1UnBound ),2),:);
TawF1pNRandG1FltCrossed = TawF1pNRandG1Flt( TawF1pNRandG1Flt(:,6)~=0,:);
TawF1pNRandG1FltReflected = TawF1pNRandG1Flt( TawF1pNRandG1Flt(:,6)==0,:);
NucOffRatioG1=size(TawF1pNRandG1OffNuc,1)/length(TawF1pNRandG1);
NucOnRatioG1=size(TawF1pNRandG1OnNuc,1)/length(TawF1pNRandG1);
CrossRatioG1=length(TawF1pNRandG1FltCrossed)/length(TawF1pNRandG1);
RefRatioG1=length(TawF1pNRandG1FltReflected)/length(TawF1pNRandG1);
BoundRatioG1=1-(length(TawF1pNRandG1UnBound)/length(TawF1pNRandG1));

TawF1pNRandG2UnBound=TawF1pNRandG2(TawF1pNRandG2(:,5)<0.4,:);
TawF1pNRandG2Nuc=TawF1pNRandG2UnBound( any( isinf( TawF1pNRandG2UnBound ),2),:);
TawF1pNRandG2OffNuc=TawF1pNRandG2Nuc(find( 0<TawF1pNRandG2Nuc (:,6)),:);
TawF1pNRandG2OnNuc=TawF1pNRandG2Nuc(find( TawF1pNRandG2Nuc (:,6)<0),:);
TawF1pNRandG2Flt=TawF1pNRandG2UnBound( ~any( isinf( TawF1pNRandG2UnBound ),2),:);
TawF1pNRandG2FltCrossed = TawF1pNRandG2Flt( TawF1pNRandG2Flt(:,6)~=0,:);
TawF1pNRandG2FltReflected = TawF1pNRandG2Flt( TawF1pNRandG2Flt(:,6)==0,:);
NucOffRatioG2=size(TawF1pNRandG2OffNuc,1)/length(TawF1pNRandG2);
NucOnRatioG2=size(TawF1pNRandG2OnNuc,1)/length(TawF1pNRandG2);
CrossRatioG2=length(TawF1pNRandG2FltCrossed)/length(TawF1pNRandG2);
RefRatioG2=length(TawF1pNRandG2FltReflected)/length(TawF1pNRandG2);
BoundRatioG2=1-(length(TawF1pNRandG2UnBound)/length(TawF1pNRandG2));

TawF1pNRandG3UnBound=TawF1pNRandG3(TawF1pNRandG3(:,5)<0.4,:);
TawF1pNRandG3Nuc=TawF1pNRandG3UnBound( any( isinf( TawF1pNRandG3UnBound ),2),:);
TawF1pNRandG3OffNuc=TawF1pNRandG3Nuc(find( 0<TawF1pNRandG3Nuc (:,6)),:);
TawF1pNRandG3OnNuc=TawF1pNRandG3Nuc(find( TawF1pNRandG3Nuc (:,6)<0),:);
TawF1pNRandG3Flt=TawF1pNRandG3UnBound( ~any( isinf( TawF1pNRandG3UnBound ),2),:);
TawF1pNRandG3FltCrossed = TawF1pNRandG3Flt( TawF1pNRandG3Flt(:,6)~=0,:);
TawF1pNRandG3FltReflected = TawF1pNRandG3Flt( TawF1pNRandG3Flt(:,6)==0,:);
NucOffRatioG3=size(TawF1pNRandG3OffNuc,1)/length(TawF1pNRandG3);
NucOnRatioG3=size(TawF1pNRandG3OnNuc,1)/length(TawF1pNRandG3);
CrossRatioG3=length(TawF1pNRandG3FltCrossed)/length(TawF1pNRandG3);
RefRatioG3=length(TawF1pNRandG3FltReflected)/length(TawF1pNRandG3);
BoundRatioG3=1-(length(TawF1pNRandG3UnBound)/length(TawF1pNRandG3));

NucOffRatioF1pN=[NucOffRatioG1; NucOffRatioG2; NucOffRatioG3];
NucOffRatioF1pN=[NucOffRatioF1pN; mean(NucOffRatioF1pN); std(NucOffRatioF1pN)/sqrt(3)]

NucOnRatioF1pN=[NucOnRatioG1; NucOnRatioG2; NucOnRatioG3];
NucOnRatioF1pN=[NucOnRatioF1pN; mean(NucOnRatioF1pN); std(NucOnRatioF1pN)/sqrt(3)]

CrossRatioF1pN=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF1pN=[CrossRatioF1pN; mean(CrossRatioF1pN); std(CrossRatioF1pN)/sqrt(3)]

RefRatioF1pN=[RefRatioG1; RefRatioG2; RefRatioG3];
RefRatioF1pN=[RefRatioF1pN; mean(RefRatioF1pN); std(RefRatioF1pN)/sqrt(3)]

BoundRatioF1pN=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF1pN=[BoundRatioF1pN; mean(BoundRatioF1pN); std(BoundRatioF1pN)/sqrt(3)]
%%

TawF1pNBound=TawF1pN(0.4<=TawF1pN(:,5),:);
EFret=[];
EFretCont=[];
for c=1:length(TawF1pNBound)
    i=TawF1pNBound(c,1);
    j=TawF1pNBound(c,2);
    k=TawF1pNBound(c,3);
    m=TawF1pNBound(c,4);
    SignalGreen=DATASortedF1pNComp{i}{j}(k:m,4);
    SignalFret=DATASortedF1pNFretComp{i}{j}(k:m,4);
    SignalFretCorrected=SignalFret-0.1*SignalGreen-0.1*SignalFret;
    FretE=SignalFretCorrected./(SignalFretCorrected+SignalGreen);
    EFretF1pN{c}=[DATASortedF1pNComp{i}{j}(k:m,2) SignalGreen SignalFretCorrected FretE];
    AveEFretF1pN(c,1)=mean(EFretF1pN{c}(:,4));
end

TawF1pNBoundCont=TawF1pNCont(0.4<=TawF1pNCont(:,5),:);
for c=1:length(TawF1pNBoundCont)
    i=TawF1pNBoundCont(c,1);
    j=TawF1pNBoundCont(c,2);
    k=TawF1pNBoundCont(c,3);
    m=TawF1pNBoundCont(c,4);
    SignalGreen=DATASortedF1pNComp{i}{j}(k:m,4);
    SignalFret=DATASortedF1pNFretComp{i}{j}(k:m,4);
    SignalFretCorrected=SignalFret-0.1*SignalGreen;
    EFretContF1pN{c}=SignalFretCorrected./(SignalFretCorrected+SignalGreen);
    AveEFretContFlpN(c,1)=mean(EFretContF1pN{c});
end

%%

TawF1pNContRand=TawF1pNCont(randperm(length(TawF1pNCont)),:);
LAnlzF1pNCont=size(TawF1pNContRand(:,1));
TawF1pNContRandG1=TawF1pNContRand(1:LAnlzF1pNCont(1)/3,:);
TawF1pNContRandG2=TawF1pNContRand(LAnlzF1pNCont(1)/3:2*LAnlzF1pNCont(1)/3,:);
TawF1pNContRandG3=TawF1pNContRand(2*LAnlzF1pNCont(1)/3:end,:);

TawF1pNContRandG1UnBound=TawF1pNContRandG1(TawF1pNContRandG1(:,5)<0.4,:);
TawF1pNContRandG1Flt=TawF1pNContRandG1UnBound( ~any( isinf( TawF1pNContRandG1UnBound ),2),:);
TawF1pNContRandG1FltCrossed = TawF1pNContRandG1Flt( TawF1pNContRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF1pNContRandG1FltCrossed)/length(TawF1pNContRandG1);
BoundRatioG1=1-(length(TawF1pNContRandG1UnBound)/length(TawF1pNContRandG1));

TawF1pNContRandG2UnBound=TawF1pNContRandG2(TawF1pNContRandG2(:,5)<0.4,:);
TawF1pNContRandG2Flt=TawF1pNContRandG2UnBound( ~any( isinf( TawF1pNContRandG2UnBound ),2),:);
TawF1pNContRandG2FltCrossed = TawF1pNContRandG2Flt( TawF1pNContRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF1pNContRandG2FltCrossed)/length(TawF1pNContRandG2);
BoundRatioG2=1-(length(TawF1pNContRandG2UnBound)/length(TawF1pNContRandG2));

TawF1pNContRandG3UnBound=TawF1pNContRandG3(TawF1pNContRandG3(:,5)<0.4,:);
TawF1pNContRandG3Flt=TawF1pNContRandG3UnBound( ~any( isinf( TawF1pNContRandG3UnBound ),2),:);
TawF1pNContRandG3FltCrossed = TawF1pNContRandG3Flt( TawF1pNContRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF1pNContRandG3FltCrossed)/length(TawF1pNContRandG3);
BoundRatioG3=1-(length(TawF1pNContRandG3UnBound)/length(TawF1pNContRandG3));

CrossRatioF1pNCont=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF1pNCont=[CrossRatioF1pNCont; mean(CrossRatioF1pNCont); std(CrossRatioF1pNCont)/sqrt(3)]

BoundRatioF1pNCont=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF1pNCont=[BoundRatioF1pNCont; mean(BoundRatioF1pNCont); std(BoundRatioF1pNCont)/sqrt(3)]

%% Analyzing Nuc Bypass Probability
TawF3pNRand=TawF3pN(randperm(length(TawF3pN)),:);
LAnlzF3pN=size(TawF3pNRand(:,1));
TawF3pNRandG1=TawF3pNRand(1:LAnlzF3pN(1)/3,:);
TawF3pNRandG2=TawF3pNRand(LAnlzF3pN(1)/3:2*LAnlzF3pN(1)/3,:);
TawF3pNRandG3=TawF3pNRand(2*LAnlzF3pN(1)/3:end,:);

TawF3pNRandG1UnBound=TawF3pNRandG1(TawF3pNRandG1(:,5)<0.4,:);
TawF3pNRandG1Flt=TawF3pNRandG1UnBound( ~any( isinf( TawF3pNRandG1UnBound ),2),:);
TawF3pNRandG1FltCrossed = TawF3pNRandG1Flt( TawF3pNRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF3pNRandG1FltCrossed)/length(TawF3pNRandG1);
BoundRatioG1=1-(length(TawF3pNRandG1UnBound)/length(TawF3pNRandG1));

TawF3pNRandG2UnBound=TawF3pNRandG2(TawF3pNRandG2(:,5)<0.4,:);
TawF3pNRandG2Flt=TawF3pNRandG2UnBound( ~any( isinf( TawF3pNRandG2UnBound ),2),:);
TawF3pNRandG2FltCrossed = TawF3pNRandG2Flt( TawF3pNRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF3pNRandG2FltCrossed)/length(TawF3pNRandG2);
BoundRatioG2=1-(length(TawF3pNRandG2UnBound)/length(TawF3pNRandG2));

TawF3pNRandG3UnBound=TawF3pNRandG3(TawF3pNRandG3(:,5)<0.4,:);
TawF3pNRandG3Flt=TawF3pNRandG3UnBound( ~any( isinf( TawF3pNRandG3UnBound ),2),:);
TawF3pNRandG3FltCrossed = TawF3pNRandG3Flt( TawF3pNRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF3pNRandG3FltCrossed)/length(TawF3pNRandG3);
BoundRatioG3=1-(length(TawF3pNRandG3UnBound)/length(TawF3pNRandG3));

CrossRatioF3pN=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF3pN=[CrossRatioF3pN; mean(CrossRatioF3pN); std(CrossRatioF3pN)/sqrt(3)]

BoundRatioF3pN=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF3pN=[BoundRatioF3pN; mean(BoundRatioF3pN); std(BoundRatioF3pN)/sqrt(3)]
%%
TawF3pNContRand=TawF3pNCont(randperm(length(TawF3pNCont)),:);
LAnlzF3pNCont=size(TawF3pNContRand(:,1));
TawF3pNContRandG1=TawF3pNContRand(1:LAnlzF3pNCont(1)/3,:);
TawF3pNContRandG2=TawF3pNContRand(LAnlzF3pNCont(1)/3:2*LAnlzF3pNCont(1)/3,:);
TawF3pNContRandG3=TawF3pNContRand(2*LAnlzF3pNCont(1)/3:end,:);

TawF3pNContRandG1UnBound=TawF3pNContRandG1(TawF3pNContRandG1(:,5)<0.4,:);
TawF3pNContRandG1Flt=TawF3pNContRandG1UnBound( ~any( isinf( TawF3pNContRandG1UnBound ),2),:);
TawF3pNContRandG1FltCrossed = TawF3pNContRandG1Flt( TawF3pNContRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF3pNContRandG1FltCrossed)/length(TawF3pNContRandG1);
BoundRatioG1=1-(length(TawF3pNContRandG1UnBound)/length(TawF3pNContRandG1));

TawF3pNContRandG2UnBound=TawF3pNContRandG2(TawF3pNContRandG2(:,5)<0.4,:);
TawF3pNContRandG2Flt=TawF3pNContRandG2UnBound( ~any( isinf( TawF3pNContRandG2UnBound ),2),:);
TawF3pNContRandG2FltCrossed = TawF3pNContRandG2Flt( TawF3pNContRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF3pNContRandG2FltCrossed)/length(TawF3pNContRandG2);
BoundRatioG2=1-(length(TawF3pNContRandG2UnBound)/length(TawF3pNContRandG2));

TawF3pNContRandG3UnBound=TawF3pNContRandG3(TawF3pNContRandG3(:,5)<0.4,:);
TawF3pNContRandG3Flt=TawF3pNContRandG3UnBound( ~any( isinf( TawF3pNContRandG3UnBound ),2),:);
TawF3pNContRandG3FltCrossed = TawF3pNContRandG3Flt( TawF3pNContRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF3pNContRandG3FltCrossed)/length(TawF3pNContRandG3);
BoundRatioG3=1-(length(TawF3pNContRandG3UnBound)/length(TawF3pNContRandG3));

CrossRatioF3pNCont=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF3pNCont=[CrossRatioF3pNCont; mean(CrossRatioF3pNCont); std(CrossRatioF3pNCont)/sqrt(3)]

BoundRatioF3pNCont=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF3pNCont=[BoundRatioF3pNCont; mean(BoundRatioF3pNCont); std(BoundRatioF3pNCont)/sqrt(3)]
%%
TawF3pNBound=TawF3pN(0.4<=TawF3pN(:,5),:);
EFret=[];
EFretCont=[];
for c=1:length(TawF3pNBound)
    i=TawF3pNBound(c,1);
    j=TawF3pNBound(c,2);
    k=TawF3pNBound(c,3);
    m=TawF3pNBound(c,4);
    SignalGreen=DATASortedF3pNComp{i}{j}(k:m,4);
    SignalFret=DATASortedF3pNFretComp{i}{j}(k:m,4);
    SignalFretCorrected=SignalFret-0.1*SignalGreen-0.1*SignalFret;
    FretE=SignalFretCorrected./(SignalFretCorrected+SignalGreen);
    EFretF3pN{c}=[DATASortedF3pNComp{i}{j}(k:m,2) SignalGreen SignalFretCorrected FretE];
    AveEFretF3pN(c,1)=mean(EFretF3pN{c}(:,2));
end

TawF3pNBoundCont=TawF3pNCont(0.4<=TawF3pNCont(:,5),:);
for c=1:length(TawF3pNBoundCont)
    i=TawF3pNBoundCont(c,1);
    j=TawF3pNBoundCont(c,2);
    k=TawF3pNBoundCont(c,3);
    m=TawF3pNBoundCont(c,4);
    SignalGreen=DATASortedF3pNComp{i}{j}(k:m,4);
    SignalFret=DATASortedF3pNFretComp{i}{j}(k:m,4);
    SignalFretCorrected=SignalFret-0.1*SignalGreen;
    EFretContF3pN{c}=SignalFretCorrected./(SignalFretCorrected+SignalGreen);
    AveEFretContF3pN(c,1)=mean(EFretContF3pN{c});
end
%%
TawF10pNRand=TawF10pN(randperm(length(TawF10pN)),:);
LAnlzF10pN=size(TawF10pNRand(:,1));
TawF10pNRandG1=TawF10pNRand(1:LAnlzF10pN(1)/3,:);
TawF10pNRandG2=TawF10pNRand(LAnlzF10pN(1)/3:2*LAnlzF10pN(1)/3,:);
TawF10pNRandG3=TawF10pNRand(2*LAnlzF10pN(1)/3:end,:);

TawF10pNRandG1UnBound=TawF10pNRandG1(TawF10pNRandG1(:,5)<0.4,:);
TawF10pNRandG1Flt=TawF10pNRandG1UnBound( ~any( isinf( TawF10pNRandG1UnBound ),2),:);
TawF10pNRandG1FltCrossed = TawF10pNRandG1Flt( TawF10pNRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF10pNRandG1FltCrossed)/length(TawF10pNRandG1);
BoundRatioG1=1-(length(TawF10pNRandG1UnBound)/length(TawF10pNRandG1));

TawF10pNRandG2UnBound=TawF10pNRandG2(TawF10pNRandG2(:,5)<0.4,:);
TawF10pNRandG2Flt=TawF10pNRandG2UnBound( ~any( isinf( TawF10pNRandG2UnBound ),2),:);
TawF10pNRandG2FltCrossed = TawF10pNRandG2Flt( TawF10pNRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF10pNRandG2FltCrossed)/length(TawF10pNRandG2);
BoundRatioG2=1-(length(TawF10pNRandG2UnBound)/length(TawF10pNRandG2));

TawF10pNRandG3UnBound=TawF10pNRandG3(TawF10pNRandG3(:,5)<0.4,:);
TawF10pNRandG3Flt=TawF10pNRandG3UnBound( ~any( isinf( TawF10pNRandG3UnBound ),2),:);
TawF10pNRandG3FltCrossed = TawF10pNRandG3Flt( TawF10pNRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF10pNRandG3FltCrossed)/length(TawF10pNRandG3);
BoundRatioG3=1-(length(TawF10pNRandG3UnBound)/length(TawF10pNRandG3));

CrossRatioF10pN=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF10pN=[CrossRatioF10pN; mean(CrossRatioF10pN); std(CrossRatioF10pN)/sqrt(3)]

BoundRatioF10pN=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF10pN=[BoundRatioF10pN; mean(BoundRatioF10pN); std(BoundRatioF10pN)/sqrt(3)]
%%
TawF10pNContRand=TawF10pNCont(randperm(length(TawF10pNCont)),:);
LAnlzF10pNCont=size(TawF10pNContRand(:,1));
TawF10pNContRandG1=TawF10pNContRand(1:LAnlzF10pNCont(1)/3,:);
TawF10pNContRandG2=TawF10pNContRand(LAnlzF10pNCont(1)/3:2*LAnlzF10pNCont(1)/3,:);
TawF10pNContRandG3=TawF10pNContRand(2*LAnlzF10pNCont(1)/3:end,:);

TawF10pNContRandG1UnBound=TawF10pNContRandG1(TawF10pNContRandG1(:,5)<0.4,:);
TawF10pNContRandG1Flt=TawF10pNContRandG1UnBound( ~any( isinf( TawF10pNContRandG1UnBound ),2),:);
TawF10pNContRandG1FltCrossed = TawF10pNContRandG1Flt( TawF10pNContRandG1Flt(:,6)~=0,:);
CrossRatioG1=length(TawF10pNContRandG1FltCrossed)/length(TawF10pNContRandG1);
BoundRatioG1=1-(length(TawF10pNContRandG1UnBound)/length(TawF10pNContRandG1));

TawF10pNContRandG2UnBound=TawF10pNContRandG2(TawF10pNContRandG2(:,5)<0.4,:);
TawF10pNContRandG2Flt=TawF10pNContRandG2UnBound( ~any( isinf( TawF10pNContRandG2UnBound ),2),:);
TawF10pNContRandG2FltCrossed = TawF10pNContRandG2Flt( TawF10pNContRandG2Flt(:,6)~=0,:);
CrossRatioG2=length(TawF10pNContRandG2FltCrossed)/length(TawF10pNContRandG2);
BoundRatioG2=1-(length(TawF10pNContRandG2UnBound)/length(TawF10pNContRandG2));

TawF10pNContRandG3UnBound=TawF10pNContRandG3(TawF10pNContRandG3(:,5)<0.4,:);
TawF10pNContRandG3Flt=TawF10pNContRandG3UnBound( ~any( isinf( TawF10pNContRandG3UnBound ),2),:);
TawF10pNContRandG3FltCrossed = TawF10pNContRandG3Flt( TawF10pNContRandG3Flt(:,6)~=0,:);
CrossRatioG3=length(TawF10pNContRandG3FltCrossed)/length(TawF10pNContRandG3);
BoundRatioG3=1-(length(TawF10pNContRandG3UnBound)/length(TawF10pNContRandG3));

CrossRatioF10pNCont=[CrossRatioG1; CrossRatioG2; CrossRatioG3];
CrossRatioF10pNCont=[CrossRatioF10pNCont; mean(CrossRatioF10pNCont); std(CrossRatioF10pNCont)/sqrt(3)]

BoundRatioF10pNCont=[BoundRatioG1; BoundRatioG2; BoundRatioG3];
BoundRatioF10pNCont=[BoundRatioF10pNCont; mean(BoundRatioF10pNCont); std(BoundRatioF10pNCont)/sqrt(3)]
%%
