clearvars -except  AnlzComp 
clc
cdoriginal=cd;
%%
[f_name1, path_name, filter_index] = uigetfile('.csv','Select first image of sequence');
cd(path_name);
DATAF1pN=readmatrix(f_name1,'Sheet','F=1pN');
DATAF3pN=readmatrix(f_name1,'Sheet','F=3pN');
DATAF10pN=readmatrix(f_name1,'Sheet','F=10pN');
cd(cdoriginal);
%%
[AnlzF1pN, MSDTableF1pN, MSDAveF1pN,DATASortedF1pN] = MSDAnls(DATAF1pN);
[AnlzF3pN, MSDTableF3pN, MSDAveF3pN,DATASortedF3pN] = MSDAnls(DATAF3pN);
[AnlzF10pN, MSDTableF10pN, MSDAveF10pN,DATASortedF10pN] = MSDAnls(DATAF10pN);
%%
ImbFlt1pN=[0.038	0.22];
[DATAF1pNTotal, DATAF1pNFiltered, DATAF1pNImmobile, DATAF1pNMobile]=FilterAlpha(AnlzF1pN,DATASortedF1pN,ImbFlt1pN);
[AnlzF1pNTotal, MSDTableF1pNTotal, MSDAveF1pNTotal,DATASortedF1pNTotal, SortLengthIndTotal] = MSDAnls2(DATAF1pNTotal);
[AnlzF1pNFiltered, MSDTableF1pNFiltered, MSDAveF1pNFiltered,DATASortedF1pNFiltered, SortLengthIndFlt2] = MSDAnls2(DATAF1pNFiltered);
[AnlzF1pNImmobile, MSDTableF1pNImmobile, MSDAveF1pNImmobile,DATASortedF1pNImmobile, SortLengthIndImmobile] = MSDAnls2(DATAF1pNImmobile);
[AnlzF1pNMobile, MSDTableF1pNMobile, MSDAveF1pNMobile,DATASortedF1pNMobile, SortLengthIndMobile] = MSDAnls2(DATAF1pNMobile);
%%
ImbFlt3pN=[0.025	0.22];
[DATAF3pNTotal, DATAF3pNFiltered, DATAF3pNImmobile, DATAF3pNMobile]=FilterAlpha(AnlzF3pN,DATASortedF3pN,ImbFlt1pN);
[AnlzF3pNTotal, MSDTableF3pNTotal, MSDAveF3pNTotal,DATASortedF3pNTotal, SortLengthIndTotal] = MSDAnls2(DATAF3pNTotal);
[AnlzF3pNFiltered, MSDTableF3pNFiltered, MSDAveF3pNFiltered,DATASortedF3pNFiltered, SortLengthIndFlt2] = MSDAnls2(DATAF3pNFiltered);
[AnlzF3pNImmobile, MSDTableF3pNImmobile, MSDAveF3pNImmobile,DATASortedF3pNImmobile, SortLengthIndImmobile] = MSDAnls2(DATAF3pNImmobile);
[AnlzF3pNMobile, MSDTableF3pNMobile, MSDAveF3pNMobile,DATASortedF3pNMobile, SortLengthIndMobile] = MSDAnls2(DATAF3pNMobile);
%%
ImbFlt10pN=[0.021	0.20];
[DATAF10pNTotal, DATAF10pNFiltered, DATAF10pNImmobile, DATAF10pNMobile]=FilterAlpha(AnlzF10pN,DATASortedF10pN,ImbFlt1pN);
[AnlzF10pNTotal, MSDTableF10pNTotal, MSDAveF10pNTotal,DATASortedF10pNTotal, SortLengthIndTotal] = MSDAnls2(DATAF10pNTotal);
[AnlzF10pNFiltered, MSDTableF10pNFiltered, MSDAveF10pNFiltered,DATASortedF10pNFiltered, SortLengthIndFlt2] = MSDAnls2(DATAF10pNFiltered);
[AnlzF10pNImmobile, MSDTableF10pNImmobile, MSDAveF10pNImmobile,DATASortedF10pNImmobile, SortLengthIndImmobile] = MSDAnls2(DATAF10pNImmobile);
[AnlzF10pNMobile, MSDTableF10pNMobile, MSDAveF10pNMobile,DATASortedF10pNMobile, SortLengthIndMobile] = MSDAnls2(DATAF10pNMobile);

%% Randomaizing each data set to calculate the mobile ratio
%LAnlzF1pN=length(DATASortedF1pN);
AnlzF1pNTotalRand=AnlzF1pNTotal(randperm(length(AnlzF1pNTotal)),:);
LAnlzF1pN=size(AnlzF1pNTotalRand);
AnlzF1pNFltRandG1=AnlzF1pNTotalRand(1:LAnlzF1pN(1)/3,:);
AnlzF1pNFltRandG2=AnlzF1pNTotalRand(LAnlzF1pN(1)/3:2*LAnlzF1pN(1)/3,:);
AnlzF1pNFltRandG3=AnlzF1pNTotalRand(2*LAnlzF1pN(1)/3:end,:);

loc = ImbFlt(2) <=  AnlzF1pNFltRandG1(:,6) ;
AnlzF1pNFltRandG1Diff=AnlzF1pNFltRandG1(loc,:);
loc = AnlzF1pNFltRandG1(:,6) < ImbFlt(2) ;
AnlzF1pNFltRandG1Imb1=AnlzF1pNFltRandG1(loc,:);
loc = AnlzF1pNFltRandG1Imb1(:,9) < ImbFlt(1) ;
AnlzF1pNFltRandG1Imb=AnlzF1pNFltRandG1Imb1(loc,:);
RatioG1=length(AnlzF1pNFltRandG1Diff(:,1))/(length(AnlzF1pNFltRandG1Diff(:,1))+length(AnlzF1pNFltRandG1Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF1pNFltRandG2(:,6) ;
AnlzF1pNFltRandG2Diff=AnlzF1pNFltRandG2(loc,:);
loc = AnlzF1pNFltRandG2(:,6) < ImbFlt(2) ;
AnlzF1pNFltRandG2Imb1=AnlzF1pNFltRandG2(loc,:);
loc = AnlzF1pNFltRandG2Imb1(:,9) < ImbFlt(1) ;
AnlzF1pNFltRandG2Imb=AnlzF1pNFltRandG2Imb1(loc,:);
RatioG2=length(AnlzF1pNFltRandG2Diff(:,1))/(length(AnlzF1pNFltRandG2Diff(:,1))+length(AnlzF1pNFltRandG2Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF1pNFltRandG3(:,6) ;
AnlzF1pNFltRandG3Diff=AnlzF1pNFltRandG3(loc,:);
loc = AnlzF1pNFltRandG3(:,6) < ImbFlt(2) ;
AnlzF1pNFltRandG3Imb1=AnlzF1pNFltRandG3(loc,:);
loc = AnlzF1pNFltRandG3Imb1(:,9) < ImbFlt(1) ;
AnlzF1pNFltRandG3Imb=AnlzF1pNFltRandG3Imb1(loc,:);
RatioG3=length(AnlzF1pNFltRandG3Diff(:,1))/(length(AnlzF1pNFltRandG3Diff(:,1))+length(AnlzF1pNFltRandG3Imb(:,1)));


RatioMobileF1pN=[RatioG1; RatioG2; RatioG3];
RatioMobileF1pN=[RatioDif1pN; mean(RatioDif1pN); std(RatioDif1pN)/sqrt(3)]
%%
AnlzF3pNTotalRand=AnlzF3pNTotal(randperm(length(AnlzF3pNTotal)),:);
LAnlzF3pN=size(AnlzF3pNTotalRand);
AnlzF3pNFltRandG1=AnlzF3pNTotalRand(1:LAnlzF3pN(1)/3,:);
AnlzF3pNFltRandG2=AnlzF3pNTotalRand(LAnlzF3pN(1)/3:2*LAnlzF3pN(1)/3,:);
AnlzF3pNFltRandG3=AnlzF3pNTotalRand(2*LAnlzF3pN(1)/3:end,:);

loc = ImbFlt(2) <=  AnlzF3pNFltRandG1(:,6) ;
AnlzF3pNFltRandG1Diff=AnlzF3pNFltRandG1(loc,:);
loc = AnlzF3pNFltRandG1(:,6) < ImbFlt(2) ;
AnlzF3pNFltRandG1Imb1=AnlzF3pNFltRandG1(loc,:);
loc = AnlzF3pNFltRandG1Imb1(:,9) < ImbFlt(1) ;
AnlzF3pNFltRandG1Imb=AnlzF3pNFltRandG1Imb1(loc,:);
RatioG1=length(AnlzF3pNFltRandG1Diff(:,1))/(length(AnlzF3pNFltRandG1Diff(:,1))+length(AnlzF3pNFltRandG1Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF3pNFltRandG2(:,6) ;
AnlzF3pNFltRandG2Diff=AnlzF3pNFltRandG2(loc,:);
loc = AnlzF3pNFltRandG2(:,6) < ImbFlt(2) ;
AnlzF3pNFltRandG2Imb1=AnlzF3pNFltRandG2(loc,:);
loc = AnlzF3pNFltRandG2Imb1(:,9) < ImbFlt(1) ;
AnlzF3pNFltRandG2Imb=AnlzF3pNFltRandG2Imb1(loc,:);
RatioG2=length(AnlzF3pNFltRandG2Diff(:,1))/(length(AnlzF3pNFltRandG2Diff(:,1))+length(AnlzF3pNFltRandG2Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF3pNFltRandG3(:,6) ;
AnlzF3pNFltRandG3Diff=AnlzF3pNFltRandG3(loc,:);
loc = AnlzF3pNFltRandG3(:,6) < ImbFlt(2) ;
AnlzF3pNFltRandG3Imb1=AnlzF3pNFltRandG3(loc,:);
loc = AnlzF3pNFltRandG3Imb1(:,9) < ImbFlt(1) ;
AnlzF3pNFltRandG3Imb=AnlzF3pNFltRandG3Imb1(loc,:);
RatioG3=length(AnlzF3pNFltRandG3Diff(:,1))/(length(AnlzF3pNFltRandG3Diff(:,1))+length(AnlzF3pNFltRandG3Imb(:,1)));


RatioMobileF3pN=[RatioG1; RatioG2; RatioG3];
RatioMobileF3pN=[RatioDiF3pN; mean(RatioDiF3pN); std(RatioDiF3pN)/sqrt(3)]
%% 
AnlzF10pNTotalRand=AnlzF10pNTotal(randperm(length(AnlzF10pNTotal)),:);
LAnlzF10pN=size(AnlzF10pNTotalRand);
AnlzF10pNFltRandG1=AnlzF10pNTotalRand(1:LAnlzF10pN(1)/3,:);
AnlzF10pNFltRandG2=AnlzF10pNTotalRand(LAnlzF10pN(1)/3:2*LAnlzF10pN(1)/3,:);
AnlzF10pNFltRandG3=AnlzF10pNTotalRand(2*LAnlzF10pN(1)/3:end,:);

loc = ImbFlt(2) <=  AnlzF10pNFltRandG1(:,6) ;
AnlzF10pNFltRandG1Diff=AnlzF10pNFltRandG1(loc,:);
loc = AnlzF10pNFltRandG1(:,6) < ImbFlt(2) ;
AnlzF10pNFltRandG1Imb1=AnlzF10pNFltRandG1(loc,:);
loc = AnlzF10pNFltRandG1Imb1(:,9) < ImbFlt(1) ;
AnlzF10pNFltRandG1Imb=AnlzF10pNFltRandG1Imb1(loc,:);
RatioG1=length(AnlzF10pNFltRandG1Diff(:,1))/(length(AnlzF10pNFltRandG1Diff(:,1))+length(AnlzF10pNFltRandG1Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF10pNFltRandG2(:,6) ;
AnlzF10pNFltRandG2Diff=AnlzF10pNFltRandG2(loc,:);
loc = AnlzF10pNFltRandG2(:,6) < ImbFlt(2) ;
AnlzF10pNFltRandG2Imb1=AnlzF10pNFltRandG2(loc,:);
loc = AnlzF10pNFltRandG2Imb1(:,9) < ImbFlt(1) ;
AnlzF10pNFltRandG2Imb=AnlzF10pNFltRandG2Imb1(loc,:);
RatioG2=length(AnlzF10pNFltRandG2Diff(:,1))/(length(AnlzF10pNFltRandG2Diff(:,1))+length(AnlzF10pNFltRandG2Imb(:,1)));

loc = ImbFlt(2) <=  AnlzF10pNFltRandG3(:,6) ;
AnlzF10pNFltRandG3Diff=AnlzF10pNFltRandG3(loc,:);
loc = AnlzF10pNFltRandG3(:,6) < ImbFlt(2) ;
AnlzF10pNFltRandG3Imb1=AnlzF10pNFltRandG3(loc,:);
loc = AnlzF10pNFltRandG3Imb1(:,9) < ImbFlt(1) ;
AnlzF10pNFltRandG3Imb=AnlzF10pNFltRandG3Imb1(loc,:);
RatioG3=length(AnlzF10pNFltRandG3Diff(:,1))/(length(AnlzF10pNFltRandG3Diff(:,1))+length(AnlzF10pNFltRandG3Imb(:,1)));


RatioMobileF10pN=[RatioG1; RatioG2; RatioG3];
RatioMobileF10pN=[RatioDiF10pN; mean(RatioDiF10pN); std(RatioDiF10pN)/sqrt(3)]
%% Creating a matrix pertaining the sorted dwell times under each category

T1pNMobile=AnlzF1pNMobile(:,2);
T1pNImmobile=AnlzF1pNImmobile(:,2);

T3pNMobile=AnlzF3pNMobile(:,2);
T3pNImmobile=AnlzF3pNImmobile(:,2);

T10pNMobile=AnlzF10pNMobile(:,2);
T10pNImmobile=AnlzF10pNImmobile(:,2);

%% Testing whether the dwell time CPD fits to a single exponential function 

Tau=T1pNMobile; % selecting which Tau data set to analyze 
TauRand=Tau(randperm(length(Tau)));
TauG1=TauRand(1:length(TauRand)/3);
TauG2=TauRand(length(TauRand)/3:2*length(TauRand)/3);
TauG3=TauRand(2*length(TauRand)/3:end);
figure(1)
hold on
[f,x]=ecdf(Tau);
fFlt=f(f<0.99); % choosing the CPD range to use for fitting
lFlt=length(fFlt);
xFlt=x(1:LFlt);
plot(xFlt, lFlt);
F2 = fittype ( @(A1,K1,x) (1- A1*exp(-K1*(x))),'independent','x');
[fitted_curve2,gof1] = fit(xFlt, fFlt, F2,'lower',[0, 0],'upper',[1, 1000],'StartPoint', [0.5, 10]);
plot(fitted_curve2,xFlt,fFlt)
TawFitAve=[fitted_curve2.K1.^-1 fitted_curve2.A1];

figure(2)
hold on
[fG1,xG1]=ecdf(TauG1);
fFltG1=fG1(fG1<0.99); % choosing the CPD range to use for fitting
lFltG1=length(fFltG1);
xFltG1=xG1(1:lFltG1);
[fitted_curve2,gof2] = fit(xFltG1, fFltG1, F2,'lower',[0, 0],'upper',[1 1000],'StartPoint', [0.5 10]);
plot(fitted_curve2,xFltG1,fFltG1)
TauFitG1=[fitted_curve2.K1.^-1 fitted_curve2.A1];

[fG2,xG2]=ecdf(TauG2);
fFltG2=fG2(fG2<0.99); % choosing the CPD range to use for fitting
lFltG2=length(fFltG2);
xFltG2=xG2(1:lFltG2);
[fitted_curve2,gof2] = fit(xFltG2, fFltG2, F2,'lower',[0, 0],'upper',[1 1000],'StartPoint', [0.5 10]);
plot(fitted_curve2,xFltG2,fFltG2)
TauFitG2=[fitted_curve2.K1.^-1 fitted_curve2.A1];

[fG3,xG3]=ecdf(TauG3);
fFltG3=fG3(fG3<0.99); % choosing the CPD range to use for fitting
lFltG3=length(fFltG3);
xFltG3=xG3(1:lFltG3);
[fitted_curve2,gof2] = fit(xFltG3, fFltG3, F2,'lower',[0, 0],'upper',[1 1000],'StartPoint', [0.5 10]);
plot(fitted_curve2,xFltG3,fFltG3)
TauFitG3=[fitted_curve2.K1.^-1 fitted_curve2.A1];

TauFitSExp=[TauFitG1; TauFitG2; TauFitG3]
TauFitSExp(4,:)=[mean(TauFitSExp(1:3,1)) mean(TauFitSExp(1:3,2))];
TauFitSExp(5,:)=[std(TauFitSExp(1:3,1))/sqrt(3) std(TauFitSExp(1:3,2))/sqrt(3)];
TauFitSExp(6,:)=TawFitAve;

%% Testing whether the dwell time CPD fits to a double exponential function 

Tau=T1pNMobile; % selecting which Tau data set to analyze 
TauRand=Tau(randperm(length(Tau)));
TauG1=TauRand(1:length(TauRand)/3);
TauG2=TauRand(length(TauRand)/3:2*length(TauRand)/3);
TauG3=TauRand(2*length(TauRand)/3:end);
figure(1)
hold on
[f,x]=ecdf(Tau);
fFlt=f(f<0.99); % choosing the CPD range to use for fitting
lFlt=length(fFlt);
xFlt=x(1:LFlt);
plot(xFlt, lFlt);
F = fittype ( @(A1,K1,K2,x) (1 - A1*exp(-K1*(x)) - ((1-A1)*exp(-K2*(x)))),'independent','x');
[fitted_curve1,gof1] = fit(xFlt, lFlt, F,'lower',[0, 0, 0],'upper',[1, 1000, 1000],'StartPoint', [0.5, 5, 25]);
plot(fitted_curve1,xFlt,lFlt)
TauAve=[fitted_curve1.K1.^-1 fitted_curve1.A1 fitted_curve1.K2.^-1];
if TauAve(3)<TauAve(1)
    TauAveSorted=[TauAve(3) TauAve(1) 1-TauAve(2)];
else
    TauAveSorted=[TauAve(1) TauAve(3) TauAve(2)];
end

% -------------------------------------
figure(2)
hold on
[fG1,xG1]=ecdf(TauG1);
fFltG1=fG1(fG1<0.99); % choosing the CPD range to use for fitting
lFltG1=length(fFltG1);
xFltG1=xG1(1:lFltG1);
[fitted_curve1,gof1] = fit(xFltG1, lFltG1, F,'lower',[0, 0, 0],'upper',[1, 1000, 1000]);
plot(fitted_curve1,xFltG1, lFltG1)
TauFitG1=[fitted_curve1.K1.^-1 fitted_curve1.A1 fitted_curve1.K2.^-1];
if TauFitG1(3)<TauFitG1(1)
    TauFitG1Sorted=[TauFitG1(3) TauFitG1(1) 1-TauFitG1(2)];
else
    TauFitG1Sorted=[TauFitG1(1) TauFitG1(3) TauFitG1(2)];
end

[fG2,xG2]=ecdf(TauG2);
fFltG2=fG2(fG2<0.99); % choosing the CPD range to use for fitting
lFltG2=length(fFltG2);
xFltG2=xG2(1:lFltG2);
[fitted_curve1,gof1] = fit(xFltG2, lFltG2, F,'lower',[0, 0, 0],'upper',[1, 1000, 1000]);
plot(fitted_curve1,xFltG2, lFltG2)
TauFitG2=[fitted_curve1.K1.^-1 fitted_curve1.A1 fitted_curve1.K2.^-1];
if TauFitG2(3)<TauFitG2(1)
    TauFitG2Sorted=[TauFitG2(3) TauFitG2(1) 1-TauFitG2(2)];
else
    TauFitG2Sorted=[TauFitG2(1) TauFitG2(3) TauFitG2(2)];
end

[fG3,xG3]=ecdf(TauG3);
fFltG3=fG3(fG3<0.99); % choosing the CPD range to use for fitting
lFltG3=length(fFltG3);
xFltG3=xG3(1:lFltG3);
[fitted_curve1,gof1] = fit(xFltG3, lFltG3, F,'lower',[0, 0, 0],'upper',[1, 1000, 1000]);
plot(fitted_curve1,xFltG3, lFltG3)
TauFitG3=[fitted_curve1.K1.^-1 fitted_curve1.A1 fitted_curve1.K2.^-1];
if TauFitG3(3)<TauFitG3(1)
    TauFitG3Sorted=[TauFitG3(3) TauFitG3(1) 1-TauFitG3(2)];
else
    TauFitG3Sorted=[TauFitG3(1) TauFitG3(3) TauFitG3(2)];
end

TauFitDExp=[TauFitG1Sorted; TauFitG2Sorted; TauFitG3Sorted];
TauFitDExp(4,:)=[mean(TauFitDExp(1:3,1)) mean(TauFitDExp(1:3,2)) mean(TauFitDExp(1:3,3)) ];
TauFitDExp(5,:)=[std(TauFitDExp(1:3,1))/sqrt(3) std(TauFitDExp(1:3,2))/sqrt(3) std(TauFitDExp(1:3,3))/sqrt(3) ];
TauFitDExp(6,:)=TauAveSorted;
%% scatter plot of the individaul mobile traces. 
figure()
hold on
D=AnlzF1pNMobile(:,7); % the estimated diffusion coefficient for each trace.
XPlot=1;
XplotD=XPlot*ones(length(D),1);
swarmchart(XplotD,D,'k')
plot(XPlot,mean(D),'r*')
DRand=D(randperm(length(D)));
DG1=DRand(1:length(DRand)/3);
DG2=DRand(length(DRand)/3:2*length(DRand)/3);
DG3=DRand(2*length(DRand)/3:end);
DMeanRand=[mean(DG1); mean(DG2); mean(DG3)];
DMeanRand=[DMeanRand; mean(DMeanRand); std(DMeanRand)/sqrt(3)];
%% Analyzing the average MSD to extract Alpha 

AveTime=MSDAveF1pNMobile(:,1);
AveTimeFlt=AveTime(AveTime<1);
AveLength=length(AveTimeFlt);
AveTimeFlt=AveTime(4:AveLength,1);
AveMSD=MSDAveF1pNMobile(4:Ave1pNLength,2);
AveSTD=MSDAveF1pNMobile(4:Ave1pNLength,3);
x_vector=[AveTimeFlt' fliplr(AveTimeFlt')];
y_vector=[(AveMSD'+AveSTD') fliplr((AveMSD'-AveSTD'))];
figure(1)
patch = fill(x_vector, y_vector, 'k');
hold on
plot(AveTimeFlt,AveMSD,'k.')
set(gca,'XScale','log')
set(gca,'YScale','log')
F3 = fittype ( @(a,b,x) (a*x+b),'independent','x');
[fitted_curve,gof] = fit(log(AveTimeFlt), log(AveMSD), F3);
Alpha=fitted_curve.a
DAve=exp(fitted_curve.b)/2;
plot(AveTimeFlt,AveMSD,'k')
plot(AveTimeFlt,2*DAve*AveTimeFlt.^Alpha,'r')
