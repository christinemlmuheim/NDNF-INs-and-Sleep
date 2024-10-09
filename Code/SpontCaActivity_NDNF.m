%%%% CM, October 2024, analysis script for spontanous calcium activity in NDNF INs in L1 
% of the Visual cortex. Analsysis script used for manuscript submitted to Natuer communications 
% October 2024, uploaded on Githut

%%% aim: explore calcium activity at baseline during various sleep stages
%%% and times and after increased sleep pressure. 

clc;
clear all;
close all hidden;
addpath('/Users/christinemuheim/Dropbox/NDNF_SleepPhenotyping/Functions')
%% ------ load the raw data ----
Data_path='/Users/christinemuheim/Dropbox/NDNF_SleepPhenotyping/RawData/SleepPhenotyping/SynchronizedFiles/'; 


%% GCaMP6f
fileC{1}='CM2308_m1_Calcium_synched.mat';
fileC{2}='CM2308_m4_Calcium_synched.mat';
fileC{3}='CM2308_m3_Calcium_synched.mat';

%%% data tables contain: timestamps, day (BL or SD), scored sleep stage,
%%% ROI calcium data as exported by IDPS (average ROI gray value per frame)

ROIsizeMax=200; % to creat an empty matrix with the maximal dimenions
%% ----- find BL and Sleep times, then find all bouts for NREM, REM and Wake

%%% loop over each file to import
%%% create a few empty variables first
numanim=3;
mat=NaN(numanim,ROIsizeMax);
RelRate=mat;
AbsRate=mat;
AbsRateBL=mat;
RelAmp=mat;
MeanBL=mat;
MeanSL=mat;

RelRateWake_Hourly=mat;
RelRateNREM_Hourly=mat;
RelRateREM_Hourly=mat;

RelRateWake_HourlySpont=mat;
RelRateNREM_HourlySpont=mat;
RelRateREM_HourlySpont=mat;
sr=10; % sampling rate for Ca imaging, 10FPS

epoch=900; % epoch size for EEG data in 4s epochs, 900 is 1h, 1800 would be 2h, 5400 would be 6h
for i=1:numanim
    clear  Ca
    Ca=load([Data_path,fileC{i}]);
    %%% animal ID
    AnimalID_SA=extractBefore(fileC{i},'_Calcium');
    clear ts 
    indxBL=find(strcmp(Ca.EEG_Ca_table.Condition,'BL'));
    indxSD=find(strcmp(Ca.EEG_Ca_table.Condition,'SD'));
    BL=Ca.EEG_Ca_table(indxBL,:);
    SD=Ca.EEG_Ca_table(indxSD,:);

    %%% find calcium frames during BL for each vigilance state
    %BL
    indxW=strcmp(BL.Stage,'W') | strcmp(BL.Stage,'W*');
    indxNR=strcmp(BL.Stage,'NR') | strcmp(BL.Stage,'N*');
    indxR=strcmp(BL.Stage,'R') | strcmp(BL.Stage,'R*');

    %SD
    indxWsd=strcmp(SD.Stage,'W') | strcmp(SD.Stage,'W*');
    indxNRsd=strcmp(SD.Stage,'NR') | strcmp(SD.Stage,'N*');
    indxRsd=strcmp(SD.Stage,'R') | strcmp(SD.Stage,'R*');

    
    %%get the sleep stages stats
    %BL
    BL.Stage=strrep(BL.Stage,'W*','W');
    BL.Stage=strrep(BL.Stage,'S*','S');
    BL.Stage=strrep(BL.Stage,'R*','R');


    T=groupcounts(BL,'Stage');
    BLstats(:,i)=T.Percent;

    %SD
    SD.Stage=strrep(SD.Stage,'W*','W');
    SD.Stage=strrep(SD.Stage,'S*','S');
    SD.Stage=strrep(SD.Stage,'R*','R');


    T=groupcounts(SD,'Stage');
    SDstats(:,i)=T.Percent;

    %get the stage specific tables, NaN all non wanted entries 
    %BL
    BLnew=timetable2table(BL(:,4:end));
    WAKE_BL=table2array(BLnew(:,vartype('numeric')));
    NREM_BL=table2array(BLnew(:,vartype('numeric')));
    REM_BL=table2array(BLnew(:,vartype('numeric')));
    
    %%% set other states to NaN
    WAKE_BL(~indxW,:)=NaN;
    NREM_BL(~indxNR,:)=NaN;
    REM_BL(~indxR,:)=NaN;

    %SD
    SDnew=timetable2table(SD(:,4:end));
    WAKE_SD=table2array(SDnew(:,vartype('numeric')));
    NREM_SD=table2array(SDnew(:,vartype('numeric')));
    REM_SD=table2array(SDnew(:,vartype('numeric')));
    
    %%% set other states to NaN

    WAKE_SD(~indxWsd,:)=NaN;
    NREM_SD(~indxNRsd,:)=NaN;
    REM_SD(~indxRsd,:)=NaN;

    %%% find the hour range, to get the rates per hour
    HourOfDay_BL=hour(BL.Time);
    HourOfDay_SD=hour(SD.Time);

    %%% find how many hours there are 
    [BL_hours,nBL,idx_BLhours]=RunLength_M(HourOfDay_BL);% find duration of of bouts
    [SD_hours,nSD,idx_SDhours]=RunLength_M(HourOfDay_SD);% find duration of of bouts
   
    Hours_TotBL(1:size(BL_hours,1),i)=BL_hours;
    Hours_TotSD(1:size(SD_hours,1),i)=SD_hours;
    %loop through every hour and sum calcium levels per vigilance state
    % normalized to the duration of the stage for that hour, there might be a lot of empty 
    %% Baseline
    for ii=1:1:(size(BL_hours,1))
        clear startepoch endepoch frameWake frameNREM frame REM medW medNR medR
        startepoch=(idx_BLhours(ii));
        endepoch=(idx_BLhours(ii)+(nBL(ii))-1);
        frameWake=(sum(indxW(startepoch:endepoch,:)));
        frameNREM=(sum(indxNR(startepoch:endepoch,:)));
        frameREM=(sum(indxR(startepoch:endepoch,:)));
       
        WAKEhourly_BL=WAKE_BL(startepoch:endepoch,:);
        NREMhourly_BL=NREM_BL(startepoch:endepoch,:);
        REMhourly_BL=REM_BL(startepoch:endepoch,:);
        
        %sum all calcium in that hour per state
 
        totW=(sum(WAKEhourly_BL(:,1:end),1,'omitnan'));
        totNR=(sum(NREMhourly_BL(:,1:end),1,'omitnan'));
        totR=(sum(REMhourly_BL(:,1:end),1,'omitnan'));

 
        medW=median(WAKEhourly_BL(:,1:end),1,'omitnan');
        medNR=median(NREMhourly_BL(:,1:end),1,'omitnan');
        medR=median(REMhourly_BL(:,1:end),1,'omitnan');
    
 
        totW=totW./(frameWake/sr);
        totNR=totNR./(frameNREM/sr);
        totR=totR./(frameREM/sr);


        CalcSum_W_BLhourly(ii,1:size(totW,2))=totW;      % summed and normlaized calcium per hour and state
        CalcSum_NR_BLhourly(ii,1:size(totNR,2))=totNR;      
        CalcSum_R_BLhourly(ii,1:size(totR,2))=totR;      

    end

    %% SleepDeprivation
    clear ii 
    for ii=1:1:(size(SD_hours,1))
        clear startepoch endepoch frameWake frameNREM frame REM medW medNR medR
        startepoch=(idx_SDhours(ii));
        endepoch=(idx_SDhours(ii)+(nSD(ii))-1);
        frameWake=(sum(indxWsd(startepoch:endepoch,:)));
        frameNREM=(sum(indxNRsd(startepoch:endepoch,:)));
        frameREM=(sum(indxRsd(startepoch:endepoch,:)));
       
        WAKEhourly_SD=WAKE_SD(startepoch:endepoch,:);
        NREMhourly_SD=NREM_SD(startepoch:endepoch,:);
        REMhourly_SD=REM_SD(startepoch:endepoch,:);

        %sum all calcium in that hour per state
        totW=(sum(WAKEhourly_SD(:,1:end),1,'omitnan'));
        totNR=(sum(NREMhourly_SD(:,1:end),1,'omitnan'));
        totR=(sum(REMhourly_SD(:,1:end),1,'omitnan'));   
     
        medW=median(WAKEhourly_SD(:,1:end),1,'omitnan');
        medNR=median(NREMhourly_SD(:,1:end),1,'omitnan');
        medR=median(REMhourly_SD(:,1:end),1,'omitnan');

        totW=totW./(frameWake/sr);
        totNR=totNR./(frameNREM/sr);
        totR=totR./(frameREM/sr);

        CalcSum_W_SDhourly(ii,1:size(totW,2))=totW;      
        CalcSum_NR_SDhourly(ii,1:size(totNR,2))=totNR;      
        CalcSum_R_SDhourly(ii,1:size(totR,2))=totR;      

    end
    %%%% save the summed hourly calcium for each animal,
    %BL
    CalcSum_Wake_BL(1:(size(CalcSum_W_BLhourly,1)),1:(size(CalcSum_W_BLhourly,2)),i)=CalcSum_W_BLhourly;
    CalcSum_NREM_BL(1:(size(CalcSum_NR_BLhourly,1)),1:(size(CalcSum_NR_BLhourly,2)),i)=CalcSum_NR_BLhourly;
    CalcSum_REM_BL(1:(size(CalcSum_R_BLhourly,1)),1:(size(CalcSum_R_BLhourly,2)),i)=CalcSum_R_BLhourly;
    %SD
    CalcSum_Wake_SD(1:(size(CalcSum_W_SDhourly,1)),1:(size(CalcSum_W_SDhourly,2)),i)=CalcSum_W_SDhourly;
    CalcSum_NREM_SD(1:(size(CalcSum_NR_SDhourly,1)),1:(size(CalcSum_NR_SDhourly,2)),i)=CalcSum_NR_SDhourly;
    CalcSum_REM_SD(1:(size(CalcSum_R_SDhourly,1)),1:(size(CalcSum_R_SDhourly,2)),i)=CalcSum_R_SDhourly;

    % total absolut rate for the entire time
    %BL
    AbsRate_BL(i,1:(size(CalcSum_W_BLhourly,2)),1)=sum(CalcSum_W_BLhourly,1,'omitnan');
    AbsRate_BL(i,1:(size(CalcSum_NR_BLhourly,2)),2)=sum(CalcSum_NR_BLhourly,1,'omitnan');
    AbsRate_BL(i,1:(size(CalcSum_R_BLhourly,2)),3)=sum(CalcSum_R_BLhourly,1,'omitnan');
 
    %SD
    AbsRate_SD(i,1:(size(CalcSum_W_SDhourly,2)),1)=sum(CalcSum_W_SDhourly,1,'omitnan');
    AbsRate_SD(i,1:(size(CalcSum_NR_SDhourly,2)),2)=sum(CalcSum_NR_SDhourly,1,'omitnan');
    AbsRate_SD(i,1:(size(CalcSum_R_SDhourly,2)),3)=sum(CalcSum_R_SDhourly,1,'omitnan');

 
    clear CalcSum_W_BLhourly CalcSum_NR_BLhourly CalcSum_R_BLhourly
    clear CalcSum_W_SDhourly CalcSum_NR_SDhourly CalcSum_R_SDhourly

end
 
%%% remove the 0, those are not real 0
CalcSum_Wake_BL(CalcSum_Wake_BL==0)=NaN;
CalcSum_NREM_BL(CalcSum_NREM_BL==0)=NaN;
CalcSum_REM_BL(CalcSum_REM_BL==0)=NaN;

CalcSum_Wake_SD(CalcSum_Wake_SD==0)=NaN;
CalcSum_NREM_SD(CalcSum_NREM_SD==0)=NaN;
CalcSum_REM_SD(CalcSum_REM_SD==0)=NaN;



%% separte them into NREM and REM active ROIs and plot in 6h intervals
% make means of the median activity to find in what state each ROI is more active acros 24h BL

CalcSum_Wake_BLmed=mean(CalcSum_Wake_BL,1,'omitnan');
CalcSum_NREM_BLmed=mean(CalcSum_NREM_BL,1,'omitnan');
CalcSum_REM_BLmed=mean(CalcSum_REM_BL,1,'omitnan');

MedianMat=[CalcSum_Wake_BLmed;CalcSum_NREM_BLmed;CalcSum_REM_BLmed]


clear i ii x
for ii=1:size(MedianMat,3)
    [~,x]=(max(MedianMat(:,:,ii))); % problem with NaN here, need to remove the later
    a=find(isnan(MedianMat(1,:,ii))); % to make sure only ROIs with actual values are integreated here
    x(:,a)=NaN; %find Nan location and make sure that the max matrix is also NaN there
    MaxMat(1,1:(size(x,2)),ii)=x;
end

%get the counts and make a pie chart
B=histcounts(MaxMat);
labels=({'Wake','NREM','REM'});

F1=figure  %% Figure 1B, manuscript
pie(B)
legend(labels)
title('Prefered Stage')

%create and index for Wake prefered, NREM prefered or REM prefered (based
%on 24h baseline median)

%Create empty matrix to populate
WakeON=NaN(size(MedianMat,1),size(MedianMat,2), numanim);
NREMON=NaN(size(MedianMat,1),size(MedianMat,2), numanim);
REMON=NaN(size(MedianMat,1),size(MedianMat,2), numanim);

clear i

for i=1:numanim
indxWpref=find(MaxMat(:,:,i)==1);
WakeON(:,indxWpref,i)=MedianMat(:,indxWpref,i)

indxNpref=find(MaxMat(:,:,i)==2);
NREMON(:,indxNpref,i)=MedianMat(:,indxNpref,i)

indxRpref=find(MaxMat(:,:,i)==3);
REMON(:,indxRpref,i)=MedianMat(:,indxRpref,i)
end

%%% remove missing values
A=(size(WakeON,2))*size(WakeON,3);
WakeON_rs=rmmissing(reshape(WakeON,[3,A]),2);
NREMON_rs=rmmissing(reshape(NREMON,[3,A]),2);
REMON_rs=rmmissing(reshape(REMON,[3,A]),2);


%% Scatter plot with median,  remove outliers first, values >2 SD removed
WakeON_rs_rm=rmoutliers(WakeON_rs');
NREMON_rs_rm=rmoutliers(NREMON_rs');
REMON_rs_rm=rmoutliers(REMON_rs');
clear x x1 x2
x(:,1)=0.5 + (0.5+0.5)*rand((size(WakeON_rs_rm,1)),1);
x(:,2)=2 + (0.5+0.5)*rand((size(WakeON_rs_rm,1)),1);
x(:,3)=3.5 + (0.5+0.5)*rand((size(WakeON_rs_rm,1)),1);

x1(:,1)=0.5 + (0.5+0.5)*rand((size(NREMON_rs_rm,1)),1);
x1(:,2)=2 + (0.5+0.5)*rand((size(NREMON_rs_rm,1)),1);
x1(:,3)=3.5 + (0.5+0.5)*rand((size(NREMON_rs_rm,1)),1);

x2(:,1)=0.5 + (0.5+0.5)*rand((size(REMON_rs_rm,1)),1);
x2(:,2)=2 + (0.5+0.5)*rand((size(REMON_rs_rm,1)),1);
x2(:,3)=3.5 + (0.5+0.5)*rand((size(REMON_rs_rm,1)),1);

c=[0.2 0.4 0.8]; % color for the fills
sz=60; %size of the circles

% median
medW=median(WakeON_rs_rm);
medN=median(NREMON_rs_rm);
medR=median(REMON_rs_rm);

F2a=figure
subplot(131)
scatter(x,WakeON_rs_rm,sz,'filled','MarkerFaceColor','#0072BD','MarkerEdgeColor','k'); hold on
l1=line([0.5,1.5],[medW(1),medW(1)],'Color','k','LineWidth',2); hold on
l2=line([2,3],[medW(2),medW(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4.5],[medW(3),medW(3)],'Color','k','LineWidth',2); hold on
ylim([-0.2 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
title('Wake-ON')
ylabel('Median Ca^{2+} activity at BL')
subplot(132)
scatter(x1,NREMON_rs_rm,sz,'filled','MarkerFaceColor','#EDB120','MarkerEdgeColor','k')
l1=line([0.5,1.5],[medN(1),medN(1)],'Color','k','LineWidth',2); hold on
l2=line([2,3],[medN(2),medN(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4.5],[medN(3),medN(3)],'Color','k','LineWidth',2); hold on
ylim([-0.2 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
title('NREM-ON')
subplot(133)
scatter(x2,REMON_rs_rm,sz,'filled','MarkerFaceColor','#77AC30','MarkerEdgeColor','k')
l1=line([0.5,1.5],[medR(1),medR(1)],'Color','k','LineWidth',2); hold on
l2=line([2,3],[medR(2),medR(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4.5],[medR(3),medR(3)],'Color','k','LineWidth',2); hold on
ylim([-0.2 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
title('REM-ON')

% run Anderson Darling test for nomality
clear i
for i=1:size(WakeON_rs_rm,2);
    clear p1 p2 p3 h
[h,p1] = adtest(WakeON_rs_rm(:,i));
[h,p2] = adtest(NREMON_rs_rm(:,i));
[h,p3] = adtest(REMON_rs_rm(:,i));
pWBL24(:,i)=p1
pNBL24(:,i)=p2
pRBL24(:,i)=p3
end

% Kruskal Wallis anova, not all datasets are parametric, aka the KS tests
[pw,pw,statsw]=kruskalwallis(WakeON_rs_rm);
[hn,pn,statsn]=kruskalwallis(NREMON_rs_rm);
[hr,pr,statsr]=kruskalwallis(REMON_rs_rm);

% protection from multiple comparsion with Tukey-Kramer, Dunn-Sidak or
% Bonferroni!!
pw1=multcompare(statsw,'CriticalValueType','bonferroni')
pn1=multcompare(statsn,'CriticalValueType','bonferroni')
pr1=multcompare(statsr,'CriticalValueType','bonferroni')

% rem multcompare gives exactly 0, update how many pos. you want in the
% output
fmt=['%5d%5d' repmat('%12.4f',1,3) '%12.3e\n'];
fprintf(fmt,pw1.')
fprintf(fmt,pn1.')
fprintf(fmt,pr1.')

%% Split the BL data in 6h intervals based on the median activity
% to compare to SD later on
%ZT0-6
CalcSum_Wake_BLmed_Int1=median(CalcSum_Wake_BL(1:6,:,:),1,'omitnan');
CalcSum_NREM_BLmed_Int1=median(CalcSum_NREM_BL(1:6,:,:),1,'omitnan');
CalcSum_REM_BLmed_Int1=median(CalcSum_REM_BL(1:6,:,:),1,'omitnan');

MedianMat_Int1=[CalcSum_Wake_BLmed_Int1;CalcSum_NREM_BLmed_Int1;CalcSum_REM_BLmed_Int1];

%ZT7-12
CalcSum_Wake_BLmed_Int2=median(CalcSum_Wake_BL(7:12,:,:),1,'omitnan');
CalcSum_NREM_BLmed_Int2=median(CalcSum_NREM_BL(7:12,:,:),1,'omitnan');
CalcSum_REM_BLmed_Int2=median(CalcSum_REM_BL(7:12,:,:),1,'omitnan');

MedianMat_Int2=[CalcSum_Wake_BLmed_Int2;CalcSum_NREM_BLmed_Int2;CalcSum_REM_BLmed_Int2];

%ZT13-18
CalcSum_Wake_BLmed_Int3=median(CalcSum_Wake_BL(12:18,:,:),1,'omitnan');
CalcSum_NREM_BLmed_Int3=median(CalcSum_NREM_BL(12:18,:,:),1,'omitnan');
CalcSum_REM_BLmed_Int3=median(CalcSum_REM_BL(12:18,:,:),1,'omitnan');

MedianMat_Int3=[CalcSum_Wake_BLmed_Int3;CalcSum_NREM_BLmed_Int3;CalcSum_REM_BLmed_Int3];

%ZT19-24
CalcSum_Wake_BLmed_Int4=median(CalcSum_Wake_BL(19:24,:,:),1,'omitnan');
CalcSum_NREM_BLmed_Int4=median(CalcSum_NREM_BL(19:24,:,:),1,'omitnan');
CalcSum_REM_BLmed_Int4=median(CalcSum_REM_BL(19:24,:,:),1,'omitnan');

MedianMat_Int4=[CalcSum_Wake_BLmed_Int4;CalcSum_NREM_BLmed_Int4;CalcSum_REM_BLmed_Int4];


%WAKE-ON, make a matrix with Wake, NREM and REM median activity for Wake ON
WakeON_Int1=NaN(size(MedianMat,1),size(MedianMat,2), numanim);
NREMON_Int1=NaN(size(MedianMat,1),size(MedianMat,2), numanim);
REMON_Int1=NaN(size(MedianMat,1),size(MedianMat,2), numanim);

WakeON_Int2=WakeON_Int1;
WakeON_Int3=WakeON_Int1;
WakeON_Int4=WakeON_Int1;

NREMON_Int2=NREMON_Int1;
NREMON_Int3=NREMON_Int1;
NREMON_Int4=NREMON_Int1;

REMON_Int2=REMON_Int1;
REMON_Int3=REMON_Int1;
REMON_Int4=REMON_Int1;

clear i

for i=1:numanim
indxWpref=find(MaxMat(:,:,i)==1);
WakeON_Int1(:,indxWpref,i)=MedianMat_Int1(:,indxWpref,i);
WakeON_Int2(:,indxWpref,i)=MedianMat_Int2(:,indxWpref,i);
WakeON_Int3(:,indxWpref,i)=MedianMat_Int3(:,indxWpref,i);
WakeON_Int4(:,indxWpref,i)=MedianMat_Int4(:,indxWpref,i);

indxNpref=find(MaxMat(:,:,i)==2);
NREMON_Int1(:,indxNpref,i)=MedianMat_Int1(:,indxNpref,i);
NREMON_Int2(:,indxNpref,i)=MedianMat_Int2(:,indxNpref,i);
NREMON_Int3(:,indxNpref,i)=MedianMat_Int3(:,indxNpref,i);
NREMON_Int4(:,indxNpref,i)=MedianMat_Int4(:,indxNpref,i);

indxRpref=find(MaxMat(:,:,i)==3);
REMON_Int1(:,indxRpref,i)=MedianMat_Int1(:,indxRpref,i);
REMON_Int2(:,indxRpref,i)=MedianMat_Int2(:,indxRpref,i);
REMON_Int3(:,indxRpref,i)=MedianMat_Int3(:,indxRpref,i);
REMON_Int4(:,indxRpref,i)=MedianMat_Int4(:,indxRpref,i);
end

% plot how each ROI changes activity during Wake NREM REM depending if they
% are WakeON NREMON or REMON
A=(size(WakeON_Int1,2))*size(WakeON_Int1,3);
WakeON_Int1_rs=rmmissing(reshape(WakeON_Int1,[3,A]),2);
NREMON_Int1_rs=rmmissing(reshape(NREMON_Int1,[3,A]),2);
REMON_Int1_rs=rmmissing(reshape(REMON_Int1,[3,A]),2);
clear A
A=(size(WakeON_Int2,2))*size(WakeON_Int2,3);
WakeON_Int2_rs=rmmissing(reshape(WakeON_Int2,[3,A]),2);
NREMON_Int2_rs=rmmissing(reshape(NREMON_Int2,[3,A]),2);
REMON_Int2_rs=rmmissing(reshape(REMON_Int2,[3,A]),2);

clear A
A=(size(WakeON_Int3,2))*size(WakeON_Int3,3);
WakeON_Int3_rs=rmmissing(reshape(WakeON_Int3,[3,A]),2);
NREMON_Int3_rs=rmmissing(reshape(NREMON_Int3,[3,A]),2);
REMON_Int3_rs=rmmissing(reshape(REMON_Int3,[3,A]),2);

clear A
A=(size(WakeON_Int4,2))*size(WakeON_Int4,3);
WakeON_Int4_rs=rmmissing(reshape(WakeON_Int4,[3,A]),2);
NREMON_Int4_rs=rmmissing(reshape(NREMON_Int4,[3,A]),2);
REMON_Int4_rs=rmmissing(reshape(REMON_Int4,[3,A]),2);



%% Split the SD data in 6h intervals based on the median activity
% to compare to BL later on
%ZT0-6
CalcSum_Wake_SDmed_Int1=median(CalcSum_Wake_SD(1:6,:,:),1,'omitnan');
CalcSum_NREM_SDmed_Int1=median(CalcSum_NREM_SD(1:6,:,:),1,'omitnan');
CalcSum_REM_SDmed_Int1=median(CalcSum_REM_SD(1:6,:,:),1,'omitnan');

MedianMatSD_Int1=[CalcSum_Wake_SDmed_Int1;CalcSum_NREM_SDmed_Int1;CalcSum_REM_SDmed_Int1];

%ZT7-12
CalcSum_Wake_SDmed_Int2=median(CalcSum_Wake_SD(7:12,:,:),1,'omitnan');
CalcSum_NREM_SDmed_Int2=median(CalcSum_NREM_SD(7:12,:,:),1,'omitnan');
CalcSum_REM_SDmed_Int2=median(CalcSum_REM_SD(7:12,:,:),1,'omitnan');

MedianMatSD_Int2=[CalcSum_Wake_SDmed_Int2;CalcSum_NREM_SDmed_Int2;CalcSum_REM_SDmed_Int2];

%ZT13-18
CalcSum_Wake_SDmed_Int3=median(CalcSum_Wake_SD(12:18,:,:),1,'omitnan');
CalcSum_NREM_SDmed_Int3=median(CalcSum_NREM_SD(12:18,:,:),1,'omitnan');
CalcSum_REM_SDmed_Int3=median(CalcSum_REM_SD(12:18,:,:),1,'omitnan');

MedianMatSD_Int3=[CalcSum_Wake_SDmed_Int3;CalcSum_NREM_SDmed_Int3;CalcSum_REM_SDmed_Int3];

%ZT19-24
CalcSum_Wake_SDmed_Int4=median(CalcSum_Wake_SD(19:24,:,:),1,'omitnan');
CalcSum_NREM_SDmed_Int4=median(CalcSum_NREM_SD(19:24,:,:),1,'omitnan');
CalcSum_REM_SDmed_Int4=median(CalcSum_REM_SD(19:24,:,:),1,'omitnan');

MedianMatSD_Int4=[CalcSum_Wake_SDmed_Int4;CalcSum_NREM_SDmed_Int4;CalcSum_REM_SDmed_Int4];


%make a matrix with wake, NREM and REM median activity for Wake ON, NREMON
%and REMON for the SD file
WakeONSD_Int1=NaN(size(MedianMatSD_Int1,1),size(MedianMatSD_Int1,2), numanim);
NRONSD_Int1=NaN(size(MedianMatSD_Int1,1),size(MedianMatSD_Int1,2), numanim);
RONSD_Int1=NaN(size(MedianMatSD_Int1,1),size(MedianMatSD_Int1,2), numanim);

WakeONSD_Int2=WakeONSD_Int1;
WakeONSD_Int3=WakeONSD_Int1;
WakeONSD_Int4=WakeONSD_Int1;

NRONSD_Int2=NRONSD_Int1;
NRONSD_Int3=NRONSD_Int1;
NRONSD_Int4=NRONSD_Int1;

RONSD_Int2=RONSD_Int1;
RONSD_Int3=RONSD_Int1;
RONSD_Int4=RONSD_Int1;

clear i indxWpref indxNpref indxRpref

for i=1:numanim
indxWpref=find(MaxMat(:,:,i)==1);
WakeONSD_Int1(:,indxWpref,i)=MedianMatSD_Int1(:,indxWpref,i);
WakeONSD_Int2(:,indxWpref,i)=MedianMatSD_Int2(:,indxWpref,i);
WakeONSD_Int3(:,indxWpref,i)=MedianMatSD_Int3(:,indxWpref,i);
WakeONSD_Int4(:,indxWpref,i)=MedianMatSD_Int4(:,indxWpref,i);

indxNpref=find(MaxMat(:,:,i)==2);
NRONSD_Int1(:,indxNpref,i)=MedianMatSD_Int1(:,indxNpref,i);
NRONSD_Int2(:,indxNpref,i)=MedianMatSD_Int2(:,indxNpref,i);
NRONSD_Int3(:,indxNpref,i)=MedianMatSD_Int3(:,indxNpref,i);
NRONSD_Int4(:,indxNpref,i)=MedianMatSD_Int4(:,indxNpref,i);

indxRpref=find(MaxMat(:,:,i)==3);
RONSD_Int1(:,indxRpref,i)=MedianMatSD_Int1(:,indxRpref,i);
RONSD_Int2(:,indxRpref,i)=MedianMatSD_Int2(:,indxRpref,i);
RONSD_Int3(:,indxRpref,i)=MedianMatSD_Int3(:,indxRpref,i);
RONSD_Int4(:,indxRpref,i)=MedianMatSD_Int4(:,indxRpref,i);
end

% plot how each ROI changes activity during Wake NREM REM depending if they
% are WakeONSD NRONSD or RONSD
A=(size(WakeONSD_Int1,2))*size(WakeONSD_Int1,3);
WakeONSD_Int1_rs=rmmissing(reshape(WakeONSD_Int1,[3,A]),2);
NRONSD_Int1_rs=rmmissing(reshape(NRONSD_Int1,[3,A]),2);
RONSD_Int1_rs=rmmissing(reshape(RONSD_Int1,[3,A]),2);
clear A
A=(size(WakeONSD_Int2,2))*size(WakeONSD_Int2,3);
WakeONSD_Int2_rs=rmmissing(reshape(WakeONSD_Int2,[3,A]),2);
NRONSD_Int2_rs=rmmissing(reshape(NRONSD_Int2,[3,A]),2);
RONSD_Int2_rs=rmmissing(reshape(RONSD_Int2,[3,A]),2);

clear A
A=(size(WakeONSD_Int3,2))*size(WakeONSD_Int3,3);
WakeONSD_Int3_rs=rmmissing(reshape(WakeONSD_Int3,[3,A]),2);
NRONSD_Int3_rs=rmmissing(reshape(NRONSD_Int3,[3,A]),2);
RONSD_Int3_rs=rmmissing(reshape(RONSD_Int3,[3,A]),2);

clear A
A=(size(WakeONSD_Int4,2))*size(WakeONSD_Int4,3);
WakeONSD_Int4_rs=rmmissing(reshape(WakeONSD_Int4,[3,A]),2);
NRONSD_Int4_rs=rmmissing(reshape(NRONSD_Int4,[3,A]),2);
RONSD_Int4_rs=rmmissing(reshape(RONSD_Int4,[3,A]),2);


%% plot the data from Int2  ZT6-12 after SD 
% scatter plot with median, (Signal is not continous) a
% remove outliers first, all values > 2 SD are removed
WakeON_Int2_rs_rm=rmoutliers(WakeON_Int2_rs');
NREMON_Int2_rs_rm=rmoutliers(NREMON_Int2_rs');
REMON_Int2_rs_rm=rmoutliers(REMON_Int2_rs');

WakeONSD_Int2_rs_rm=rmoutliers(WakeONSD_Int2_rs');
NREMONSD_Int2_rs_rm=rmoutliers(NRONSD_Int2_rs');
REMONSD_Int2_rs_rm=rmoutliers(RONSD_Int2_rs');

%Anderson-Darling test for Normality:
clear i
for i=1:size(WakeON_Int2_rs_rm,2);
    clear p1 p2 p3 h
[h,p1] = adtest(WakeON_Int2_rs_rm(:,i));
[h,p2] = adtest(NREMON_Int2_rs_rm(:,i));
[h,p3] = adtest(REMON_Int2_rs_rm(:,i));
pWBL(:,i)=p1
pNBL(:,i)=p2
pRBL(:,i)=p3
end

clear i
for i=1:size(WakeONSD_Int2_rs_rm,2);
    clear p1 p2 p3 h
[h,p1] = adtest(WakeONSD_Int2_rs_rm(:,i));
[h,p2] = adtest(NREMONSD_Int2_rs_rm(:,i));
[h,p3] = adtest(REMONSD_Int2_rs_rm(:,i));
pWSD(:,i)=p1
pNSD(:,i)=p2
pRSD(:,i)=p3
end


clear x x1 x2 xa x1a x2a
x(:,1)=0.5 + (0.25+0.25)*rand((size(WakeON_Int2_rs_rm,1)),1);
x(:,2)=2   + (0.25+0.25)*rand((size(WakeON_Int2_rs_rm,1)),1);
x(:,3)=3.5 + (0.25+0.25)*rand((size(WakeON_Int2_rs_rm,1)),1);

x1(:,1)=0.5 + (0.25+0.25)*rand((size(NREMON_Int2_rs_rm,1)),1);
x1(:,2)=2   + (0.25+0.25)*rand((size(NREMON_Int2_rs_rm,1)),1);
x1(:,3)=3.5 + (0.25+0.25)*rand((size(NREMON_Int2_rs_rm,1)),1);

x2(:,1)=0.5 + (0.25+0.25)*rand((size(REMON_Int2_rs_rm,1)),1);
x2(:,2)=2   + (0.25+0.25)*rand((size(REMON_Int2_rs_rm,1)),1);
x2(:,3)=3.5 + (0.25+0.25)*rand((size(REMON_Int2_rs_rm,1)),1);

xa(:,1)=1   + (0.25+0.25)*rand((size(WakeONSD_Int2_rs_rm,1)),1);
xa(:,2)=2.5 + (0.25+0.25)*rand((size(WakeONSD_Int2_rs_rm,1)),1);
xa(:,3)=4   + (0.25+0.25)*rand((size(WakeONSD_Int2_rs_rm,1)),1);

x1a(:,1)=1   + (0.25+0.25)*rand((size(NREMONSD_Int2_rs_rm,1)),1);
x1a(:,2)=2.5 + (0.25+0.25)*rand((size(NREMONSD_Int2_rs_rm,1)),1);
x1a(:,3)=4   + (0.25+0.25)*rand((size(NREMONSD_Int2_rs_rm,1)),1);

x2a(:,1)=1   + (0.25+0.25)*rand((size(REMONSD_Int2_rs_rm,1)),1);
x2a(:,2)=2.5 + (0.25+0.25)*rand((size(REMONSD_Int2_rs_rm,1)),1);
x2a(:,3)=4 + (0.25+0.25)*rand((size(REMONSD_Int2_rs_rm,1)),1);

sz=60; %size of the circles

% median
medW_I2=median(WakeON_Int2_rs_rm);
medN_I2=median(NREMON_Int2_rs_rm);
medR_I2=median(REMON_Int2_rs_rm);

medWSD_I2=median(WakeONSD_Int2_rs_rm);
medNSD_I2=median(NREMONSD_Int2_rs_rm);
medRSD_I2=median(REMONSD_Int2_rs_rm);

F5a=figure % Figure 2, Manuscripts
subplot(131)
s1=scatter(x,WakeON_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#808080'); hold on
l1=line([0.5,1],[medW_I2(1),medW_I2(1)],'Color','k','LineWidth',2); hold on
l2=line([2,2.5],[medW_I2(2),medW_I2(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4],[medW_I2(3),medW_I2(3)],'Color','k','LineWidth',2); hold on
s2=scatter(xa,WakeONSD_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#A2142F'); hold on
l1a=line([1,1.5],[medWSD_I2(1),medWSD_I2(1)],'Color','r','LineWidth',2); hold on
l2a=line([2.5,3],[medWSD_I2(2),medWSD_I2(2)],'Color','r','LineWidth',2); hold on
l3a=line([4,4.5],[medWSD_I2(3),medWSD_I2(3)],'Color','r','LineWidth',2); hold off
ylim([-0.3 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
legend([s1(1) s2(1)],'BL','SD')
legend('boxoff')
title('Wake')
ylabel('median Ca^{2+} activity ZT6-12')
subplot(132)
scatter(x1,NREMON_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#808080')
l1=line([0.5,1],[medN_I2(1),medN_I2(1)],'Color','k','LineWidth',2); hold on
l2=line([2,2.5],[medN_I2(2),medN_I2(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4],[medN_I2(3),medN_I2(3)],'Color','k','LineWidth',2); hold on
scatter(x1a,NREMONSD_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#A2142F')
l1a=line([1,1.5],[medNSD_I2(1),medNSD_I2(1)],'Color','r','LineWidth',2); hold on
l2a=line([2.5,3],[medNSD_I2(2),medNSD_I2(2)],'Color','r','LineWidth',2); hold on
l3a=line([4,4.5],[medNSD_I2(3),medNSD_I2(3)],'Color','r','LineWidth',2); hold on
ylim([-0.3 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
title('NREM')
subplot(133)
scatter(x2,REMON_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#808080')
l1=line([0.5,1],[medR_I2(1),medR_I2(1)],'Color','k','LineWidth',2); hold on
l2=line([2,2.5],[medR_I2(2),medR_I2(2)],'Color','k','LineWidth',2); hold on
l3=line([3.5,4],[medR_I2(3),medR_I2(3)],'Color','k','LineWidth',2); hold on
scatter(x2a,REMONSD_Int2_rs_rm,sz,'filled','MarkerFaceColor','none','MarkerEdgeColor','#A2142F')
l1a=line([1,1.5],[medRSD_I2(1),medRSD_I2(1)],'Color','r','LineWidth',2); hold on
l2a=line([2.5,3],[medRSD_I2(2),medRSD_I2(2)],'Color','r','LineWidth',2); hold on
l3a=line([4,4.5],[medRSD_I2(3),medRSD_I2(3)],'Color','r','LineWidth',2); hold on
ylim([-0.3 0.6])
xlim([0 5])
set(gca,'XTick',[1 2.5 4],'XTickLabel',{'Wake','NREM','REM'})
title('REM')

%%% STATS, run tests for Int2 for the SD vs. BL data, some groups are not
%%% parametric, therefor use a Kruskall-Wallis, followed by a Tukey-Kramer
% also needs a grouping vector.
BLmat=[{'BLW'},{'BLN'},{'BLR'}];
SDmat=[{'SDW'},{'SDN'},{'SDR'}];

% for Wake-ON ROIs
condW=[(repmat(BLmat,size(WakeON_Int2_rs_rm,1),1));(repmat(SDmat,size(WakeONSD_Int2_rs_rm,1),1))];
condW=condW(:);
Data_Wake=[WakeON_Int2_rs_rm;WakeONSD_Int2_rs_rm];
Data_Wake=Data_Wake(:);

[pw,hw,statsw]=kruskalwallis(Data_Wake,condW);
pw1=multcompare(statsw,'CriticalValueType','tukey-kramer')

% for NREM-ON ROIs
condN=[(repmat(BLmat,size(NREMON_Int2_rs_rm,1),1));(repmat(SDmat,size(NREMONSD_Int2_rs_rm,1),1))];
condN=condN(:);
Data_NREM=[NREMON_Int2_rs_rm;NREMONSD_Int2_rs_rm];
Data_NREM=Data_NREM(:);

[pn,hn,statsn]=kruskalwallis(Data_NREM,condN);
pn1=multcompare(statsn,'CriticalValueType','tukey-kramer')

% for REM-ON ROIs
condR=[(repmat(BLmat,size(REMON_Int2_rs_rm,1),1));(repmat(SDmat,size(REMONSD_Int2_rs_rm,1),1))];
condR=condR(:);
Data_REM=[REMON_Int2_rs_rm;REMONSD_Int2_rs_rm];
Data_REM=Data_REM(:);

[pr,hr,statsr]=kruskalwallis(Data_REM,condR);
pr1=multcompare(statsr,'CriticalValueType','tukey-kramer')
