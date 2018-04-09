%%
% This script generate ROC curves to evaluate the classification
% performance of the windshield wiper and the wiper plus radar (product) 

clear all;close all;clc

%% settings
nsampInfo = 10000;    % number of information sample
dSamp = 12;          % down-sample, e.g., 1/8, smaller the value, the more accurate

varInfo = 0.0001;     % how noisy the sensor is (e.g., 0: perfect) default 0.01....
varPos = 0.0001;     % decay as the distance between the windshield wiper measurement and the source of rain, increases...

vehicleDataU = csvread('./data/camera_combined_aditya.csv',1,0);
vehicleDataU(:,1) = round(vehicleDataU(:,1));
% this eliminates duplicate rows
[~,vidt,vcdt]= unique(vehicleDataU(:,1:6),'rows');
vehicleData = vehicleDataU(vidt,:);

% column numbers for table
tStepIx = 1;
deviceIdIdx = 2 ;
tripIx = 3 ; 
latIx = 4 ; 
lonIx = 5 ; 
wiperIx = 6 ;
cameraIx = 7;

% load netcdf data
ncdisp('/Users/mdbartos/Data/combined_nc/data_20140811.nc');
gageData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','gage');
radarData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','radar');
wiperData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','wiper');
timeData = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','time');

% covert timeData to milliseconds (begining from 12:00:00:000 am)
for i = 1:length(timeData)
    ns = uint64(timeData(i));
    wholeSecs = floor(double(ns)/1e9);
    fracSecs = double(ns - uint64(wholeSecs)*1e9)/1e9;
    timeVal(i,:) = datevec(datetime(wholeSecs,'ConvertFrom','posixTime','Format','HH:mm:ss.SSSSSSSSS') + seconds(fracSecs));
    timeinSec(i) = timeVal(i,6)+timeVal(i,5)*60+60*60*timeVal(i,4);
end
%set radar data to zero if it's NaN
radarData(isnan(radarData)) = 0;
wiperData(isnan(wiperData)) = 0;
lonNet = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','longitude');
latNet = ncread('/Users/mdbartos/Data/combined_nc/data_20140811.nc','latitude');

% GPS locations to positions in [0,1]x[0,1] 
vehicleData(:,5)=(vehicleData(:,5)-min(lonNet))./(max(lonNet) - min(lonNet));
vehicleData(:,4)=(vehicleData(:,4)-min(latNet))./(max(latNet) - min(latNet));

% generate timeseries data {radar, windshield wiper}
for i = 1:size(radarData,3)
    radarTSeries{i} = radarData(:,:,i);
%     radarTSeries{i} = imresize(radarTSeriesU{i}, 0.05, 'bicubic');
    wiperTSeries{i} = wiperData(:,:,i);    
end
lon_size = size(radarTSeries{i},1);     % number of grids to represent longitude axis, e.g., 800
lat_size = size(radarTSeries{i},2);     % number of grids to represent longitude axis, e.g., 299

figure,
spy(radarTSeries{200}')

% generate information samples
hSet2 = haltonset(1,'Skip',1e3,'Leap',1e2);
hScrambled2 = scramble(hSet2,'RR2');
sampInfo = net(hScrambled2,nsampInfo) ;

maxes = [];
for i=1:length(radarTSeries)
    maxes(i) = max(max(radarTSeries{i}));    
end
% find the maximum radar intensity value (for normalization purpose)
grandmax = max(maxes);

% random point generation to find the circumcenter of the positions of at
% least 2 vehicles
hSet3 = haltonset(2,'Skip',1e3,'Leap',1e2);
hScrambled3 = scramble(hSet3,'RR2');
pos = net(hScrambled3,10000);

% find all unique time steps (made at every 2-3 seconds using Aditya's
% data)
timeHorz = unique(vehicleData(:,tStepIx));

for u1 = 1:length(timeHorz)     % for each time from the 'timeHorz'
    curStep2 = timeHorz(u1);      % current time step

    % interpolate radar data to cover the time horizon
    for b1 = 1:length(timeinSec)
        if timeinSec(b1) <curStep2
            curStep = b1;
        end
    end
    % find time frame where there is windshield-wiper measurement 
    wiperOnIdx = find(vehicleData(:,tStepIx) == curStep2);
    % add to the last column of 'vehicleData', time synchronized radar
    % data
    for y1 =1:size(wiperOnIdx,1)
        vehicleData(wiperOnIdx(y1),8)=radarTSeries{curStep}(max(1,round((vehicleData(wiperOnIdx(y1),4)*size(radarTSeries{1},1)))),max(1,round((vehicleData(wiperOnIdx(y1),5))*size(radarTSeries{1},2))))/grandmax;
    end
    % number of vehicles that appear at each time step
    nCar(u1) = length(wiperOnIdx);

    % each time step --> radar time index (between 1-212)
    nRad(u1) = curStep;

end

%% roc curve (wiper vs radar against camera) generated with aditya's data
% this is non-time aggregated--non-space aggregated version

[f01,t01,~,auc01] = perfcurve(vehicleData(:,cameraIx),vehicleData(:,wiperIx),1);
[f02,t02,~,auc02] = perfcurve(vehicleData(:,cameraIx),vehicleData(:,8),1);
f0 = figure('position',[100 100  800 800],'Color',[1 1 1]);
p1 = plot(f01, t01,'k-','LineWidth',2);hold on  % product-WW 
p2 = plot(f02, t02,'k--','LineWidth',2);hold on % radar-WW 

legend([p1,p2],{strcat('wiper, AUC=',sprintf('%.3f',auc01)),strcat('radar, AUC=',sprintf('%.3f',auc02))},'location','southeast')

p4 = line([0 1],[0 1],'LineStyle','-.','Color','k');
p4.Annotation.LegendInformation.IconDisplayStyle = 'off';

ylabel('True Positive Rate (Sensitivity)');
xlabel('False Positive Rate (1-Specificity)');
% title('Product vs Radar: Time&Space-aggregated  (Tian)');

axis('equal');
axis([0 1 0 1]);
set(gca,'FontSize',16);

%% 
% for quick evaluation, for the experiment, we only consider instances
% where at least two vehicles become sufficienly near to each other:

% first find time instances where there multiple devices are in the vicinity
[~,idx] = sort(nCar);
timeinstances = idx(:,end-100:end);

% set the threshold for the maximum distance between two vehicles

rad = 0.030864 * 1.39;   % the value must be small enough; however if it is
% two small, it is nearly impossible to find such vehicle rendezvous event
% 0.001 approximately translates to 0.0324 miles = 0.0519 km = 50.19 m
% 1 miles = 0.030864, thus if rad = 0.030864 * 1.45 = 1.45 mile


% brief description of the following for--loop
% for each time step,
% using the pseudo random sample points, find the location of the center of
% the smallest circle with radius 'rad' which contains at least 2
% vehicles, and record the measurement information from those vehicles
for i = 1:length(timeinstances)
    addtolist1{i} = [];
    a1 = vehicleData((vehicleData(:,tStepIx) == timeHorz(timeinstances(i))),4:5);
    tmp1 = vehicleData((vehicleData(:,tStepIx) == timeHorz(timeinstances(i))),:);
    for b1 = 1:size(pos,1)
        distcoll = [];
        for b2 = 1:size(a1,1)
            distcoll = [distcoll norm(pos(b1,:)-a1(b2,:))]; 
        end
        rec(b1) = length(find(distcoll <= rad));
        rec1{b1} = find(distcoll <= rad);
    end    
    [numavlcar,maxidxrec] = max(rec);
    samppntmx(i,:) = pos(maxidxrec,:);
    radTime =nRad(timeinstances(i));
    
    % if in each circle there is at least 2 vehicles, save
    %       1) information of the device/trip within the circle
    %       2) 
    if numavlcar >=2
        addtolist1{i} = tmp1(rec1{maxidxrec},:);
        radtolist1{i} = radTime;
        cntrod{i} = samppntmx(i,:);
    end
end

%% table creation 
% createTblFull = [radar | product | camera (gt) | device_id ] (nx4)
%% table creation 
% createTblFull = [radar | product | camera (gt) | device_id ] (nx4)
cnt = 0;
createTblFull = [];
% for each time 
for i = 1:length(addtolist1)
    % if there is a point (centered at the point, you can draw a circle 
    % with radius 'rad' that encircle at least two circle
    if ~isempty(addtolist1{i})
        cnt = cnt + 1;
        % skip the empty entries
        addtolist2{cnt}=addtolist1{i};
        cntrod2{cnt} = cntrod{i};
        
        % read the radar value at the point, at the given time
        radtolist2(cnt)=radarTSeries{radtolist1{i}}(round(cntrod{i}(1)*size(radarTSeries{1},1)),round(cntrod{i}(2)*size(radarTSeries{1},2)))/grandmax;
        % read 'wiper | camera' measurement
        addtolist3 = addtolist2{cnt}(:,6:7);
        % read the position of the centroid
        cnter = cntrod{i};
        % beyasian filter
        % generates product using radar and wiper measurement,
        % where '0.1' was as an estimate of rain intensity if wiper reading
        % is 'on' and '0' otherwise
        for j = 1:size(addtolist3,1)
            % detection likelihood
            mdp(j) = (mvnpdf(addtolist2{cnt}(j,4:5),cnter,eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
            % precipitation estimate = binary(wiper measurements) * 0.1
            % (likelihood of rain)
            prepestimate = 0.01*addtolist3(j,1); 
            for m = 1:length(sampInfo)
                  infoLhd(m)=normpdf(sampInfo(m),prepestimate,varInfo)/normpdf(0,0,varInfo);
            end
            for m = 1:length(sampInfo)
                  infoLhdU(m)=normpdf(sampInfo(m),radtolist2(cnt),varInfo)/normpdf(0,0,varInfo);
            end
            addtolist3(j,1) = (infoLhd * mdp(j)+ infoLhdU * (1-mdp(j)))*sampInfo;
        end
        createTbl{cnt} = [ones(size(addtolist2{cnt},1),1)*radtolist2(cnt) addtolist3 addtolist2{cnt}(:,2)];
        createTblFull = [createTblFull;ones(size(addtolist2{cnt},1),1)*radtolist2(cnt) addtolist3 addtolist2{cnt}(:,2)];
    end
end

%% without cross-validation
% createTblFull(:,3): camera
% createTblFull(:,2): product
% createTblFull(:,1): radar

% generate fprs, tprs, aucs
[fpr1,tpr1,~,auc1] = perfcurve(createTblFull(:,3),createTblFull(:,2),1);
[fpr2,tpr2,~,auc2] = perfcurve(createTblFull(:,3),createTblFull(:,1),1);

% generate roc curves
% f0 = figure('position',[100 100  800 800],'Color',[1 1 1]);
% p1 = plot(fpr1, tpr1,'k-','LineWidth',2);hold on  % product-WW 
% p2 = plot(fpr2, tpr2,'k--','LineWidth',2);hold on % radar-WW 
% 
% legend([p1,p2],{strcat('product, AUC=',sprintf('%.3f',auc1)),strcat('radar, AUC=',sprintf('%.3f',auc2))},'location','southeast')
% 
% p4 = line([0 1],[0 1],'LineStyle','-.','Color','k');
% p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% ylabel('True Positive Rate (Sensitivity)');
% xlabel('False Positive Rate (1-Specificity)');
% 
% axis('equal');
% axis([0 1 0 1]);
% set(gca,'FontSize',16);

%% generate results using leave-one-out-cross-validation
% LOOCV (I): for the experiment, leave 1 vehicle out, and genearte product
% using other vehicles, and compare it against the ground truth (camera
% data) generated from the left-out vehicle
% (I): this is a time-independent version of LOOCV

% find all the device id's in the scene
allvehiclesID = unique(createTblFull(:,4));
createTblEAug= [];

% try exclude one device at a time
for i = 1:length(allvehiclesID)
    % 
    excludedVecID = allvehiclesID(i);
    for j = 1:length(createTbl)
        flag = 0;
        for k = 1:size(createTbl{j},1)
            if createTbl{j}(k,4) == excludedVecID
                createTblE{i}{j} = createTbl{j}(setdiff(1:size(createTbl{j},1),k),:);
                % use ground truth from the left-out vehicle to evaluate
                % all other data
                createTblE{i}{j}(:,3) = ones(size(createTblE{i}{j}(:,3),1),1)*createTbl{j}(k,3);
                flag = 1;
            end
        end
        if flag == 0
            createTblE{i}{j} = createTbl{j};
        end
        createTblEAug = [createTblEAug;createTblE{i}{j}];
    end
end
[fpr1,tpr1,~,auc1] = perfcurve((createTblEAug(:,3)),createTblEAug(:,2),1);
[fpr2,tpr2,~,auc2] = perfcurve((createTblEAug(:,3)),createTblEAug(:,1),1);

% f0 = figure('position',[100 100 800 800],'Color',[1 1 1]);
% p1 = plot(fpr1, tpr1,'k-','LineWidth',2);hold on  % product-WW 
% p2 = plot(fpr2, tpr2,'k--','LineWidth',2);hold on % radar-WW 
% 
% legend([p1,p2],{strcat('product, AUC=',sprintf('%.3f',auc1)),strcat('radar, AUC=',sprintf('%.3f',auc2))},'location','southeast')
% 
% p4 = line([0 1],[0 1],'LineStyle','-.','Color','k');
% p4.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% ylabel('True Positive Rate (Sensitivity)');
% xlabel('False Positive Rate (1-Specificity)');
% 
% axis('equal');
% axis([0 1 0 1]);
% set(gca,'FontSize',16);

%%
% LOOCV (II)
% at each time step, leave 1 vehicle out, and genearte product
% using other vehicles, and compare it against the ground truth
% (II) this version is time-dependent

compdata = [];
for i = 1:length(createTbl)
    for j = 1:size(createTbl{i},1)
        compdata = [compdata;[mean(createTbl{i}(setdiff(1:size(createTbl{i},1),j),1:2),1) createTbl{i}(j,3)]];
    end
end

[fpr1,tpr1,~,auc1] = perfcurve(compdata(:,3),compdata(:,2),1);
[fpr2,tpr2,~,auc2] = perfcurve(compdata(:,3),compdata(:,1),1);

f0 = figure('position',[100 100  800 800],'Color',[1 1 1]);
p1 = plot(fpr1, tpr1,'k-','LineWidth',2);hold on  % product-WW 
p2 = plot(fpr2, tpr2,'k--','LineWidth',2);hold on % radar-WW 

legend([p1,p2],{strcat('product, AUC=',sprintf('%.3f',auc1)),strcat('radar, AUC=',sprintf('%.3f',auc2))},'location','southeast')

p4 = line([0 1],[0 1],'LineStyle','-.','Color','k');
p4.Annotation.LegendInformation.IconDisplayStyle = 'off';

ylabel('True Positive Rate (Sensitivity)');
xlabel('False Positive Rate (1-Specificity)');
% title('Product vs Radar: Time&Space-aggregated  (Tian)');

axis('equal');
axis([0 1 0 1]);
set(gca,'FontSize',16);

%% display locations where vehicles meet (at common times)
f0 = figure('position',[100 100  800 299],'Color',[1 1 1]);
for i = 1:length(addtolist2)
    plot(addtolist2{i}(1,4),addtolist2{i}(1,5),'o','MarkerSize',6);hold on;
    plot(addtolist2{i}(2,4),addtolist2{i}(2,5),'o','MarkerSize',6);hold on;
    plot(cntrod2{i}(:,1),cntrod2{i}(:,2),'s','MarkerSize',6);hold on;
    line([addtolist2{i}(1,4) cntrod2{i}(:,1)],[addtolist2{i}(1,5) cntrod2{i}(:,2)],'Color','k');hold on;
    line([cntrod2{i}(:,1) addtolist2{i}(2,4)],[cntrod2{i}(:,2) addtolist2{i}(2,5)],'Color','k');hold on;
end
axis([0 1 0 1]);