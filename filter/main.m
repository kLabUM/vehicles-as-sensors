% Spatial distribution estimation (SPE) using Bayesian Approach (MMLE + SIR filter)
% Hyongju Park & Matt Bartos
clear all;close all;clc
%% 

nsampInfo = 100;     % number of information sample
dSamp = 12;          % down-sample, e.g., 1/8, smaller the value, the more accurate

varInfo = 0.001;     % how noisy the sensor is (e.g., 0: perfect) default 0.01....
varPos = 0.01;       % decay as the distance between the windshield wiper measurement and the source of rain, increases...
weight = 0.9;		 % weight on prior (e.g., if weight=1, it forgets about the prior measurements and only consider the new radar measurement, default: 0.9)

vehicleData = csvread('../data/20140811_cameracombined_filtered.csv', 1, 1);
deviceIdIdx = 2 - 1;
tripIx = 3 - 1;      
latIx = 4 - 1;       
lonIx = 5 - 1;       
wiperIx = 6 - 1;     
gpsSpeedIx = 7 - 1;
cameraIx = 8 - 1;
radarIx = 9 - 1;
gageIx = 10 - 1; 
tStepIx = 11 - 1;
yStepIx = 12 - 1;
xStepIx = 13 - 1;

allVehicleID = unique(vehicleData(:,deviceIdIdx))';
numExcludedVehicles = 0;
excludedVehicleID = allVehicleID(1,randsample(1:length(allVehicleID),numExcludedVehicles));

ncdisp('../data/data_20140811.nc');
gageData = ncread('../data/data_20140811.nc','gage');
radarData = ncread('../data/data_20140811.nc','radar');
wiperData = ncread('../data/data_20140811.nc','wiper');

%set radar data to zero if it's NaN
radarData(isnan(radarData)) = 1e-9;
wiperData(isnan(wiperData)) = 0;
lonNet = ncread('../data/data_20140811.nc','longitude');
latNet = ncread('../data/data_20140811.nc','latitude');

%% load data
% change scale (GPS locations -> [0, 1]x[0, 1]) % for the sake of convenience
lonNetScaled=(lonNet-min(lonNet))/(max(lonNet) - min(lonNet));
latNetScaled=(latNet-min(latNet))/(max(latNet) - min(latNet));

%radius of rain detection for each vehicle (scaled based upon our scaled area [0, 1]x[0, 1])
radD = 0.1;

% generate timeseries data {radar, windshield wiper}
for i = 1:size(radarData,3)
    radarTSeries{i} = radarData(:,:,i);
    wiperTSeries{i} = wiperData(:,:,i);    
end

% extract effective radar measurements (with non-NaNs)
radarIdx = [];
for i = 1:length(radarTSeries)
    if ~isempty(find(radarTSeries{i}, 1))
        radarIdx(end+1) = i;
    end
end
for i = 1:length(radarIdx)
    radar_nz(i) = length(find(radarTSeries{radarIdx(i)}));
end
for i = 1:length(radar_nz)
    radar_nz2(i) = radarIdx(i);
end

% number of runs for a finite time horizeon
nRuns = max(unique(vehicleData(:,end-2)));

% generate matrix M to plot the radar measurements
M(:,1) = repmat(lonNetScaled,size(radarTSeries{1},2),1);
lvec = [];
for i= 1:size(radarTSeries{1},2)
    lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{1},1),1)];
end
M(:,2) = lvec;
% normalize
M(:,3) = radarTSeries{1}(:);
[qx,qy,qz] = drawNoFigure(M,radarTSeries{1},dSamp);

%%
% clear M for later use...
clear M;
M = [qx(:) qy(:) qz(:)];
% let ground truth be 1 for all regions...
gTruth = ones(size(M,1),1);

% generate samples for information vector
hSet2 = haltonset(1,'Skip',1e3,'Leap',1e2);
hScrambled2 = scramble(hSet2,'RR2');

maxes = []

for i=1:nRuns
    maxes(i) = max(max(radarTSeries{i}));    
end

grandmax = max(maxes);

% find maximum information...
sampInfo = net(hScrambled2,nsampInfo) ;

% apply different weights to priors:
%   - prvWgt0: prior from radar measurements
%   - prvWgt1: uniform prior (no information)
% prvWgt0 = [];
% for i = 1:size(M,1)
%     p0 = mvnpdf(sampInfo,M(i,3),varInfo)/sum(mvnpdf(sampInfo,M(i,3),varInfo)); 
%     prvWgt0 = [prvWgt0 p0];
% end
prvWgt1 = ones(size(sampInfo,1),size(M,1))/size(sampInfo,1);
prvWgt = prvWgt1;

% generate samples for locations (use M)
sampPos = M(:,1:2);

% expected value for the information, e.g., rain intensity...
particleWgt = sampInfo'*prvWgt;

% update M using the weighted prior
M(:,3) =particleWgt'/max(sampInfo);

%%
% ++++++++++++PLOT++++++++++++
[qx,qy,qz] = drawFigure(M,radarTSeries{1},dSamp);
% ++++++++++++PLOT++++++++++++

% save current data
savData{1,1} = sampPos;
savData{1,2} = particleWgt;
nRuns = max(unique(vehicleData(:,end-2)));

%% create netcdf
% 
nccreate('../data/product_20140811.nc','lat', 'Dimensions', {'lat', length(qx(1,:))});
ncwrite('product_20140811.nc','lat', linspace(min(latNet), max(latNet), length(qx(1,:))));
 
nccreate('../data/product_20140811.nc','lon', 'Dimensions', {'lon', length(qy(:,1))});
ncwrite('product_20140811.nc','lon', linspace(min(lonNet), max(lonNet), length(qy(:,1))));
 
nccreate('../data/product_20140811.nc','normmax');
ncwrite('product_20140811.nc','normmax', grandmax);
 
nccreate('../data/product_20140811.nc','combined',...
         'Dimensions', {'time', inf, 'lat', length(qx(1,:)), 'lon', length(qy(:,1))});

%% run filter

k5 = 0;
for curStep = 150:nRuns
    curStep
    clear Mt;
    Mt(:,1) = repmat(lonNetScaled,size(radarTSeries{curStep},2),1);
    lvec = [];
    for i= 1:size(radarTSeries{curStep},2)
        lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{curStep},1),1)];
    end
    Mt(:,2) = lvec;
    % normalize
    Mt(:,3) = radarTSeries{curStep}(:)/grandmax;
    [qx,qy,qz] = drawNoFigure(Mt,radarTSeries{curStep},dSamp);
    clear Mt;
    Mt = [qx(:) qy(:) qz(:)];    
    
	% if there is a non-zero, non-NaN radar measurement, we will update the previous belief...
    if ismember(curStep,radar_nz2)
        clear M;
        % DO SENSOR FUSION...
        % generate matrix M to plot the radar measurements
        M(:,1) = repmat(lonNetScaled,size(radarTSeries{curStep},2),1);
        lvec = [];
        for i= 1:size(radarTSeries{curStep},2)
            lvec = [lvec;repmat(latNetScaled(i),size(radarTSeries{curStep},1),1)];
        end
        M(:,2) = lvec;
        % normalize
        M(:,3) = radarTSeries{curStep}(:)/grandmax;
        [qx,qy,qz] = drawFigure(M,radarTSeries{curStep},dSamp);
        filestr = sprintf('./img_05/r%03d.png', curStep);
        saveas(gcf, filestr);        
        clear M;
        M = [qx(:) qy(:) qz(:)];
        prvWgt0 = [];
        for i = 1:size(M,1)
            p0 = mvnpdf(sampInfo,M(i,3),varInfo)/sum(mvnpdf(sampInfo,M(i,3),varInfo)); 
            prvWgt0 = [prvWgt0 p0];
        end
		% current belief = weight * (PDF from radar measurement) * (1 - weight) * prior belief
        prvWgt = weight*prvWgt0 + (1-weight)*prvWgt;
        for i = 1:size(sampPos,1)
            particleWgtRad{curStep}(:,i) = sampInfo' * prvWgt0(:,i);
        end        
    else
        for i = 1:size(sampPos,1)
            tmp = rand(size(sampInfo));
            particleWgtRad{curStep}(:,i) = sampInfo' * tmp/sum(tmp);
        end
    end
    
	% find time frame where there is windshield-wiper measurement 
    wiperOnIdx = find(vehicleData(:,tStepIx) == curStep);
    if ~isempty(wiperOnIdx)
		% find unique vehicle IDs
        vehicleID = unique(vehicleData(:,deviceIdIdx));
        wiperInfo = [];        
        vehicleDataEff = vehicleData(wiperOnIdx,:);
		% scale GPS locations to fit [0, 1] x [0, 1]
        vehicleDataEff(:,lonIx)=(vehicleDataEff(:,lonIx)-min(lonNet))/(max(lonNet) - min(lonNet));
        vehicleDataEff(:,latIx)=(vehicleDataEff(:,latIx)-min(latNet))/(max(latNet) - min(latNet));        
		% find the closest points from the downsampled map to vehicle locations
        clear nearPos;
        for i = 1:size(vehicleDataEff,1)  
            [~,nearPos(i)] = nearestPntDist([vehicleDataEff(i,lonIx) vehicleDataEff(i,latIx)],[qx(:),qy(:)]);
        end

		% [winshield wiper data, vehicle positions on the map]
        vehicleDataEff = [vehicleDataEff nearPos'];
		
		% find unique vehicle positions
        [nearPosUniq, idv] = unique(nearPos);
        clear nearPosMed;
        k = 0;
		k1 = 0;
        for i = 1:length(nearPosUniq)
			% find windshield wiper value(5), and vehicleID(1) at unique vehicle positions...
            nearPosUniqWiper{i,1}=vehicleDataEff(find(vehicleDataEff(:,end) == nearPosUniq(i)),wiperIx);
            nearPosUniqWiper{i,2}=vehicleDataEff(find(vehicleDataEff(:,end) == nearPosUniq(i)),deviceIdIdx);

			% for each vehicle (with different IDs)
            for j = 1:length(vehicleID)
                % if the vehicle is not excluded from the experiment...
                if ~ismember(vehicleID(j),excludedVehicleID)
                    % take median of the windshield-wiper measurements made by each vehicle
                    if ~isempty(find(nearPosUniqWiper{i,2}==vehicleID(j)))
                        k = k+1;
                        t_val = median(nearPosUniqWiper{i,1}(find(nearPosUniqWiper{i,2}==vehicleID(j)),:));
                        if isnan(t_val)
                            t_val = 0;
                        end
                        nearPosMed(k,:) = [nearPosUniq(i) vehicleID(j) min(t_val,2)];
                        
                    end
                else
                    if ~isempty(find(nearPosUniqWiper{i,2}==vehicleID(j)))
                        invalid = 0;
                        k1 = k1+1;
                        t_val = median(nearPosUniqWiper{i,1}(find(nearPosUniqWiper{i,2}==vehicleID(j)),:));
                        if isnan(t_val)
                            t_val = 0;
                        end
                        
                        nearPosMedEx(k1,:) = [nearPosUniq(i) vehicleID(j) min(t_val,2)];
                        
                    else
                        invalid = 1;
                    end                    
                end
            end
        end
        
        clear M;

        M(:,1) = qx(:);M(:,2) = qy(:);
		
		% find indice for unique vehicle positions with windshield-wiper measurements
        if exist('nearPosMed')
            rainIdx = nearPosMed(:,1);
        end

		maxRainIntensity = 1;		% average rain intensity that is captured by windshield wiper 	(level 2)
		% obtain rain detection likelihoods from the wind-shield wiper measurements
        for i_2 = 1:size(M,1)
            misDetectLhd = 1;
            if exist('nearPosMed')
                [out1, out2] = nearestPntDist(M(i_2,1:2),M(rainIdx,1:2));
                if out1 <= radD
                        if ((nearPosMed(out2,3) > 0.1) & (nearPosMed(out2,3) < 0.5))
                            if radarTSeries{curStep}(i_2) >= 0.025
                                gTruth(i_2) = Mt(i_2,3);
                            else
                                gTruth(i_2) = 0.01;
                            end
                        elseif ((nearPosMed(out2,3) >= 0.5) & (nearPosMed(out2,3) <= 1.5))
                            if radarTSeries{curStep}(i_2) >= 0.025
                                gTruth(i_2) = Mt(i_2,3);
                            else
                                gTruth(i_2) = 0.025;
                            end
                        elseif ((nearPosMed(out2,3) >= 1.5) & (nearPosMed(out2,3) <= 2.5))
                            if radarTSeries{curStep}(i_2) >= 0.025
                                gTruth(i_2) = Mt(i_2,3);
                            else
                                gTruth(i_2) = 0.05;
                            end
                        elseif ((nearPosMed(out2,3) >= 2.5) & (nearPosMed(out2,3) <= 3))
                            if radarTSeries{curStep}(i_2) >= 0.025
                                gTruth(i_2) = Mt(i_2,3);
                            else
                                gTruth(i_2) = 0.05;
                            end
                        elseif nearPosMed(out2,3) < 0.1
                            gTruth(i_2) = 0.0;
                        else
                            gTruth(i_2) = Mt(i_2,3);
                        end
                    misDetectLhd = misDetectLhd * (1-mvnpdf(M(i_2,1:2),M(rainIdx(out2),1:2),eye(2)*varPos)/mvnpdf([0 0],[0 0],eye(2)*varPos));
                else
                    misDetectLhd = 1;
                    gTruth(i_2)=Mt(i_2,3);
                end   
            end
            detectLhd(i_2) = 1-misDetectLhd;
        end

		% obtain observation likelihoods
        for i = 1:size(sampPos,1)
            for j = 1:length(sampInfo)
                  infoLhd{i}(j)=normpdf(sampInfo(j),gTruth(i,:),varInfo)/normpdf(0,0,varInfo);
            end

            for j = 1:length(sampInfo)
                  infoLhdU{i}(j)=normpdf(sampInfo(j),Mt(i,3),varInfo)/normpdf(0,0,varInfo);
            end
            infoLhd{i} = infoLhd{i} * detectLhd(i)+ infoLhdU{i} * (1-detectLhd(i));
            % dirac-delta function approximation
            if all(infoLhd{i} == 0)
                [~,idx_tmp] = min((repmat(gTruth(i,:),size(sampInfo,1),1) - sampInfo).^2);
                infoLhd{i}(idx_tmp) = 1;
            end
        end
		% particle filtering
        for i = 1:size(sampPos,1)
            wgt2(:,i) = sirFilter(infoLhd{i}',prvWgt(:,i));
            prvWgt(:,i) = wgt2(:,i);                
        end
        for i = 1:size(sampPos,1)
            particleWgt(:,i) = sampInfo' * wgt2(:,i);
        end
        clear M;        
        M(:,1) = qx(:);
        M(:,2) = qy(:);        
        M(:,3) =particleWgt'./max(sampInfo);    
        set(gca,'color','none');
        [qx,qy,qz] = drawFigure(M,radarTSeries{curStep},dSamp);
        hold on;
        plot(M(nearPosUniq,1),M(nearPosUniq,2),'cd');
        A = reshape(permute(qz, [2 1]), [1 size(permute(qz, [2 1]))]);
        ncwrite('../data/product_20140811.nc','combined', A, [curStep 1 1]);
        clf;
        close;
    else
        clear M;
        M(:,1) = qx(:);
        M(:,2) = qy(:);        
        M(:,3) = Mt(:,3);
        set(gca,'color','none');
        [qx,qy,qz] = drawFigure(M,radarTSeries{curStep},dSamp);
        hold on;
        plot(M(nearPosUniq,1),M(nearPosUniq,2),'x');
        A = reshape(permute(qz, [2 1]), [1 size(permute(qz, [2 1]))]);
        ncwrite('../data/product_20140811.nc','combined', A, [curStep 1 1]);
        clf;
        close;
        invalid = 1;
    end   