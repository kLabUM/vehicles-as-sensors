%% This script is needed to generate fused precipitation maps for different days
clear all;close all;clc

% select date of the storm events on 2014: 6/12, 6/28, 8/11

strdate = '0811';
% strdate = '0612';
% strdate = '0628';

dSamp = 12;         % down-sample rate, e.g., larger the value less accurate
weight = 1;         % weights applied to the radar prior, e.g., 0.1, 0.5, 0.9

strnc = strcat('../data/data_2014',strdate,'.nc');

ncdisp(strnc);
radarData = ncread(strnc,'radar');
lonNet = ncread(strnc,'longitude');
latNet = ncread(strnc,'latitude');

% set radar data to zero if it's NaNs

radarData(isnan(radarData)) = 1e-9;

strcsv = strcat('../data/camera_combined_filtered_2014',strdate,'_1min.csv');
vehicleData = csvread(strcsv,1,1);

updatedmap = mapUpdate(weight,dSamp,vehicleData,radarData,lonNet,latNet);

for i = 1:length(updatedmap)
    drawFigure(updatedmap{i,1},radarData(:,:,i),dSamp);     % radar
    drawFigure(updatedmap{i,2},radarData(:,:,i),dSamp);     % product
    hold on;
    if ~isempty(updatedmap{i,3})
        plot(updatedmap{i,3}(:,1),updatedmap{i,3}(:,2),'ok')    % vehicle locations
        hold on;
    end
end


