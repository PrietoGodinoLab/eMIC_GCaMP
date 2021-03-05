%To identify the number of waves and their mean amplitudes, 
% this code averages ?F/F data across segments uses "findpeaks" to extract 
%the locations and amplitudes of local maxima in each sample. 

clc;
clear;
clearvars;
close all;

% load file & Mean the cods
load('deltaFbase_3Darray_means_conv.mat');

ROI = size(convCods,1);
timePts = size(convCods,2);
larvae = size(convCods,3);


% variables to store peak info in
peakAmps = cell(size(convCods,1),1);
peakLocs = cell(size(convCods,1),1);
peakWidths = cell(size(convCods,1),1);
peakProms = cell(size(convCods,1),1);
minPeakH = zeros(size(convCods,1),1);
minPeakProm = zeros(size(convCods,1),1);
minPeakW = zeros(size(convCods,1),1);

%% GET PEAK DATA

% do it manually for each one
i = 1; % % sample index; 1-10 are CTRL, 11-20 are MUT and 21-30 are RES

% find best values for each sample
peakH = 0.7;
peakP = 0;
peakW = 0;

toPlot = convCods(i,1:1300);
[peakAmps{i}, peakLocs{i}, peakWidths{i}, peakProms{i}] = findpeaks(toPlot,'MinPeakHeight', peakH, 'MinPeakProminence', peakP, 'MinPeakWidth', peakW);

minPeakH(i) = peakH;
minPeakProm(i) = peakP;
minPeakW(i) = peakW;


figure,
hold on
grid on
findpeaks(toPlot, 'MinPeakHeight', peakH, 'MinPeakProminence', peakP, 'MinPeakWidth', peakW); 
xticks([0:31.5789:size(toPlot,2)])
xticklabels ({})
set(gcf,'InvertHardCopy','Off')
set(gcf,'color','w')
hold off


disp (['# of peaks found = ', num2str(size(peakAmps{i},2))]) 

%% get means of peak amplitudes, widths & prominences

for i = 1:size(convCods,1)
    meanAmps(i,1) = mean(peakAmps{i});
    meanWidths(i,1) = mean(peakWidths{i});
    meanProms(i,1) = mean(peakProms{i});
    maxWIdths(i,1) = max(peakWidths{i});
end


%% SAVE

save('peakAmps_conv.mat','peakAmps')
save('peakLocs_conv.mat','peakLocs')
save('peakProms_conv.mat','peakProms')
save('peakWidths_conv.mat','peakWidths')
save('minPeakH_conv.mat','minPeakH')
save('minPeakProm_conv.mat','minPeakProm')
save('minPeakW_conv.mat','minPeakW')