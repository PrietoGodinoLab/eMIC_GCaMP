
%This script measures absolute differences in intensity between each 
%left-right hemi-segment pair. Absolute difference is taken into account if
%it is greater than a predefined threshold of 0.7 ?F/F per sample. 
%The threshold was established by comparing visually observed one-sided 
%activity to noise in two example samples.

clc;
clear;
clearvars; 
close all;

load('deltaFbase_3Darray.mat')

% threshold for counting a frame as asynchronous
thresh = 0.7; % selected after looking at MUT2 T3 vs. CTRL3 T3

% how many continuous frames above threshold constitute one event
frPerEvent = 1;

segm = size(deltaCods,1); %j
timePts = size(deltaCods,2);
larvae = size(deltaCods,3); %k

% array for holding baseline-subtracted cods
deltaCodsSubtrBase = zeros(segm,timePts,larvae);

for k = 1:larvae
    for j = 1:segm
        
      data = deltaCods(j,:,k);
      
        if min(data) < 0
            data = data + abs(min(data));
        elseif min(data) > 0
            data = data - min(data);
        end
      
        deltaCodsSubtrBase(j,:,k) = data;
    end
end

% array for holding absolute differences between Left and Right hemi-segments
absDiff = zeros(segm/2,timePts,larvae);
absDiffmean = zeros(segm/2,larvae);
for k = 1:larvae
    for j = 1:segm/2-1
        
        Lsegm = deltaCods(j,:,k);
        Rsegm = deltaCods(j+segm/2,:,k);
        
        LRdiff = abs(Lsegm-Rsegm);
        
        absDiff(j,:,k) = LRdiff;
        absDiffmean (j,k) = mean(LRdiff);
    end
end

absDiffThresh = zeros(segm/2,timePts,larvae);
absDiffThresh = absDiff >= thresh;

% get the average and sum of timepoints in which the absolute difference is
% above threshold for each sample
S = squeeze(sum(absDiffThresh,2));
avgS = mean(S,1);
sumS = sum(S,1);

% array for holding the number of events with differences in intensity 
%between Left and Right hemi-segments
absDiffThreshEvents = zeros(segm/2,larvae);
maxfrPerEvent = zeros(segm/2,larvae);

for k = 1:larvae
    for j = 1:segm/2-1
        
       connected = absDiffThresh(j,:,k);
       connectedObj = bwconncomp(connected);
       
       numConnFrames = zeros(size(connectedObj.PixelIdxList,2),1);
       
       for c = 1:size(connectedObj.PixelIdxList,2)
           
            numConnFrames(c,:) = size(connectedObj.PixelIdxList{1,c},1);
       end
    
        if size(numConnFrames,1) > 0
            maxfrPerEvent(j,k) = max(numConnFrames);
        else maxfrPerEvent(j,k) = 0;
        end
            
    events = sum(numConnFrames >= frPerEvent);
    absDiffThreshEvents(j,k) = events;  

    end
end


% get the average and sum of events in which the absolute difference is
% above threshold for each sample
avgEvents = mean(absDiffThreshEvents,1);
sumEvents =  sum(absDiffThreshEvents,1);


