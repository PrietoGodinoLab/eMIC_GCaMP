
%?F/F data from 'wdw' frames flanking each peak location is extracted from 
%each segment, and the location of the maximum peak identified, together 
%with its width at half-height. For each wave, the slope of a linear fit 
%for peak locations as a function of VNC segment is used to determine its 
%type, with a negative slope identifying a forward wave and a positive slope 
%identifying a backward wave. Peak widths were used to determine the starting 
%and ending points of a wave.

clear;
clearvars; 
close all;

load('deltaFbase_3Darray.mat')
load('peakLocs_conv.mat')

Fs = 0.47; % sampling frequency in Hz, i.e. volume rate
wdw = 3; % # of frames to take left and right of peak Loc to do cross correllations; 3 is a good one
start = 1; % from which segment (1-3 is T1-T3, 4-11 is A1-A8/9)
stop = 11; % to which segment (11 is A8/9)

segm = size(deltaCods,1);
timePts = size(deltaCods,2);
larvae = size(deltaCods,3);



% cell arrays to put stuff in
wavePeaks = cell(larvae,1);
waveStarts = cell(larvae,1);
waveEnds = cell(larvae,1);

% get chunks of size wdw*2+1 around each peak location and place traces in arrays
% get peak location and width at half-height for each chunk
for k = 1:larvae % sample index; 1-10 are CTRL, 11-20 are MUT and 21-30 are RES
    for j = start:stop % segment index; 1-11 are Left T1-A8/9 and 12-22 are Right T1-A8/9


        L = deltaCods(j,:,k);
        R = deltaCods(j+segm/2,:,k);
        
        avgLR = (L+R)/2;

        for i = 1:size(peakLocs{k,1},2)

            loc =  peakLocs{k,1}(1,i);
            
            %to account for locs at the edges of the signal
            if loc <= wdw
               loc = wdw + 1;
            elseif loc+wdw >= timePts
               loc = timePts - wdw - 1;
            else loc = loc;
            end

            chunk = avgLR(loc-wdw:loc+wdw);
            chunk = chunk + abs(min(chunk)); % bring baseline to 0
            
            [~, wPeakLocs, width] = findpeaks(chunk);
                
                if size(wPeakLocs,2) > 1
                   [~,ind] = min(abs(wPeakLocs-wdw+1)); % find the peak closest to the middle
                   wPeakLocs = wPeakLocs(1,ind);
                   width = width(1,ind);
                   wStart = wPeakLocs - width/2;
                   wEnd = wPeakLocs + width/2;
                   
                elseif size(wPeakLocs,2) < 1
                   wPeakLocs = NaN;
                   width = NaN;
                   wStart = NaN;
                   wEnd = NaN;
                else
                    wStart = wPeakLocs - width/2;
                    wEnd = wPeakLocs + width/2;
                end
             
             wavePeaks{k,1}(j-start+1,i) = wPeakLocs;                            
             waveStarts{k,1}(j-start+1,i) = wStart;
             waveEnds{k,1}(j-start+1,i) = wEnd;
             
        end
                      
    end
    
end



% find forward and backward waves
waveType = cell(larvae,1);
slope = cell(larvae,1);
fwSlope = cell(larvae,1);
bwSlope = cell(larvae,1);

for k = 1:larvae % sample index; 1-10 are CTRL, 11-20 are MUT and 21-30 are RES
    
    for i = 1:size(peakLocs{k,1},2)

        wave = wavePeaks{k,1}(:,i);
        
        % remove indeces of NaNs; X used for line fit     
        X = [1:stop-start+1]' .* ~isnan(wave); 
        X = X(find(X));    
        
        % squeeze out NaNs
        wave =  wave(~isnan(wave)); 
        
        
        if size(unique(wave),1) == 1 % if the peak is in the same location in all segments
           waveType{k,1}(i,1) = 2; % index = 2 for synchronous activity
           slope{k,1}(i,1) = 0;
           fwSlope{k,1}(i,1) = NaN;
           bwSlope{k,1}(i,1) = NaN;
           
        else
            lineFit = polyfit(X,wave,1);
            slope{k,1}(i,1) = lineFit(1,1);
            
            if lineFit(1,1) < 0
                fwSlope{k,1}(i,1) = lineFit(1,1);
                bwSlope{k,1}(i,1) = NaN;
            else
                fwSlope{k,1}(i,1) = NaN;
                bwSlope{k,1}(i,1) = lineFit(1,1);
            end
            % forward wave if slope < 0, backward wave if slope > 0
            waveType{k,1}(i,1) = lineFit(1,1) < 0; 
                % true for forward, false for backward
        end
        
    end
end


for k = 1 : larvae % sample index; 1-10 are CTRL, 11-20 are MUT and 21-30 are RES
          
    wTypes = waveType{k,1};

    wBackward (k,1) = size(find(wTypes == 0),1);
    wForward(k,1) = size(find(wTypes == 1),1);
    wSynchr(k,1) = size(find(wTypes == 2),1);
    
    sFw = fwSlope{k,1};
    sFw = sFw(~isnan(sFw));
    avgFwSlope(k,1) = mean(sFw);
    
    sBw = bwSlope{k,1};
    sBw = sBw(~isnan(sBw));
    avgBwSlope(k,1) = mean(sBw);
            
end

% store the sum of backward and forward waves
wSum = wForward +  wBackward ;


% find wave durations
waveDuration = cell(larvae,1);
for k = 1:larvae
    
    s = waveStarts{k,1};
    e = waveEnds{k,1};    
    dur = max(e,[],1) - min(s,[],1);
    waveDuration{k,1} = dur;
    
    clear s e dur
end

avgWaveDuration = cellfun(@mean,waveDuration);
avgWaveDuration = avgWaveDuration/Fs; % express values in seconds


% find backward wave durations
waveDurationBckw = zeros(larvae,1);
for k = 1:larvae
    
    t = waveType{k,1};
    d = waveDuration{k,1}';
    
    idx = find(t==0);   
    
    waveDurationBckw(k,1) = mean(d(idx));
    
    clear t d idx
end

waveDurationBckw = waveDurationBckw/Fs; % express values in seconds


% find forward wave durations
waveDurationFw = zeros(larvae,1);
for k = 1:larvae
    
    t = waveType{k,1};
    d = waveDuration{k,1}';
    
    idx = find(t==1);   
    
    waveDurationFw(k,1) = mean(d(idx));
    
    clear t d idx
end
 waveDurationFw =  waveDurationFw/Fs; % express values in seconds
