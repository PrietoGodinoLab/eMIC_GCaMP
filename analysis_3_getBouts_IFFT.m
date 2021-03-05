
%This script identifies activity bouts on ?F/F data averaged across segments
%The data is smoothened using a low-pass filter, and local maxima identified 
%using the "findpeaks" function. Smoothened traces are obtained using a 
%Short-time Fourier Transform (FT), followed by an inverse FT with a 
%low-pass filter. 

%The data is split into fragments of 100 frames (~ 3.5 min) overlapping by ½,
%each fragment is windowed using a Hamming function, its discrete FT computed, 
%and the resulting magnitudes and phases are used to reconstitute the data 
%from the lowest 5 frequency components. 

% Windowing is necessary! One assumption of FFT is that you can fit an
% integer amount of cycles in your signal, i.e. that you can stitch the end
% to the beginning of the signal and get perfect periodicity, but that's
% never the case. Windowing reduces the values at the edges of the signal
% close to zero so that they can be fitted together without making a big
% 'jump' 

clear;
clearvars; 
close all;

% load data
load('deltaFbase_3Darray_means_conv.mat')
load('peakLocs_conv.mat')

% user-input variables
Fs = 0.47; % sampling frequency in Hz, i.e. volume rate
fragLength = 100; % the data is split into fragments of this length 
freqWin = 5;% size of frequencies window to include for bouts; 5 works

% real frequencies in Hz for FFT
f = Fs*(0:fragLength/2)/fragLength; 

% get data array dimensions
larvae = size(convCods,1);
timePts = size(convCods,2);

% fragments of data to analyse
fragStarts = 1:fragLength/2:timePts-fragLength; % so that fragments overlap by 1/2
nFrags = length(fragStarts); % number of fragments
outputFrags = fragStarts+(fragLength/4); % needed to reconstitute data 'y' after FT
outputLength = fragLength/2; % needed to reconstitute data 'y' after FT

% Hamming window
wndw = hamming(fragLength+1); % uses a Hamming windowing function; the window length
                              % is an odd number

% arrays for holding frequency magnitudes and amplitudes, bouts
Mags = zeros((fragLength/2)+1,nFrags);
Amps = zeros((fragLength/2)+1,nFrags);
allWavesOnly = zeros(larvae,timePts-fragLength);
allAmps = zeros(larvae, (fragLength/2)+1);
allBoutsOnly = zeros(larvae, timePts);
allBoutLocs = cell(larvae,1);
allBoutsNum = zeros(larvae,1);

for k = 1:larvae
   
    data = convCods(k,:)';

    % array for storing the IFFT output signal
    y = zeros(timePts,1)+mean(data);


    % loopy
    for i = 1:nFrags


        % get out a fragment and zero-mean
        fragData = (data(fragStarts(i):fragStarts(i)+fragLength));
        fragMean = mean(fragData);
        fragData = fragData-fragMean;

        % get the frequency information
        X = fft(fragData.*wndw);

        % get the magnitude and the phase separately
        XMag = real(X);
        XPhs = imag(X);

        % get magnitudes from half the D=discrete Fourier Transform, 
        % i.e. the positive frequency values (and ignore the negative ones)
        Mags(:,i) = abs(X(1:(fragLength/2)+1));

        % convert magnitudes into amplitudes
        amp = abs(X(1:fragLength/2+1)/fragLength);
        amp(2:end-1) = 2*amp(2:end-1); 
        Amps(:,i) = amp;


        % set some magnitudes and phase to zero to use for IFFT
        XMagLow = XMag;
        XPhsLow = XPhs;
        zeroIdx = [1+freqWin,(fragLength+1)-(freqWin-1)];
        XMagLow(zeroIdx(1):zeroIdx(2)) = 0;
        XPhsLow(zeroIdx(1):zeroIdx(2)) = 0;

        % IFFT
        dataRecLow = ifft(complex(XMagLow,XPhsLow));

        % divide by the Hamming window
        dataRecLow = dataRecLow./wndw; 

        % add the mean back
        dataRecLow = dataRecLow+fragMean;

        % take a chunk
        y(outputFrags(i):outputFrags(i)+(outputLength-1)) = ...
            dataRecLow((fragLength/4)+1:((fragLength/4)+1)+(outputLength-1));

    end
    
    % get location and number of bouts
    [~, bLocs] = findpeaks(y, 'MinPeakHeight', max(y)/2);  
    allBoutLocs{k,1} = bLocs;
    allBoutsNum(k,:) = size(bLocs,1);
    
    % store IFFT of mean trace (only shows bouts, waves are filtered out)
    allBoutsOnly(k,:) = y;
    
    % store amps
    allAmps(k,:) = mean(Amps,2);


end

