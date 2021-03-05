%Use this script to bring baseline to zero and convolute traces to correct
%for bleaching


clc;
clear;
clearvars;
close all;

% load file
load('deltaFbase_3Darray.mat');


ROI = size(deltaCods,1);
timePts = size(deltaCods,2);
larvae = size(deltaCods,3);

% create array to store convoluted mean traces in 
convCods = zeros(larvae, timePts);

%% Mean the cods

meanCods = zeros(larvae, timePts);

for i = 1:larvae
    meanCods(i, :) = mean(deltaCods(:,:,i));
end

%% Convolute to correct baseline

% do it manually for each one
i = 1;

% remove the baseline
% comment this if you don't need it
x = linspace(0,1,1300);
base = 0.7*sin(0.2*pi*x); 

y = base;
%y = conv(meanCods(i,:),base);
toPlot = meanCods(i,:) + y(1:1300);


% bring baseline to zero
if min(toPlot) < 0
    toPlot = toPlot + abs(min(toPlot));
elseif min(toPlot) > 0
    toPlot = toPlot - min(toPlot);
end

% check if it worked
figure, 
plot(meanCods(i,:), 'k','LineWidth',2)
hold on
plot(toPlot, 'r','LineWidth',1)
plot(y,'b','LineWidth',2 )
xlim([0 size(meanCods(i,:),2)])
grid on

% store convoluted mean trace
convCods(i,:) = toPlot;

%% Generates a stacked plot of convoluted traces for each genotype

ctrl = convCods(1:10,:);
mut = convCods(11:20,:);
res = convCods(21:end,:);

% choose which data to plot
data = res;

numGrps = size(data,1);
yVal = 2; % add this to plot 
adders = (0:yVal:(numGrps-1)*yVal)';

 
    dataToPlot =data + flip(adders); % adders flipped to have T1 at the top and A8-9 at the bottom

    % plot it
    figure,
    plot(dataToPlot', 'LineWidth',1.5, 'color', 'k')
    hold on
    set(gca,'LineWidth',2)
    set(gca,'YTickLabels',[])
    set(gca,'YTick',[])
    set(gca,'YColor','none')
    ylabel('Fluorescence intensity', 'color', 'k')
    xlabel('Time (min)')
    xticks([0:300:size(dataToPlot,2)]) % for a tick every 10 min
    xticklabels ({})
    set(gcf,'InvertHardCopy','Off')
    set(gcf,'color','w')
    box off


%% SAVE

save('deltaFbase_3Darray_means_conv.mat','convCods')