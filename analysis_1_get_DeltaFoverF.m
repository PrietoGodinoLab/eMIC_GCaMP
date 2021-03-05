%This script imports calcium imaging data from .csv files with each row 
%representing one ROI and columns representing frames. 
%Next, it calculates ?F/F = (F(t) - F0) / F0, where F(t) is the 
%fluorescence intensity at a given time point and F0 is the mean 
%fluorescence intensity in a manually defined window lacking any spontaneous 
%activity, spanning 10 frames (here, 21.3 s). 


clc;
clear;
clearvars;
close all;

% load start frames for subtracting the baseline
% here, a 3 x 1 cell, 1 for each genotype
load('baseline_indeces.mat');

% User-defined variables
genotypes = 3;
movDuration = 1300; % number of frames to take

%% Import data

% Get file names and paths

fileNames = cell(genotypes,1);  
filePaths = cell(genotypes,1);  


for i = 1:genotypes
    [fileNames{i,1}, filePaths{i,1}] =  uigetfile('*.csv','MultiSelect','on');
end


%% Import cods and get deltaF/F

cods = cell(genotypes,1);
baseCods = cell(genotypes,1);


for i = 1:genotypes
        for k = 1:length(fileNames{i,1})
            importcods = readmatrix([filePaths{i,1},'/',fileNames{i,1}{1,k}])';
            importcods = importcods(2:23,1:movDuration); %gets rid of the first row which is just indeces and selects first 1300 frames
            cods{i,1}{k,1} = importcods;
              
            baseIdx = baseline{i,1}(k,1);
            base = mean(importcods(:,baseIdx:baseIdx+10),2);
            baseCods{i,1}{k,1} = base;
            if i == 1          
                deltaCods{k,1} = (cods{i,1}{k,1} - baseCods{i,1}{k,1}(k,1)) / baseCods{i,1}{k,1}(k,1);
                
            elseif i == 2
                idx = length(fileNames{1,1});
                deltaCods{idx+k,1} = (cods{i,1}{k,1} - baseCods{i,1}{k,1}(k,1)) / baseCods{i,1}{k,1}(k,1);
            elseif i == 3
                idx = length(fileNames{1,1})+length(fileNames{2,1});
                deltaCods{idx+k,1} = (cods{i,1}{k,1} - baseCods{i,1}{k,1}(k,1)) / baseCods{i,1}{k,1}(k,1);
            end
           
           
        end
end
%% save data as 3D array

deltaCods = cell2mat(deltaCods');
deltaCods = reshape(deltaCods,22,movDuration,30);

