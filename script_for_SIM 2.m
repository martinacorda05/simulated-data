
% This code generates 3 brain sources with specific phase syncrhonization for
% a certain time interval, which can be used as ground-truth to test
% time-verying connectivity approaches 

% Is is based on the work of Gordon et al. 

% Gordon S M, Franaszczuk P J, Hairston W D, Vindiola M and McDowell K 2013 
% Comparing parametric and nonparametric methods for detecting 
% phase synchronization in EEG J Neurosci Methods 212 247–58Gordon S M, 
% Franaszczuk P J, Hairston W D, Vindiola M and McDowell K 2013

%%
clc
clear all
close all

%% Parameters
% Parameters can be varied depending on the dataset properties you want to simulate 

fs = 250;  % Sampling frequency
time=1; % length of the simulated signal in seconds
t = 0:1/fs:time-1/fs; % time vector  
num_trials = 100; % number of trials 
frequencies = linspace(1, 50, 50); % vector of frequencies from 1Hz to 50Hz 
num_sources = 3; % number of sources generated
num_samples = length(t); %number of samples 
SNR=10; % signal to noise ratio
n_rep=50; % number of time to repeat the signal generation
alpha = [8, 13];  % frequency range of the syncrhonization (arbitrarly defined)
start_and_length = [0.3 0.2]; % vector of the time point representing the onset of the syncrhonization and the duration (arbitrarly defined)
% if other syncrhonization intervals are defined, you should add a row
brain_sources_50_noise={}; % cell initialization for the results 

%%
for rep=1:n_rep

% variables initialization 
brain_sources = zeros(num_sources, num_samples, num_trials);
noise= zeros(num_sources, num_samples, num_trials);

for trial = 1:num_trials
    
    s1 = zeros(1, num_samples); % source number 1 
    s2 = zeros(1, num_samples); % source number 2 
    s3 = zeros(1, num_samples); % source number 3
    
    for f = frequencies
        amplitude=1/f; %1/f amplitude effect
        phase1 = 2*pi*rand;
        phase2 = 2*pi*rand;

        s1 = s1 + amplitude * sin(2 * pi * f * t + phase1);
        if f >= alpha(1) && f <= alpha(2) 
            for period = 1:size(start_and_length, 1)
                phase_sync=phase1+pi/4; % imposed synchronization (arbitrarly defined)
                start_sync_idx = round(fs * start_and_length(period, 1)); % index where syncrhonization starts
                end_sync_idx = round(fs * (start_and_length(period, 1) + start_and_length(period, 2))); % index where syncronization ends 
                s2(start_sync_idx:end_sync_idx) = s2(start_sync_idx:end_sync_idx) + amplitude * sin(2 * pi * f * t(start_sync_idx:end_sync_idx) + phase_sync);
            end

            s2(1:start_sync_idx-1) = s2(1:start_sync_idx-1) + amplitude * sin(2 * pi * f * t(1:start_sync_idx-1) + phase2); 
            s2(end_sync_idx+1:end) = s2(end_sync_idx+1:end) + amplitude * sin(2 * pi * f * t(end_sync_idx+1:end) + phase2);
        
        else
            s2 = s2 + amplitude * sin(2 * pi * f * t + phase2); 
        end
    end
    
    % s3 is not syncrhonized 
    for f = frequencies
        amplitude=1/f; 
        phase = 2*pi*rand;
        s3 = s3 + amplitude * sin(2 * pi * f * t + phase);
    end
    
    % add random noise 
    wn=randn(1,size(s1,2)); 
    noise1=(norm(s1)/ (norm(wn) * sqrt (SNR))) * wn; 

    wn=randn(1,size(s2,2)); 
    noise2=(norm(s2)/ (norm(wn) * sqrt (SNR))) * wn; 

    wn=randn(1,size(s3,2)); 
    noise3=(norm(s3)/ (norm(wn) * sqrt (SNR))) * wn;

    noise(1,:,trial)=noise1; 
    noise(2,:,trial)=noise2; 
    noise(3,:,trial)=noise3; 


    brain_sources(1, :, trial) = s1;
    brain_sources(2, :, trial) = s2;
    brain_sources(3, :, trial) = s3; 
end

%%
brain_sources_new=brain_sources+noise; 
% reshape of the brain sources matrix (number of trials, num sources,
% number of samples) 

Y=zeros(num_trials,num_sources,num_samples); 
for i=1:num_trials 
    Y(i,:,:)=brain_sources_new(:,:,i); 
end 

brain_sources_50_noise{rep}=Y; % values are saved in a cell
end 

%% save the results 
% save('brain_sources_100','brain_sources'); 


