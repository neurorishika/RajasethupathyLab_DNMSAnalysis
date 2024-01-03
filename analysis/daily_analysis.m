% MATLAB CODE
% By Rishika Mohanta, The Rockefeller University, 2023

% CLEAR WORKSPACE
clear all; close all; clc;

% SET IMPORTANT VARIABLES
datafolder = '../data/delay_match_to_sample/round5_additional_training/';
filtername = '_DNMS_gonogo_Data_with_punishment_and_puff.mat';
odor_duration = 2.0; % in seconds
start_tone_duration = 1.00; % in seconds
puff_duration = 0.5; % in seconds

% LOAD DATA from .mat file

% find all .mat files in datafolder and ask user to select one
files = dir(fullfile(datafolder,'*.mat'));
filenames = {files.name};

% keep only files that contain filtername
idx = contains(filenames,filtername);
filenames = filenames(idx);

% sort filenames
[~,idx] = sort(filenames);
filenames = filenames(idx);

% ask user to select one file
[idx,~] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',filenames);
filename = filenames{idx};

% make simple name (remove filtername, '_', 'Date_' from filename)
simplename = strrep(filename,filtername,'');
simplename = strrep(simplename,'Date_','');
simplename = strrep(simplename,'_',' ');

% create a figs folder if it does not exist
if ~exist(strcat(datafolder,'figs/'), 'dir')
    mkdir(strcat(datafolder,'figs/'));
    figsfolder = strcat(datafolder,'figs/');
else
    figsfolder = strcat(datafolder,'figs/');
end

% create a processed folder if it does not exist
if ~exist(strcat(datafolder,'processed/'), 'dir')
    mkdir(strcat(datafolder,'processed/'));
    processedfolder = strcat(datafolder,'processed/');
else
    processedfolder = strcat(datafolder,'processed/');
end

% delete temporary variables
clear files filenames idx;

% load data
load([datafolder filename]);

% get the M field in the exper struct
DATA = exper.M;

% extract the fields of interest (columns of DATA; dont keeo last timestep
% 1: time
% 11: trial number
% 3: start tone
% 4,5: odor 1,2
% 9: water delivered
% 15: lick contact
% 13: wrong

time = DATA(1:end-1,1);
trial = DATA(1:end-1,11);
start_tone = DATA(1:end-1,3);
odor1 = DATA(1:end-1,4);
odor2 = DATA(1:end-1,5);
water = DATA(1:end-1,9);
lick = DATA(1:end-1,15);
wrong = DATA(1:end-1,13);

% delete temporary variables
clear exper DATA;

% get the unique trial numbers and sort them
trial_numbers = unique(trial);
trial_numbers = sort(trial_numbers);
max_trial_number = max(trial_numbers);


% create a matrix to save the lick rates/latency for each trial
all_lick_rates = zeros(4,length(trial_numbers));
all_lick_latency = zeros(4,length(trial_numbers));
trial_type = [];

figure;
% loop through each trial number and extract the data for that trial
for i = 1:length(trial_numbers)

    % create a array to save notable times for the current trial
    notable_times = [];
    
    % trial number sets the y-axis value
    y = max_trial_number - trial_numbers(i);
    
    % get the index where the trial number is equal to the current trial number
    idx = trial == trial_numbers(i);

    % get the time for the current trial and reorigin it
    time_trial = time(idx);

    % get the time of the start tone for the current trial
    start_tone_trial = start_tone(idx);
    start_tone_time = time_trial(start_tone_trial == 1);

    % reorigin the trial to the start tone
    time_trial = time_trial - start_tone_time(1);
    start_tone_time = start_tone_time - start_tone_time(1);
    
    % save the start tone time as a notable time
    notable_times = [notable_times; start_tone_time];

    % plot the start tone (rectangle at [start_tone_time, y] with width start_tone_duration and height 1)
    rectangle('Position',[start_tone_time,y,start_tone_duration,0.75],'FaceColor','#ff6161','EdgeColor','none');

    % get the time of the odor1 (local time)
    odor1_trial = odor1(idx);
    odor1_time = time_trial(odor1_trial == 1);
    
    % loop through each odor1 time and plot a rectangle at [odor1_time, y] with width odor_duration and height 1
    for j = 1:length(odor1_time)
        rectangle('Position',[odor1_time(j),y,odor_duration,0.75],'FaceColor','#f7ec86','EdgeColor','none');
        % save the odor1 times as a notable time
        notable_times = [notable_times; odor1_time(j)];
    end

    % get the time of the odor2 (local time)
    odor2_trial = odor2(idx);
    odor2_time = time_trial(odor2_trial == 1);

    % loop through each odor2 time and plot a rectangle at [odor2_time, y] with width odor_duration and height 1
    for j = 1:length(odor2_time)
        rectangle('Position',[odor2_time(j),y,odor_duration,0.75],'FaceColor','#fab6ee','EdgeColor','none');
        % save the odor2 times as a notable time
        notable_times = [notable_times; odor2_time(j)];
    end

    % get trial type 
    % if 2 odor 1 trials or 2 odor 2 trials, trial type is "match odor 1/2"
    % if 1 odor 1 trial and 1 odor 2 trial, trial type is "non-match odor 1" if odor 1 is first, "non-match odor 2" if odor 2 is first
    if length(odor1_time) == 2
        trial_type = [trial_type; "match odor 1"];
    elseif length(odor2_time) == 2
        trial_type = [trial_type; "match odor 2"];
    elseif length(odor1_time) == 1 && length(odor2_time) == 1
        if odor1_time < odor2_time
            trial_type = [trial_type; "non-match odor 1"];
        else
            trial_type = [trial_type; "non-match odor 2"];
        end
    else
        trial_type = [trial_type; "error"];
    end

    % get the time of the water (local time)
    water_trial = water(idx);
    water_time = time_trial(water_trial == 1);
    % check if water was delivered
    if ~isempty(water_time)
        water_time_onset = water_time(1);
        water_time_offset = water_time(end);
        % plot a rectangle at [water_time_onset, y] with width water_time_offset - water_time_onset and height 1
        rectangle('Position',[water_time_onset,y,water_time_offset - water_time_onset,0.75],'FaceColor','#a1e8cc','EdgeColor','none');
        % save the water onset time as a notable time
        notable_times = [notable_times; water_time_onset];
    end

    % if filtername has 'and_puff' in it, plot the puff at the time marked by the wrong field
    if contains(filtername,'and_puff')
        % get the time of the wrong (local time)
        wrong_trial = wrong(idx);
        wrong_time = time_trial(wrong_trial == 1);
        % check if wrong was delivered
        if ~isempty(wrong_time)
            wrong_time = wrong_time(1);
            % plot a rectangle at [wrong_time, y] with width puff_duration and height 1
            rectangle('Position',[wrong_time,y,puff_duration,0.75],'FaceColor','#ff6161','EdgeColor','none');
            % save the wrong time as a notable time
            notable_times = [notable_times; wrong_time];
        end
    end

    % save the last time of the trial as a notable time
    notable_times = [notable_times; time_trial(end)];

    % get the times of the licks (local time)
    lick_trial = lick(idx);
    delta_lick_trial = diff(lick_trial)>2;
    lick_times = time_trial(1:end-1);
    lick_times = lick_times(delta_lick_trial);

    % loop through each lick time and plot a rectangle at [lick_time, y] with width 0.05 and height 1
    for j = 1:length(lick_times)
        rectangle('Position',[lick_times(j),y,0.05,0.75],'FaceColor','black','EdgeColor','none');
    end

    % sort the notable times
    notable_times = sort(notable_times);

    % measure the lick rate between each notable time
    lick_rate = zeros(length(notable_times)-2,1);
    for j = 1:length(notable_times)-1
        lick_rate(j) = sum(lick_times > notable_times(j) & lick_times < notable_times(j+1));
        lick_rate(j) = lick_rate(j)/(notable_times(j+1) - notable_times(j));
    end

    % add the lick rates to the all_lick_rates
    all_lick_rates(1:length(lick_rate),i) = lick_rate;

    % measure the lick latency between each notable time
    lick_latency = zeros(length(notable_times)-2,1);
    for j = 1:length(notable_times)-1
        latency = min(lick_times(lick_times > notable_times(j) & lick_times < notable_times(j+1)) - notable_times(j));
        if isempty(latency)
            latency = NaN;
        end
        lick_latency(j) = latency;
    end

    % add the lick latencies to the all_lick_latencies
    all_lick_latency(1:length(lick_latency),i) = lick_latency;

    % clear temporary variables
    clear idx time_trial start_tone_trial odor1_trial odor2_trial water_trial lick_trial delta_lick_trial;

    % end of trial loop
end

% set the y-axis labels
yticks(0.5:max_trial_number-0.5);
yticklabels(flip(trial_numbers));
xlabel('Time (s)');
ylabel('Trial Number');
title(simplename,"Interpreter","none");

% save the figure
saveas(gcf,strcat(figsfolder,strrep(simplename,' ','_'),".png"));


% show the plot
hold off;

% sort all_lick_rates by trial type
[trial_type,idx] = sort(trial_type);
all_lick_rates = all_lick_rates(:,idx);

% plot the lick rates for each trial type
figure;
unique_trial_type = unique(trial_type);
for i = 1:length(unique_trial_type)
    subplot(2,2,i);
    plot(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),"Marker","o");
    % plot the mean and standard deviation of the lick rates for each trial type
    hold on;
    plot(mean(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),2),'LineWidth',2, 'Color', 'black');
    plot(mean(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),2) + std(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),0,2),'LineWidth',2,'LineStyle','--','Color', 'black');
    plot(mean(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),2) - std(all_lick_rates(:,strcmp(trial_type,unique_trial_type(i))),0,2),'LineWidth',2,'LineStyle','--','Color', 'black');
    title(unique_trial_type(i));
    xlabel('Time (s)');
    ylabel('Lick Rate (Hz)');
    ylim([0 max(max(all_lick_rates))]);
    % only 3 x-axis labels needed for non-match odor 1 and non-match odor 2
    if strcmp(unique_trial_type(i),"non-match odor 1") || strcmp(unique_trial_type(i),"non-match odor 2")
        xticks([1 2 3]);
        xticklabels({'After Start','After Odor 1','After Odor 2'});
        xlim([1,3]);
    else
        xticks([1 2 3 4]);
        if contains(filtername,'and_puff')
            xticklabels({'After Start','After Odor 1','After Odor 2','After Puff'});
        else
            xticklabels({'After Start','After Odor 1','After Odor 2','After Water'});
        end
        xlim([1,4]);
    end
end
sgtitle(simplename + " Lick Rates","Interpreter","none");

% save the figure
saveas(gcf,strcat(figsfolder,strrep(simplename,' ','_'),"_lickrate.png"));


% show the plot
hold off;

% plot the lick latency for each trial type
figure;
unique_trial_type = unique(trial_type);
for i = 1:length(unique_trial_type)
    sub = subplot(2,2,i);
    plot(all_lick_latency(:,strcmp(trial_type,unique_trial_type(i))),"Marker","o");
    title(unique_trial_type(i));
    xlabel('Time (s)');
    ylabel('Lick Latency (s)');
    ylim([0 max(max(all_lick_latency))]);
    % only 3 x-axis labels needed for non-match odor 1 and non-match odor 2
    if strcmp(unique_trial_type(i),"non-match odor 1") || strcmp(unique_trial_type(i),"non-match odor 2")
        xticks([1 2 3]);
        xticklabels({'After Start','After Odor 1','After Odor 2'});
        xlim([1,3]);
    else
        xticks([1 2 3 4]);
        xticklabels({'After Start','After Odor 1','After Odor 2','After Water'});
        xlim([1,4]);
    end
end

sgtitle(simplename + " Lick Latency","Interpreter","none");

% save the figure
saveas(gcf,strcat(figsfolder,strrep(simplename,' ','_'),"_licklatency.png"));

% save lick data as a csv file
writematrix(all_lick_rates,strcat(processedfolder,strrep(simplename,' ','_'),"_lickrate.csv"));
writematrix(all_lick_latency,strcat(processedfolder,strrep(simplename,' ','_'),"_licklatency.csv"));
writematrix(trial_type,strcat(processedfolder,strrep(simplename,' ','_'),"_trialtype.csv"));


