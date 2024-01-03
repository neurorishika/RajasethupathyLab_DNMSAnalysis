% code must be compatible with Matlab 2014b

function code = Tunel
% Tunel   Code for the ViRMEn experiment Tunel.
%   code = Tunel   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT


% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)
    clc; % clear display

    % ***********************
    % USER DEFINED PARAMETERS
    % ***********************

    vr.debugMode = false; 

    % Animal Parameters
    % =================

    % REMEMBER TO CHANGE THESE PARAMETERS FOR EACH ANIMAL
    vr.name = 'test'; % Enter animal name
    vr.cagenum = 'test'; % Enter cage number

    % Experiment Parameters
    % =====================

    % REMEMBER TO CHECK THESE PARAMETERS FOR EACH EXPERIMENT
    vr.rewarded_task = 'non-match';  % Whether 'match' or 'non-match' trials are rewarded
    vr.additional_str = '';  % Additional notes to be added to the file name (optional)
    
    % PRESET OPTIONS: 
    % ---------------
    % - '' : default parameters
    % - 'water_release' : water released on left lickport to prime the rig for water delivery
    % - 'interleaved_refresher' : interleaved trials with only the rewarded stimulus (no lick requirement or punishment)
    % - 'short_blocked_habituation' : blocked trials (no lick requirement or punishment); 20 trials
    % - 'long_blocked_habituation' : blocked trials (no lick requirement or punishment); 40 trials
    % - 'short_blocked_lick_training' : blocked trials (lick requirement, no punishment); 20 trials
    % - 'long_blocked_lick_training' : blocked trials  (lick requirement, no punishment); 40 trials
    % - 'short_blocked_contrast_training' : blocked trials (lick requirement, air puff as punishment); 20 trials
    % - 'long_blocked_contrast_training' : blocked trials (lick requirement, air puff as punishment); 40 trials
    % - 'short_test' : interleaved trials (lick requirement, no punishment); 20 trials
    % - 'long_test' : interleaved trials (lick requirement, no punishment); 40 trials

    vr.use_preset = ''; % Enter the name of the preset to use (leave empty for default)

    % DETAILED DESCRIPTION OF THE EXPERIMENT PARAMETERS
    % -------------------------------------------------
    % - vr.experimentName : name of the experiment
    % - vr.outputPath : path to the folder where the data will be saved (data and experiment files will be saved in separate files)
    % - vr.lick_requirement : whether the animal needs to lick to get reward
    % - vr.puff_punishment : whether to deliver air puff as punishment
    % - vr.actuator_active : whether to use the actuator or not
    % - vr.rewarded_task : whether 'match' or 'non-match' trials are rewarded
    % - vr.blocked_vs_interleaved : whether the experiment will be blocked or interleaved
    % - vr.ntrials : total number of trials
    % - vr.blockSize : number of trials per block (for blocked trials)
    % - vr.matched_bias : probability of matched trials 
    %       FOR BLOCKED TRIALS : 0 for only non-matched trials, 1 for only matched trials, 0.5 for 50/50
    %       FOR INTERLEAVED TRIALS : any number between 0 and 1 for the probability of matched trials
    % - vr.stim_bias : probability of stim 2 first trials
    %       FOR BLOCKED TRIALS : 0 for stim 1 first, 1 for stim 2 first, 0.5 for 50/50
    %       FOR INTERLEAVED TRIALS : any number between 0 and 1 for the probability of stim 2 first trials
    % - vr.lickportBias : probability of left vs right lickport (0.5 for either side, 0 for left, 1 for right)
    % - vr.startTone_to_stim1_delay : time between the start tone and the first stim
    % - vr.stim_wait : time the odor is presented for (also defined in arduino code)
    % - vr.decision_period_delay : time into the odor presentation the decision period starts
    % - vr.decision_period : time the decision period is for
    % - vr.reward_period : time the reward period is for
    % - vr.ITI : time the intertrial interval is for
    % - vr.stim1_to_stim2_delay_mean : mean of the Gaussian distribution from which the delay will be drawn
    % - vr.stim1_to_stim2_delay_std : standard deviation of the Gaussian distribution from which the delay will be drawn
    % - vr.probe_trials : probability of probe trials (0 to 1; 0.5 means 50% of trials are probe trials)
    % - vr.integrate_licks : true if lick signal is integrated, false if only onsets are counted


    % ==================
    % DEFAULT PARAMETERS
    % ==================

    vr.outputPath = 'C:\Users\VR Rig #2\Documents\MATLAB\Rishika_data\';    % Path to the folder where the data will be saved (data and experiment files will be saved in separate files)
    

    % EXPERIMENT PARAMETERS
    vr.experimentName = '';                                                 % Name of the experiment
    vr.lick_requirement = false;                                            % Require licking to get reward
    vr.puff_punishment = false;                                             % Deliver air puff as punishment instead of just not giving water
    
    % TRIAL PARAMETERS
    vr.blocked_vs_interleaved = 'blocked';                                  % Define if the experiment will be blocked or interleaved
    vr.ntrials = 20;                                                        % Define the total number of trials
    vr.blockSize = 5;                                                       % for blocked trials, this is the number of trials per block
    vr.matched_bias = 0.5;                                                  % Define the probability of matched trials
    vr.stim_bias = 0.5;                                                     % Define the probability of stim 2 first trials
    vr.probe_trials = 0;                                                    % probability of probe trials (0 to 1; 0.5 means 50% of trials are probe trials)

    % RIG PARAMETERS
    vr.lickportBias = 0;                                                    % lickport bias (0.5 for either side, 0 for left, 1 for right)
    vr.integrate_licks = false;                                             % true if lick signal is integrated, false if only onsets are counted
    vr.actuator_active = false;                                             % Whether to use the actuator or not
    
    % TIME PARAMETERS
    vr.startTone_to_stim1_delay = 1.3;                                      % Define the time between the start tone and the first stim
    vr.stim_wait = 3.4;                                                     % How long the odor is presented for (also defined in arduino code)
    vr.decision_period_delay = 0.9;                                         % How long into the odor presentation the decision period starts
    vr.decision_period = 1;                                                 % How long the decision period is for
    vr.reward_period = 4;                                                   % How long the reward period is for
    vr.ITI = 15;                                                            % How long the intertrial period is for
    vr.stim1_to_stim2_delay_mean = 1;                                       % Define the mean of the Gaussian distribution from which the delay will be drawn
    vr.stim1_to_stim2_delay_std = 0.0;                                      % Define the standard deviation of the Gaussian distribution from which the delay will be drawn

    % ---------------------
    % PRESET CONFIGURATIONS
    % ---------------------

    % water release
    if strcmp(vr.use_preset,'water_release')
        vr.experimentName = 'rig_test';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 10;
        vr.blockSize = 5;
        if strcmp(vr.rewarded_task,'match')
            vr.matched_bias = 1;
        elseif strcmp(vr.rewarded_task,'non-match')
            vr.matched_bias = 0;
        end
        vr.reward_period = 10;
    end

    % interleaved refresher
    if strcmp(vr.use_preset,'interleaved_refresher')
        vr.experimentName = 'refresher';
        vr.blocked_vs_interleaved = 'interleaved';
        vr.ntrials = 6;
        if strcmp(vr.rewarded_task,'match')
            vr.matched_bias = 1;
        elseif strcmp(vr.rewarded_task,'non-match')
            vr.matched_bias = 0;
        end
    end

    % short blocked habituation
    if strcmp(vr.use_preset,'short_blocked_habituation')
        vr.experimentName = 'habituation';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 20;
        vr.blockSize = 5;
    end

    % long blocked habituation
    if strcmp(vr.use_preset,'long_blocked_habituation')
        vr.experimentName = 'habituation';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 40;
        vr.blockSize = 5;
    end

    % short blocked lick training
    if strcmp(vr.use_preset,'short_blocked_lick_training')
        vr.experimentName = 'lick_training';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 20;
        vr.blockSize = 5;
        vr.lick_requirement = true;
    end

    % long blocked lick training
    if strcmp(vr.use_preset,'long_blocked_lick_training')
        vr.experimentName = 'lick_training';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 40;
        vr.blockSize = 5;
        vr.lick_requirement = true;
    end

    % short blocked contrast training
    if strcmp(vr.use_preset,'short_blocked_contrast_training')
        vr.experimentName = 'contrast_training';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 20;
        vr.blockSize = 5;
        vr.lick_requirement = true;
        vr.puff_punishment = true;
    end

    % long blocked contrast training
    if strcmp(vr.use_preset,'long_blocked_contrast_training')
        vr.experimentName = 'contrast_training';
        vr.blocked_vs_interleaved = 'blocked';
        vr.ntrials = 40;
        vr.blockSize = 5;
        vr.lick_requirement = true;
        vr.puff_punishment = true;
    end

    % short test
    if strcmp(vr.use_preset,'short_test')
        vr.experimentName = 'test';
        vr.blocked_vs_interleaved = 'interleaved';
        vr.ntrials = 20;
        vr.lick_requirement = true;
    end

    % long test
    if strcmp(vr.use_preset,'long_test')
        vr.experimentName = 'test';
        vr.blocked_vs_interleaved = 'interleaved';
        vr.ntrials = 40;
        vr.lick_requirement = true;
    end



    % ***************************
    % DEFINE THE OUTPUT STRUCTURE
    % ***************************


    % Define the output file name (FORMAT: name_cagenum_mmddyyyy_<experiment_name>_<lick requirement>_<puff punishment>.mat)
    if vr.lick_requirement == true
        vr.str_lickreq = '_lick'; % Add this to file name if lick requirement is true
    else
        vr.str_lickreq = ''; % Add this to file name if lick requirement is true
    end

    if vr.puff_punishment == true
        vr.str_puff = '_puff'; % Add this to file name if puff punishment is true
    else
        vr.str_puff = ''; % Add this to file name if puff punishment is true
    end

    % add DNMS/DMS to the experiment name
    if strcmp(vr.rewarded_task,'match')
        vr.experimentName = strcat('DMS_gonogo_',vr.experimentName);
    elseif strcmp(vr.rewarded_task,'non-match')
        vr.experimentName = strcat('DNMS_gonogo_',vr.experimentName);
    end


    vr.outputFile = strcat(vr.name,'_',vr.cagenum,'_',datestr(now,'mmddyyyy'),'_',vr.experimentName,vr.str_lickreq,vr.str_puff,vr.additional_str);

    % Save the ViRMEn experiment structure to a separate file
    exper = vr.exper;

    % check if output file already exists
    if exist([vr.outputPath vr.outputFile '_Experiment.mat'],'file')
        % recusively increment number and check again
        i = 1;
        while exist([vr.outputPath vr.outputFile '_Experiment_' num2str(i) '.mat'],'file')
            i = i+1;
        end
        vr.name_addendum = strcat('_',num2str(i));
        disp(strcat('Experiment already exists, saving with additional marker ',vr.name_addendum,'. Please check and rename after completion.'));
    else
        vr.name_addendum = '';
    end
    save([vr.outputPath vr.outputFile '_Experiment' vr.name_addendum '.mat'],'exper');
    
    % save the experiment metadata as a separate file (using a structure)
    experiment_metadata.outputFile = vr.outputFile;
    experiment_metadata.outputPath = vr.outputPath;
    experiment_metadata.name = vr.name;
    experiment_metadata.cagenum = vr.cagenum;
    experiment_metadata.experimentName = vr.experimentName;
    experiment_metadata.lick_requirement = vr.lick_requirement;
    experiment_metadata.puff_punishment = vr.puff_punishment;
    experiment_metadata.blocked_vs_interleaved = vr.blocked_vs_interleaved;
    experiment_metadata.ntrials = vr.ntrials;
    experiment_metadata.blockSize = vr.blockSize;
    experiment_metadata.matched_bias = vr.matched_bias;
    experiment_metadata.stim_bias = vr.stim_bias;
    experiment_metadata.lickportBias = vr.lickportBias;
    experiment_metadata.startTone_to_stim1_delay = vr.startTone_to_stim1_delay;
    experiment_metadata.stim1_to_stim2_delay_mean = vr.stim1_to_stim2_delay_mean;
    experiment_metadata.stim1_to_stim2_delay_std = vr.stim1_to_stim2_delay_std;
    experiment_metadata.decision_period_delay = vr.decision_period_delay;
    experiment_metadata.decision_period = vr.decision_period;
    experiment_metadata.reward_period = vr.reward_period;
    experiment_metadata.ITI = vr.ITI;
    experiment_metadata.probe_trials = vr.probe_trials;
    experiment_metadata.integrate_licks = vr.integrate_licks;
    experiment_metadata.str_punish = vr.str_punish;
    experiment_metadata.str_puff = vr.str_puff;
    experiment_metadata.additional_str = vr.additional_str;

    % check if metadata file already exists
    if exist([vr.outputPath vr.outputFile '_Metadata' vr.name_addendum '.mat'],'file')
        % recusively increment number and check again
        i = 1;
        while exist([vr.outputPath vr.outputFile '_Metadata' vr.name_addendum '_' num2str(i) '.mat'],'file')
            i = i+1;
        end
        vr.addl_name_addendum = strcat('_',num2str(i));
        disp(strcat('Metadata file already exists, saving with additional marker ',vr.addl_name_addendum,'. Please check and rename after completion.'));
    else
        vr.addl_name_addendum = '';
    end
    save([vr.outputPath vr.outputFile '_Metadata' vr.name_addendum vr.addl_name_addendum '.mat'],'experiment_metadata');
    

    % ****************************
    % INITIALIZATION OF EXPERIMENT
    % ****************************

    % Get data from  input channels via NI-DAQ

    % Start analog session to recieve inputs from two lickports
    vr.analogSess = daq.createSession('ni');
    addAnalogInputChannel(vr.analogSess,'Dev1',[0 1], 'Voltage'); %device object, device ID, channel numbers, measurement type (voltage)

    % Start digital session to control stimulus delivery
    vr.DigitalSess = daq.createSession('ni');
    addDigitalChannel(vr.DigitalSess,'dev1', 'Port0/Line0:7', 'OutputOnly');
    addDigitalChannel(vr.DigitalSess,'dev1', 'Port1/Line0:1', 'OutputOnly');
    outputSingleScan(vr.DigitalSess,[0 0 0 0 0 0 0 0 0 0]);

    % Define triggering variables for the stimuli
    vr.newtrial = 0;                            % used to initiate a new trial
    vr.tone = 0;                                % used to trigger a start tone (implemented in arduino code)
    vr.Stim1_valve = 0;                         % used to trigger odor 1 and sound 1 (implemented in arduino code)
    vr.Stim2_valve = 0;                         % used to trigger odor 2 and sound 2 (implemented in arduino code)
    vr.AirON = 0;                               % used to trigger airflow for puff punishment
    vr.Extend_actuator = 0;                     % used to trigger actuator extension
    vr.Retract_actuator = 0;                    % used to trigger actuator retraction
    vr.Deliver_water_left = 0;                  % used to trigger water delivery to left lickport
    vr.Deliver_water_right = 0;                 % used to trigger water delivery to right lickport
    vr.synchPulse = 0;                          % used to trigger a synch pulse for imaging
    vr.puffPun = 0;                             % used to trigger a puff of air as punishment

    % Define variables which enable execution of proper trial chronology
    vr.tone_presented = 0;               % check if tone has been presented
    vr.first_stim_started = 0;     % check if stim 1 has been presented
    vr.delay_started = 0;          % check if delay has started
    vr.second_stim_begin = 0;      % check if stim 2 has begun
    vr.ITI_started = 0;           % check if ITI has happened
    vr.decision_started = 0;       % check if decision period has begun
    vr.reward_started = 0;         % check if reward period has begun
    vr.stim1_first = 0;             % 1 if stim 1 was presented first in the trial
    vr.stim2_first = 0;             % 1 if stim 2 was presented first in the trial
    vr.stim1_second = 0;            % 1 if stim 1 was presented second in the trial
    vr.stim2_second = 0;            % 1 if stim 2 was presented second in the trial
    vr.end_trial_early = 0;       % 1 if trial was ended early
    vr.wrong = 0;                  % 1 if the animal made a wrong decision
    vr.correct = 0;               % 1 if the animal made a correct decision

    % Define variables which keep track of the trial chronology
    vr.stopwatch = 0; % variable to keep track of time
    vr.timepoint_to_track = 0; % variable to mark the timepoint to track
    vr.current_trial = 1; % variable to keep track of the current trial

    % Define variables to store the time of important events
    vr.t_tone = 0; % time of the start tone
    vr.t_stim1_start = 0; % time of the start of stim 1
    vr.t_delay_start = 0; % time of the start of the delay
    vr.t_stim2_start = 0; % time of the start of stim 2
    vr.t_decision_start = 0; % time of the start of the decision period
    vr.t_reinforcement_start = 0; % time of the start of the reinforcement period
    vr.t_ITI_start = 0; % time of the start of the intertrial interval

    % Define variables to track whether reinforcement was delivered
    vr.reward_delivered = 0; % 1 if reward was delivered
    vr.punishment_delivered = 0; % 1 if punishment was delivered

    % Define variables to store lick data
    vr.licks_right = 0;
    vr.licks_left = 0;
    vr.total_licks = 0;
    vr.lick_signal = [0 0];
    vr.last_lick_signal = inputSingleScan(vr.analogSess);
    vr.nogo_licks = 0;
    vr.go_licks = 0;

    % HOW TRIALS ARE DEFINED
    % We create the vector of different conditions of length which is eaqual to the number of trials where:
    % - Condition 1 ---> Stim1 - Delay - Stim2 (Stim 1 Non-matched)
    % - Condition 2 ---> Stim2 - Delay - Stim1 (Stim 2 Non-matched)
    % - Condition 3 ---> Stim1 - Delay - Stim1 (Stim 1 Matched)
    % - Condition 4 ---> Stim2 - Delay - Stim2 (Stim 2 Matched)
    
    % Blocked vs Interleaved

    switch vr.blocked_vs_interleaved
        case 'blocked'
            % BLOCK TRIALS
            
            vr.matched_bias = round(vr.matched_bias*2)/2; % round to the closest of 0, 0.5 or 1
            vr.stim_bias = round(vr.stim_bias*2)/2; % round to the closest of 0, 0.5 or 1

            % Define the conditions which will be used in the experiment
            if and(vr.matched_bias == 0.5, vr.stim_bias == 0.5)
                vr.rand = [1 2 3 4];
            elseif and(vr.matched_bias == 0.5, vr.stim_bias == 0)
                vr.rand = [1, 3];
            elseif and(vr.matched_bias == 0.5, vr.stim_bias == 1)
                vr.rand = [2, 4];
            elseif and(vr.matched_bias == 0, vr.stim_bias == 0.5)
                vr.rand = [1, 2];
            elseif and(vr.matched_bias == 1, vr.stim_bias == 0.5)
                vr.rand = [3, 4];
            elseif and(vr.matched_bias == 0, vr.stim_bias == 0)
                vr.rand = [1];
            elseif and(vr.matched_bias == 0, vr.stim_bias == 1)
                vr.rand = [2];
            elseif and(vr.matched_bias == 1, vr.stim_bias == 0)
                vr.rand = [3];
            elseif and(vr.matched_bias == 1, vr.stim_bias == 1)
                vr.rand = [4];
            end
            
            % make sure that the number of trials is divisible by the blockSize*4 (number of conditions)
            % otherwise it will be rounded to the closest number divisible by the blockSize*4
            if mod(vr.ntrials,vr.blockSize*length(vr.rand)) ~= 0
                vr.ntrials = ceil(vr.ntrials/(vr.blockSize*length(vr.rand)))*(vr.blockSize*length(vr.rand));
                disp(strcat('More trials needed for coverage. ',num2str(vr.ntrials),' trials will be run.'));
            end

            vr.j = 1; % iterator
            vr.conditions = []; % the conditions for each trial will be stored here

            while vr.j <= vr.ntrials
                vr.rand = vr.rand(randperm(length(vr.rand))); % shuffle the conditions
                for i = 1:length(vr.rand)
                    vr.conditions = [vr.conditions; ones(vr.blockSize,1)*vr.rand(i)]; % add the conditions to the vector
                    vr.j = vr.j + vr.blockSize; % increment the iterator
                end
            end
        
        case 'interleaved'
            % INTERLEAVE TRIALS
            
            % determine number of stim 1 vs stim 2 first trials
            vr.n_stim1_first = round(vr.stim_bias*vr.ntrials);
            vr.n_stim2_first = vr.ntrials - vr.n_stim1_first;

            % determine number of matched vs non-matched trials for each stim 1 vs stim 2 first trial
            vr.n_stim1_first_matched = round(vr.matched_bias*vr.n_stim1_first);
            vr.n_stim1_first_nonmatched = vr.n_stim1_first - vr.n_stim1_first_matched;
            vr.n_stim2_first_matched = round(vr.matched_bias*vr.n_stim2_first);
            vr.n_stim2_first_nonmatched = vr.n_stim2_first - vr.n_stim2_first_matched;

            % create the conditions vector
            vr.conditions = [ones(vr.n_stim1_first_matched,1)*3; ones(vr.n_stim1_first_nonmatched,1)*1; ones(vr.n_stim2_first_matched,1)*4; ones(vr.n_stim2_first_nonmatched,1)*2];
            
            % shuffle the conditions vector
            vr.conditions = vr.conditions(randperm(length(vr.conditions)));

    end

    % And finally within the same block create a vector of delays whose length is equal to the number
    % of trials by rounding an absolute value of a number taken from the Gaussian distribution with 
    % mean of vr.stim1_to_stim2_delay_mean and std of deviation of vr.stim1_to_stim2_delay_std
    if vr.stim1_to_stim2_delay_std <= 0
        vr.stim1_to_stim2_delays = ones(vr.ntrials,1)*vr.stim1_to_stim2_delay_mean;
    else
        vr.stim1_to_stim2_delays = ones(vr.ntrials,1)*vr.stim1_to_stim2_delay_mean + round(vr.stim1_to_stim2_delay_std*randn(vr.ntrials,1));
    end
    vr.stim1_to_stim2_delays(vr.stim1_to_stim2_delays < 0) = 0; % clip the delays to be positive

    % Probe trials for used for  shaping without decision period such that the animal is 
    vr.probetrials = [ones(round(vr.probe_trials*vr.ntrials),1); zeros(vr.ntrials-round(vr.probe_trials*vr.ntrials),1)];
    vr.probetrials = vr.probetrials(randperm(length(vr.probetrials)));

    % LOG OUTPUT
    % ==========
    disp('Experiment started');
    disp('==================');
    if vr.debugMode == true
        disp('DEBUG MODE ON');
        disp('Trial conditions:');
        disp(vr.conditions);
        disp('Trial delays:');
        disp(vr.stim1_to_stim2_delays);
        disp('Probe trials:');
        disp(vr.probetrials);
    else
        disp('DEBUG MODE OFF');
    end

    % Initiaize the Imaging system by sending a synch pulse to NI-DAQ (NEEDED?)
    if (round(vr.iterations) == 100)
        vr.synchPulse=1;
        disp('output dout')
    end    



% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
    % Annul the stimuli triggering variables at the start of each trial
    vr.Stim1_valve = 0;
    vr.Stim2_valve = 0;
    vr.tone = 0;
    vr.AirON = 0;
    vr.lick_signal = [0 0];
    vr.Extend_actuator = 0;
    vr.Retract_actuator = 0;
    vr.synchPulse = 0;

    % ***************************
    % START OF TRIAL + START TONE
    % ***************************

    % start experiment after 100 iterations
    if (round(vr.iterations) == 100) || (vr.newtrial == 1)
        % Send a synch pulse to NI-DAQ at the start of the first trial
        if vr.current_trial == 1;
            % Send a synch pulse to NI-DAQ and display that it has been sent
            vr.synchPulse=1; 
            disp('output dout');
        end
        % Display the trial number and the condition+delay (if debug mode is on)
        disp(strcat('Initiated trial: ',num2str(vr.current_trial)));
        if vr.debugMode == true
            disp(strcat('Condition: ',num2str(vr.conditions(vr.current_trial))));
            disp(strcat('Delay: ',num2str(vr.stim1_to_stim2_delays(vr.current_trial))));
        end
        % Setup the start tone and initiate the trial
        vr.tone = 1; % initiate the start tone
        disp(strcat('Start tone sent at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        vr.newtrial = 0;  % new trial started, set to 0
        vr.tone_presented = 1;  % variable to verify if tone has been presented
        vr.timepoint_to_track = vr.timeElapsed; % track the time of the start tone
        vr.t_tone = vr.timepoint_to_track; % store the time of the start tone
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; 


    % Get the analog signal from the lickports and measure the licking
    vr.lick_signal=inputSingleScan(vr.analogSess);

    % ***************
    % START OF STIM 1
    % ***************

    % If the delay between the start tone and the first stim has passed and tone has been presented
    % present the first stim depending on the condition for this particular trial
    if (round(vr.stopwatch,1) == vr.startTone_to_stim1_delay) && (vr.tone_presented == 1)
        % Present the first stim depending on the condition for this particular trial
        switch vr.conditions(vr.current_trial)
            case {1,3}
                % Stim 1 is presented first
                vr.Stim1_valve = 1; % initiate stim 1 presentation
                vr.stim1_first = 1; % tag that stim 1 was presented first
                disp(strcat('First Stimulus: Stim1 started at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
            case {2,4}
                % Stim 2 is presented first
                vr.Stim2_valve = 1; % initiate stim 2 presentation
                vr.stim2_first = 1; % tag that stim 2 was presented first
                disp(strcat('First Stimulus: Stim2 started at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        end
        vr.tone_presented = 0; % tone has been presented, set to 0
        vr.first_stim_started = 1; % variable to verify if stim 1 has been presented
        vr.timepoint_to_track = vr.timeElapsed;  % reset the timepoint to track the start of stim 1 presentation
        vr.t_stim1_start = vr.timepoint_to_track; % store the time of the start of stim 1 presentation
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track;

    % **************
    % START OF DELAY
    % **************

    % Wait for stim 1 to be presented and then start the delay period
    if (round(vr.stopwatch,1) == vr.stim_wait) && (vr.first_stim_started == 1)
        disp(strcat('Delay started for ',num2str(round(vr.stim1_to_stim2_delays(vr.current_trial),1)), ' s'));
        vr.first_stim_started = 0; % stim 1 has been presented, set to 0
        vr.delay_started = 1; % variable to verify if delay has started
        vr.timepoint_to_track = vr.timeElapsed; % reset the timepoint to track the start of the delay
        vr.t_delay_start = vr.timepoint_to_track; % store the time of the start of the delay
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; % track the time of the delay

    % ***************
    % START OF STIM 2
    % ***************

    % Wait for the delay to pass and then present the second stim depending on the condition for this particular trial
    if (round(vr.stopwatch,1) == vr.stim1_to_stim2_delays(vr.current_trial)) && (vr.delay_started == 1)
        switch vr.conditions(vr.current_trial)
            case {1,4}
                % Stim 2 is presented second
                vr.Stim2_valve = 1; % initiate stim 2 presentation
                vr.stim2_second = 1; % tag that stim 2 was presented second
                disp(strcat('Second Stimulus: Stim2 started at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
            case {2,3}
                % Stim 1 is presented second
                vr.Stim1_valve = 1; % initiate stim 1 presentation
                vr.stim1_second = 1; % tag that stim 1 was presented second
                disp(strcat('Second Stimulus: Stim1 started at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        end
        vr.delay_started = 0; % delay has passed, set to 0
        vr.second_stim_begin = 1; % variable to verify if stim 2 has begun
        vr.timepoint_to_track = vr.timeElapsed; % reset the timepoint to track the start of stim 2 presentation
        vr.t_stim2_start = vr.timepoint_to_track; % store the time of the start of stim 2 presentation
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; % track the time of the delay

    % ************************
    % START OF DECISION PERIOD
    % ************************

    % Initiation of the decision point. (handle actuators and water delivery for probe trials)
    if (round(vr.stopwatch,1) == vr.decision_period_delay) && (vr.second_stim_begin == 1)
        % Display that the decision period has started
        disp(strcat('Decision Point Started at t = ',num2str(round(vr.timeElapsed,1))));
        % Extend the actuator if it is active
        if vr.actuator_active == true
            vr.Extend_actuator = 1;
        end
        % Deliver water if it is a probe trial
        if vr.probetrials(vr.current_trial) == 1
            disp('Probe trial');
            switch vr.conditions(vr.current_trial)
                case {1,2} % Non-matched trials
                    % check if non-matched trial are rewarded
                    if strcmp(vr.rewarded_task,'non-match')
                        if vr.lickportBias == 0.5
                            vr.Deliver_water_right = 1;
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to both lickports (non-match probe trial)');
                            vr.reward_delivered = 1;
                        elseif vr.lickportBias == 0
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to left lickport (non-match probe trial)');
                            vr.reward_delivered = 1;
                        elseif vr.lickportBias == 1
                            vr.Deliver_water_right = 1;
                            disp('Delivering water to right lickport (non-match probe trial)');
                            vr.reward_delivered = 1;
                        end
                    else
                        vr.Deliver_water_left = 0;
                        vr.Deliver_water_right = 0;
                        disp('No water delivered (non-match probe trial)');
                    end
                case {3,4} % Matched trials
                    % check if matched trial are rewarded
                    if strcmp(vr.rewarded_task,'match')
                        if vr.lickportBias == 0.5
                            vr.Deliver_water_right = 1;
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to both lickports (match probe trial)');
                            vr.reward_delivered = 1;
                        elseif vr.lickportBias == 0
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to left lickport (match probe trial)');
                            vr.reward_delivered = 1;
                        elseif vr.lickportBias == 1
                            vr.Deliver_water_right = 1;
                            disp('Delivering water to right lickport (match probe trial)');
                            vr.reward_delivered = 1;
                        end
                    else
                        vr.Deliver_water_left = 0;
                        vr.Deliver_water_right = 0;
                        disp('No water delivere (match probe trial)');
                    end
            end
            % end the trial early
            vr.end_trial_early = 1;
        end
        % Update the variables to track the decision period
        vr.second_stim_begin = 0; % stim 2 was presented and decision period has begun, set to 0
        vr.decision_started = 1; % variable to verify if decision period has begun
        vr.timepoint_to_track = vr.timeElapsed; % reset the timepoint to track the start of the decision period
        vr.t_decision_start = vr.timepoint_to_track; % store the time of the start of the decision period
        % vr.timepoint_to_track = vr.timeElapsed+1; %1s delay once extended to make decision  ??????
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; % track the time of the delay

    % **********************
    % DURING DECISION PERIOD
    % **********************

    % if the decision period has started and the animal has not made a decision yet
    if (vr.decision_started == 1) && (round(vr.stopwatch,1) < vr.decision_period)
        if vr.integrate_licks == true
            if vr.lick_signal(1) > 3 % if the lick signal is above 3V the tongue is in the left lickport
                vr.licks_left = vr.licks_left + 1;
            end
            if vr.lick_signal(2) > 3 % if the lick signal is above 3V the tongue is in the right lickport
                vr.licks_right = vr.licks_right + 1;
            end
        else
            if vr.lick_signal(1) > 3 && vr.last_lick_signal(1) < 3 % lick onset on left lickport
                vr.licks_left = vr.licks_left + 1;
            end
            if vr.lick_signal(2) > 3 && vr.last_lick_signal(2) < 3 % lick onset on right lickport
                vr.licks_right = vr.licks_right + 1;
            end
        end
        vr.last_lick_signal = vr.lick_signal; % update the last lick signal
    end

    % **********************
    % START OF REINFORCEMENT
    % **********************
    
    % End of the decision point. Check for the presence of the punishment and
    % act acordingly depending on the combination of stims being present within
    % the experiment
    if (round(vr.stopwatch,1) == vr.decision_period) && (vr.decision_started == 1)
        % Display that the decision period has ended
        disp(strcat('Decision Point Ended at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        if vr.lickportBias == 0.5
            disp(strcat('Licks Left: ',num2str(vr.licks_left)));
            disp(strcat('Licks Right: ',num2str(vr.licks_right)));
            vr.total_licks = vr.licks_left + vr.licks_right;
        elseif vr.lickportBias == 0
            disp(strcat('Licks Left: ',num2str(vr.licks_left)));
            vr.total_licks = vr.licks_left;
        elseif vr.lickportBias == 1
            disp(strcat('Licks Right: ',num2str(vr.licks_right)));
            vr.total_licks = vr.licks_right;
        end
        % Check if the animal made a decision and act accordingly
        switch vr.conditions(vr.current_trial)
            case {1,2} % Non-matched trials
                % Check if non-matched trial are rewarded
                if strcmp(vr.rewarded_task,'non-match')
                    % this is a go (rewarded) trial
                    vr.go_licks = vr.go_licks + vr.total_licks;
                    % check if lick is required
                    if vr.lick_requirement == true
                        % deliver water if the animal licked
                        if vr.total_licks > 0
                            if vr.lickportBias == 0.5
                                vr.Deliver_water_right = 1;
                                vr.Deliver_water_left = 1;
                                disp('Delivering water to both lickports (non-match lick requirement met)');
                            elseif vr.lickportBias == 0
                                vr.Deliver_water_left = 1;
                                disp('Delivering water to left lickport (non-match lick requirement met)');
                            elseif vr.lickportBias == 1
                                vr.Deliver_water_right = 1;
                                disp('Delivering water to right lickport (non-match lick requirement met)');
                            end
                            vr.correct = 1;
                            vr.reward_delivered = 1;
                        else
                            vr.wrong = 1;
                            disp('No water delivered, ending trial early (non-match lick requirement NOT met)');
                            % end the trial if the animal did not lick
                            vr.end_trial_early = 1;

                        end
                    else
                        % deliver water regardless of licking
                        if vr.lickportBias == 0.5
                            vr.Deliver_water_right = 1;
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to both lickports (non-match no lick requirement)');
                        elseif vr.lickportBias == 0
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to left lickport (non-match no lick requirement)');
                        elseif vr.lickportBias == 1
                            vr.Deliver_water_right = 1;
                            disp('Delivering water to right lickport (non-match no lick requirement)');
                        end
                        vr.reward_delivered = 1;
                        % check if the animal licked
                        if vr.total_licks > 0
                            vr.correct = 1;
                        else
                            vr.wrong = 1;
                        end
                    end
                else
                    % this is a nogo (unrewarded) trial
                    vr.nogo_licks = vr.nogo_licks + vr.total_licks;
                    % check if puff punishment is active
                    if vr.puff_punishment == true
                        % deliver puff of air if the animal licked
                        if vr.total_licks > 0
                            vr.puffPun = 1;
                            disp('Puff punishment delivered, ending trial early (non-matched trial lick observed)');
                            vr.wrong = 1;
                            vr.punishment_delivered = 1;
                            % end the trial if the animal licked
                            vr.end_trial_early = 1;
                        else
                            disp('No puff punishment delivered (non-matched trial no lick observed)');
                            vr.correct = 1;
                        end
                    else
                        % deliver no air puff regardless of licking
                        disp('No puff punishment delivered (non-matched trial no puff punishment)');
                        % check if the animal licked
                        if vr.total_licks > 0
                            vr.wrong = 1;
                        else
                            vr.correct = 1;
                        end
                    end
                end
            case {3,4} % Matched trials
                % Check if matched trial are rewarded
                if strcmp(vr.rewarded_task,'match')
                    % this is a go (rewarded) trial
                    vr.go_licks = vr.go_licks + vr.total_licks;
                    % check if lick is required
                    if vr.lick_requirement == true
                        % deliver water if the animal licked
                        if vr.total_licks > 0
                            if vr.lickportBias == 0.5
                                vr.Deliver_water_right = 1;
                                vr.Deliver_water_left = 1;
                                disp('Delivering water to both lickports (match lick requirement met)');
                            elseif vr.lickportBias == 0
                                vr.Deliver_water_left = 1;
                                disp('Delivering water to left lickport (match lick requirement met)');
                            elseif vr.lickportBias == 1
                                vr.Deliver_water_right = 1;
                                disp('Delivering water to right lickport (match lick requirement met)');
                            end
                            vr.correct = 1;
                            vr.reward_delivered = 1;
                        else
                            vr.wrong = 1;
                            disp('No water delivered, ending trial early (match lick requirement NOT met)');
                            % end the trial if the animal did not lick
                            vr.end_trial_early = 1;
                        end
                    else
                        % deliver water regardless of licking
                        if vr.lickportBias == 0.5
                            vr.Deliver_water_right = 1;
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to both lickports (match no lick requirement)');
                        elseif vr.lickportBias == 0
                            vr.Deliver_water_left = 1;
                            disp('Delivering water to left lickport (match no lick requirement)');
                        elseif vr.lickportBias == 1
                            vr.Deliver_water_right = 1;
                            disp('Delivering water to right lickport (match no lick requirement)');
                        end
                        vr.reward_delivered = 1;
                        % check if the animal licked
                        if vr.total_licks > 0
                            vr.correct = 1;
                        else
                            vr.wrong = 1;
                        end
                    end
                else
                    % this is a nogo (unrewarded) trial
                    vr.nogo_licks = vr.nogo_licks + vr.total_licks;
                    % check if puff punishment is active
                    if vr.puff_punishment == true
                        % deliver puff of air if the animal licked
                        if vr.total_licks > 0
                            vr.puffPun = 1;
                            disp('Puff punishment delivered, ending trial early (matched trial lick observed)');
                            vr.wrong = 1;
                            vr.punishment_delivered = 1;
                            % end the trial if the animal licked
                            vr.end_trial_early = 1;
                        else
                            disp('No puff punishment delivered (matched trial no lick observed)');
                            vr.correct = 1;
                        end
                    else
                        % deliver no air puff regardless of licking
                        disp('No puff punishment delivered (matched trial no puff punishment)');
                        % check if the animal licked
                        if vr.total_licks > 0
                            vr.wrong = 1;
                        else
                            vr.correct = 1;
                        end
                    end
                end
        end
        vr.decision_started = 0; % decision period has ended, set to 0
        vr.reward_started = 1; % variable to verify if decision period has ended
        %  track the timepoint to track the start of the reinforcement period
        vr.timepoint_to_track = vr.timeElapsed; % reset the timepoint to track the start of the reinforcement period
        vr.t_reinforcement_start = vr.timepoint_to_track; % store the time of the start of the reinforcement period
        
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; % track the time from the end of the decision period
    
    % *************************************
    % END OF REINFORCEMENT and START OF ITI
    % *************************************

    % If the trial has ended early end the trial and start the ITI OR if the
    % trial has not ended early and the reward period has passed end the trial
    % and start the ITI
    if vr.end_trial_early == 1 || (vr.reward_started == 1 && round(vr.stopwatch,1) == vr.reward_period)
        % Display that the reward period has ended
        disp(strcat('Reward Period Ended at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        % Retract the actuator if it is active
        if vr.actuator_active == true
            vr.Retract_actuator = 1;
        end
        % reset all the lick and reward variables to their initial values
        vr.Deliver_water_left = 0;
        vr.Deliver_water_right = 0;
        vr.last_lick_signal = [0 0];
        % turn AirON on for puff punishment
        vr.AirON = 1;
        % Start the ITI
        disp(strcat('ITI started at t = ',num2str(round(vr.ITI,1)), ' s'));
        vr.reward_started = 0; % reward period has ended, set to 0
        vr.end_trial_early = 0; % trial has ended, set to 0
        vr.ITI_started = 1; % variable to verify if ITI has started
        vr.timepoint_to_track = vr.timeElapsed; % reset the timepoint to track the start of the ITI
        vr.t_ITI_start = vr.timepoint_to_track; % store the time of the start of the ITI
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track; % track the time from the start of the ITI

    % After 15 seconds ITI ends
    if ((round(vr.stopwatch,1) == vr.ITI) && (vr.ITI_started == 1))
        disp(strcat('ITI ended at t = ',num2str(round(vr.timeElapsed,1)), ' s'));
        disp('------------------------------');
        % turn puff punishment off
        vr.puffPun = 0;

        % ADD ALL TRIAL VARIABLES TO THE TRIAL MATRIX
        % TRIAL MATRIX:
        % 1 - trial number
        % 2 - condition
        % 3 - delay
        % 4 - stim 1 first
        % 5 - stim 2 first
        % 6 - stim 1 second
        % 7 - stim 2 second
        % 8 - time of the start tone
        % 9 - time of the start of stim 1
        % 10 - time of the start of the delay
        % 11 - time of the start of stim 2
        % 12 - time of the start of the decision period
        % 13 - time of the start of the reinforcement period
        % 14 - time of the start of the ITI
        % 15 - licks left
        % 16 - licks right
        % 17 - total licks
        % 18 - correct
        % 19 - wrong
        % 20 - reward delivered
        % 21 - punishment delivered
        % 22 - probe trial


        vr.TrialData(vr.current_trial,1)=vr.current_trial;
        vr.TrialData(vr.current_trial,2)=vr.conditions(vr.current_trial);
        vr.TrialData(vr.current_trial,3)=vr.stim1_to_stim2_delays(vr.current_trial);
        vr.TrialData(vr.current_trial,4)=vr.stim1_first;
        vr.TrialData(vr.current_trial,5)=vr.stim2_first;
        vr.TrialData(vr.current_trial,6)=vr.stim1_second;
        vr.TrialData(vr.current_trial,7)=vr.stim2_second;
        vr.TrialData(vr.current_trial,8)=vr.t_tone;
        vr.TrialData(vr.current_trial,9)=vr.t_stim1_start;
        vr.TrialData(vr.current_trial,10)=vr.t_delay_start;
        vr.TrialData(vr.current_trial,11)=vr.t_stim2_start;
        vr.TrialData(vr.current_trial,12)=vr.t_decision_start;
        vr.TrialData(vr.current_trial,13)=vr.t_reinforcement_start;
        vr.TrialData(vr.current_trial,14)=vr.t_ITI_start;
        vr.TrialData(vr.current_trial,15)=vr.licks_left;
        vr.TrialData(vr.current_trial,16)=vr.licks_right;
        vr.TrialData(vr.current_trial,17)=vr.total_licks;
        vr.TrialData(vr.current_trial,18)=vr.correct;
        vr.TrialData(vr.current_trial,19)=vr.wrong;
        vr.TrialData(vr.current_trial,20)=vr.reward_delivered;
        vr.TrialData(vr.current_trial,21)=vr.punishment_delivered;
        vr.TrialData(vr.current_trial,22)=vr.probetrials(vr.current_trial);

        % reset all the variables to their initial values
        vr.licks_left = 0;
        vr.licks_right = 0;
        vr.total_licks = 0;
        vr.reward_delivered = 0;
        vr.punishment_delivered = 0;
        vr.correct = 0;
        vr.wrong = 0;
        vr.stim1_first = 0;
        vr.stim2_first = 0;
        vr.stim1_second = 0;
        vr.stim2_second = 0;
        vr.t_tone = 0;
        vr.t_stim1_start = 0;
        vr.t_delay_start = 0;
        vr.t_stim2_start = 0;
        vr.t_decision_start = 0;
        vr.t_reinforcement_start = 0;
        vr.t_ITI_start = 0;

        % reset all the trial variables to their initial values
        vr.tone_presented = 0;
        vr.first_stim_started = 0;
        vr.delay_started = 0;
        vr.second_stim_begin = 0;
        vr.decision_started = 0;
        vr.reward_started = 0;
        vr.end_trial_early = 0;
        vr.ITI_started = 0;

        % update the current trial number
        vr.current_trial = vr.current_trial +1;

        % start a new trial
        vr.newtrial = 1; 
    end
    vr.stopwatch = vr.timeElapsed - vr.timepoint_to_track;

    % If all the trials have been executed end the experiment and display that
    % the experiment has ended and display the performance summary
    if (vr.current_trial == vr.ntrials + 1)
        vr.experimentEnded = true;
        disp('================');
        disp('Experiment ended');     
        % show the performance
        disp('================');
        disp('Performance:');
        disp('------------');
        % show the fraction of go and nogo licks
        disp(strcat('Go Licks: ',num2str(vr.go_licks)));
        disp(strcat('NoGo Licks: ',num2str(vr.nogo_licks)));
        % show the fraction of correct trials and the fraction of incorrect trials
        disp(strcat('Mean Correct: ',num2str(mean(vr.TrialData(:,18)))));
        disp(strcat('Mean Wrong: ',num2str(mean(vr.TrialData(:,19)))));
        % show the fraction of rewarded trials and the fraction of punished trials
        disp(strcat('Mean Rewarded: ',num2str(mean(vr.TrialData(:,20)))));
        disp(strcat('Mean Punished: ',num2str(mean(vr.TrialData(:,21)))));
    end

    % At the end of each iteration output a digital signal to NI-DAQ to trigger
    % the stimuli which were defined to be triggered in the previous if
    % statements by setting them to 1. Otherwise if 0 they will not be triggered

    % OUTPUTS:
    % 1 - Air ON
    % 2 - Tone
    % 3 - Stim 1 valve
    % 4 - Stim 2 valve
    % 5 - Deliver water left
    % 6 - Deliver water right
    % 7 - Extend actuator
    % 8 - Retract actuator
    % 9 - Synch pulse
    % 10 - Puff punishment

    outputSingleScan(vr.DigitalSess,[vr.AirON vr.tone vr.Stim1_valve vr.Stim2_valve vr.Deliver_water_left vr.Deliver_water_right vr.Extend_actuator vr.Retract_actuator vr.synchPulse vr.puffPun]);

    % Save all the relevant interation variables to the matrix IterationData
    % ITERATION MATRIX:
    % 1 - time elapsed
    % 2 - current trial
    % 3 - Air ON
    % 4 - Tone
    % 5 - Stim 1 valve
    % 6 - Stim 2 valve
    % 7 - Deliver water left
    % 8 - Deliver water right
    % 9 - Extend actuator
    % 10 - Retract actuator
    % 11 - Synch pulse
    % 12 - Puff punishment
    % 13 - Left Lickport Signal
    % 14 - Right Lickport Signal
    % 15 - Go Licks
    % 16 - NoGo Licks

    vr.IterationData(vr.iterations,1)=vr.timeElapsed;
    vr.IterationData(vr.iterations,2)=vr.current_trial;
    vr.IterationData(vr.iterations,3)=vr.AirON;
    vr.IterationData(vr.iterations,4)=vr.tone;
    vr.IterationData(vr.iterations,5)=vr.Stim1_valve;
    vr.IterationData(vr.iterations,6)=vr.Stim2_valve;
    vr.IterationData(vr.iterations,7)=vr.Deliver_water_left;
    vr.IterationData(vr.iterations,8)=vr.Deliver_water_right;
    vr.IterationData(vr.iterations,9)=vr.Extend_actuator;
    vr.IterationData(vr.iterations,10)=vr.Retract_actuator;
    vr.IterationData(vr.iterations,11)=vr.synchPulse;
    vr.IterationData(vr.iterations,12)=vr.puffPun;
    vr.IterationData(vr.iterations,13)=vr.lick_signal(1);
    vr.IterationData(vr.iterations,14)=vr.lick_signal(2);
    vr.IterationData(vr.iterations,15)=vr.go_licks;
    vr.IterationData(vr.iterations,16)=vr.nogo_licks;


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)

    % Stop all outputs to NI-DAQ by outputing 0 on every channel
    outputSingleScan(vr.DigitalSess,[0 0 0 0 0 0 0 0 0 0]);

    % Terminate both analog and digital sessions in NI-DAQ
    vr.analogSess.stop
    vr.DigitalSess.stop;

    % Make a struct out of the matrix containing variables for each iteration
    % and the matrix containing variables for each trial and save it
    exper = struct('IterationData',vr.IterationData,'TrialData',vr.TrialData);

    % check if output file already exists
    if exist([vr.outputPath vr.outputFile '_Data' vr.name_addendum '.mat'],'file')
        % recusively increment number and check again
        i = 1;
        while exist([vr.outputPath vr.outputFile '_Data' vr.name_addendum '_' num2str(i) '.mat'],'file')
            i = i+1;
        end
        vr.addl_name_addendum = strcat('_',num2str(i));
        disp(strcat('Data already exists, saving with additional marker ',vr.addl_name_addendum,'. Please check and rename after completion.'));
        
    else
        vr.addl_name_addendum = '';
    end
    save([vr.outputPath vr.outputFile '_Data' vr.name_addendum vr.addl_name_addendum '.mat'],'exper');

    % Save the iteration and trial matrices as csv files after adding headers
    IterationData = vr.IterationData;
    TrialData = vr.TrialData;

    % Add headers to the iteration matrix
    IterationData = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;IterationData];
    IterationData(1,:) = {'Time Elapsed','Current Trial','Air ON','Tone','Stim 1 Valve','Stim 2 Valve','Deliver Water Left','Deliver Water Right','Extend Actuator','Retract Actuator','Synch Pulse','Puff Punishment','Left Lickport Signal','Right Lickport Signal','Go Licks','NoGo Licks'};

    % Add headers to the trial matrix
    TrialData = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;TrialData];
    TrialData(1,:) = {'Trial Number','Condition','Delay','Stim 1 First','Stim 2 First','Stim 1 Second','Stim 2 Second','Time of the Start Tone','Time of the Start of Stim 1','Time of the Start of the Delay','Time of the Start of Stim 2','Time of the Start of the Decision Period','Time of the Start of the Reinforcement Period','Time of the Start of the ITI','Licks Left','Licks Right','Total Licks','Correct','Wrong','Reward Delivered','Punishment Delivered','Probe Trial'};

    % Save the iteration and trial matrices as csv files
    csvwrite([vr.outputPath vr.outputFile '_IterationData' vr.name_addendum vr.addl_name_addendum '.csv'],IterationData);
    csvwrite([vr.outputPath vr.outputFile '_TrialData' vr.name_addendum vr.addl_name_addendum '.csv'],TrialData);
    
    % Close everything
    close all;

