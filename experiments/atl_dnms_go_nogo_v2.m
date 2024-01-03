function code = Tunel
% Tunel   Code for the ViRMEn experiment Tunel.
%   code = Tunel   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT

% Deliver after stimuli are presented, withdraw if incorrect choice

% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)

vr.name = 'M11'; % DO NOT FORGET TO CHANGE THE NAME BEFORE EXPERIMENT AND THE DATA FILE NAME AT THE END !
vr.cagenum = 'C249376'; % Enter cage number
% C249376 M11
% C258735 F20 (death) F21 (out) F22
% C266941 F21(R) F22(R) F25
% C266148 F28 (out) F29(grabs) F30(not licking)

vr.punishment = false; % if this is set to true than some kind of punishment will be present the type of the punishment will be determined later in the decision period code

vr.puff_punishment = false; %deliver air puff for nogo licks

vr.correct_left_bias = false; % This corrects left bias by overrepresenting the trials which result in reward coming from the right lickport (i.e. will run the experiment where 80% of the trials are on the right and 20% of the trials are on the left. Setting this to true makes sense only for randomly interleaved trials 
    
vr.correct_right_bias = true; % If this is set to true it will do vice versa compared to if previous is set to true. Setting this to true makes sense only for randomly interleaved trials a

% Set the path for the experiment and data .mat files to be saved 
vr.finalPathname = 'C:\Users\VR Rig #2\Documents\MATLAB\Andrew_data\';

% create a string for the name of the file by adding vr.name to the current
% path 
vr.filename = datestr(now,'mmddyyyy');
vr.finalname = strcat(vr.name,'_',vr.cagenum,'_Date_',vr.filename);
vr.str_punish = '_with_punishment';

% save experiment .mat file
exper = vr.exper;

% Make a different name by adding with_punishment to file name if
% punishment is present
if vr.punishment == true
    save([vr.finalPathname vr.finalname vr.str_punish '_refresher.mat'],'exper');
else
    save([vr.finalPathname vr.finalname '_refresher.mat'],'exper');
end


% Start analog session to recieve inputs from two lickports
vr.analogSess = daq.createSession('ni');
addAnalogInputChannel(vr.analogSess,'Dev1',[0 1], 'Voltage'); %device object, device ID, channel numbers, measurement type (voltage)

% Start digital session to control stimulus delivery
vr.DigitalSess = daq.createSession('ni');
addDigitalChannel(vr.DigitalSess,'dev1', 'Port0/Line0:7', 'OutputOnly');
addDigitalChannel(vr.DigitalSess,'dev1', 'Port1/Line0:1', 'OutputOnly');
outputSingleScan(vr.DigitalSess,[0 0 0 0 0 0 0 0 0 0]);

% Define triggering variables for tone,odor,air and the start of new trial
vr.newtrial = 0;
vr.tone = 0;
vr.Odor1_valve = 0;
vr.Odor2_valve = 0;
vr.AirON = 0;
vr.Extend_actuator = 0;
vr.Retract_actuator = 0;
vr.Deliver_water_left = 0;
vr.Deliver_water_right = 0;
vr.synchPulse = 0;
vr.puffPun = 0;

% Define variables which enable execution of proper trial chronology
vr.tone_end = 0;
vr.odor1_end = 0;
vr.odor2_end = 0;
vr.ITI_happened = 0;
vr.odor2_begin = 0;
vr.decision = 0;
vr.decision_end = 0;
vr.decision_runing = 0;
vr.wrong_comitment = 0;
vr.wrong = 0;
vr.correct = 0;
vr.wrong_picked_left = 0;
vr.wrong_picked_right = 0;
vr.reward_runing = 0;
vr.trialset = 0;
vr.stim1_first = 0;
vr.stim2_first = 0;

% Define variables which keep track of the trial chronology
vr.move = 0;
vr.counter = 0;
vr.trialnum_true = 1;
vr.trialnum = 1;

% Define variables to store analog signals from the lickports
vr.decision_time = [];
vr.licks_right = 0;
vr.licks_left = 0;
vr.record_lick = [0 0];
vr.decision_label = 0;
vr.nogo_licks = 0;
vr.go_licks = 0;

% Define number of trials, create the vector of different conditions
% of length which is eaqual to the number of trials where:
% Condition 1 ---> Odor1 - Delay - Odor2
% Condition 2 ---> Odor2 - Delay - Odor1
% Condition 3 ---> Odor1 - Delay - Odor1
% Condition 4 ---> Odor2 - Delay - Odor2
% And finally within the same block create a vector of delays whose length
% is eaqual to the number of trials by rounding an absolute value of a number
% taken from the Gaussian distribution with mean of vr.dist_mean and standard 
% deviation of vr.dist_std


vr.trials = 10; % Enter number of trials 

% % Create experiment structure if one wants to use blocks of 5 trials
% within experiment if one wants to randomly interleave trials this should
% be comented (alternatively one can write an if statement with boolean
% variable as for punishment above or biass correction for example)

% vr.resolution = 5;
% vr.j = 1;
% vr.conditions = [];
% vr.rand = [1 2 3 4];
% while vr.j <= vr.trials
% vr.rand = vr.rand(randperm(length(vr.rand)));
% vr.temp1 = ones(5,1)*vr.rand(1);
% vr.temp2 = ones(5,1)*vr.rand(2);
% vr.temp3 = ones(5,1)*vr.rand(3);
% vr.temp4 = ones(5,1)*vr.rand(4);
% vr.conditions = [vr.conditions; vr.temp1; vr.temp2; vr.temp3; vr.temp4];
% vr.j = vr.j + vr.resolution*4;
% end
% vr.conditions(vr.trials+1:(length(vr.conditions))) = [];

% If one wants to interleave 
vr.disl = 0;
vr.disr = 0;

% Check for bias correction and if true overrepresent trials. Otherwise
% this will just create randomly interleaved trials without
% overrepresentation
if vr.correct_left_bias == true
    vr.disl = 21;
    vr.disr = 0;
end

if vr.correct_right_bias == true
    vr.disl = 0;
    vr.disr = 21;
end

vr.delta = round(vr.trials/4);
% disp(vr.delta)
vr.one = ones(vr.delta-vr.disr+vr.disl,1)*1; % For different combinations refer to 
% coments above and multiply with according numbers
vr.two = ones(vr.delta-vr.disr+vr.disl,1)*2;
vr.three = ones(vr.delta-vr.disl+vr.disr,1)*3; 
vr.four = ones(vr.delta-vr.disl+vr.disr,1)*4;

vr.conditions = [vr.one; vr.two; vr.three; vr.four];
vr.conditions = vr.conditions(randperm(length(vr.conditions)));

vr.delays = ones(vr.trials,1)*1; % Create a vector of delays
% This lines create a probe trials for shaping without decision period such
% that 10% of the trials are probe trials i.e. trials with decision period
% if the first element of vector below is set to zeros(round(vr.trials - vr.trials*0.1),1)
vr.probetrials = [ones(round(vr.trials - vr.trials*0.1),1); ones(round(vr.trials*0.1),1)];
vr.probetrials = vr.probetrials(randperm(length(vr.probetrials)));

% disp(vr.delays);
% disp(vr.conditions);
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp('Experiment started');
disp('------------------------');
if (round(vr.iterations) == 100)
    vr.synchPulse=1; %synch pulse goes to col (:,17?) of data struct
    disp('output dout')
end    

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
% Annul the stimuli triggering variables at the start of each trial
vr.Odor1_valve = 0;
vr.Odor2_valve = 0;
vr.tone = 0;
vr.AirON = 0;
vr.record_lick = [0 0];
vr.Extend_actuator = 0;
vr.Retract_actuator = 0;
vr.wrong = 0;
vr.correct = 0;
vr.decision_label = 0;
vr.synchPulse = 0;
vr.wrong_picked_left = 0;
vr.wrong_picked_right = 0;
vr.trialset = 0;
vr.stim1_first = 0;
vr.stim2_first = 0;


% Send a tone to initiate a trial. The tone length presentation is
% goverened by the arduino code in which it has been set to 2000 ms
if vr.trialnum_true == 1;
    if (round(vr.iterations) == 100)
    vr.synchPulse=1; %synch pulse goes to col (:,17?) of data struct
    disp('output dout')
    end    
end

if (round(vr.iterations) == 100) || (vr.newtrial == 1)
    disp(strcat('Initiated trial: ',num2str(vr.trialnum_true)));
%    disp(strcat('Virtual trial: ',num2str(vr.trialnum)));
    vr.tone = 1;
    vr.newtrial = 0;
    vr.tone_end = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter; 

vr.record_lick=inputSingleScan(vr.analogSess);

% After 2 seconds have passed present the first odor depending on the above
% defined condition during the current trial. The length of the odor
% presentation is regulated by the arduino code as well and it is set to
% 5000 ms - 02/05/2021
if (round(vr.move,1) == 1.3) && (vr.tone_end == 1)
    switch vr.conditions(vr.trialnum)
    case {1,3}
        vr.Odor1_valve = 1;
        vr.stim1_first = 1;
        disp('Stimulus_1: Odor1');
    case {2,4}
        vr.Odor2_valve = 1;
        vr.stim2_first = 1;
        disp('Stimulus_1: Odor2');
    end
    vr.tone_end = 0;
    vr.odor1_end = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter;

% if (round(vr.move,1) == 1.3) && (vr.tone_end == 1) && ((vr.wrong_picked_left == 1) || (vr.wrong_picked_right == 1))
%     switch vr.conditions(vr.trialnum_true)
%     case {1,3}
%         vr.Odor1_valve = 1;
%         disp('Stimulus_1: Odor1');
%     case {2,4}
%         vr.Odor2_valve = 1;
%         disp('Stimulus_1: Odor2');
%     end
%     vr.tone_end = 0;
%     vr.odor1_end = 1;
%     vr.counter = vr.timeElapsed;
%     vr.wrong_picked_left = 0;
%     vr.wrong_picked_right = 0;
% end
% vr.move = vr.timeElapsed - vr.counter; 

% % After 2 seconds have passed display that the delay has started and for
% % how long will this delay be in this particular trial
% % 1.9s is default (0.5s delay)
if (round(vr.move,1) == 3.4) && (vr.odor1_end == 1)
    vr.string = strcat('Delay: ',num2str(vr.delays(vr.trialnum)));
    vr.string = strcat(vr.string,' s');
    disp(vr.string);
    vr.odor1_end = 0;
    vr.odor2_begin = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter; 

% After the number of seconds amounting to the delay time has passed
% present the second odor depending on the condition for this particular 
% trial. The length of the odor presentation is regulated by the arduino 
% code as well and it is set to 5000 ms - 02/05/2021
if (round(vr.move,1) == vr.delays(vr.trialnum)) && (vr.odor2_begin == 1)
    switch vr.conditions(vr.trialnum)
    case {1,4}
        vr.Odor2_valve = 1;
        disp('Stimulus_2: Odor2');
    case {2,3}
        vr.Odor1_valve = 1;
        disp('Stimulus_2: Odor1');
    end
    vr.odor2_begin = 0;
    vr.odor2_end = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter;

% Initiation of the decision point. Extend the actuator display time
% elapsed and check for the presence of the probe trials
if (round(vr.move,1) == 0.9) && (vr.odor2_end == 1)
    vr.Extend_actuator = 1;
    string = num2str(round(vr.timeElapsed,1));
    disp(strcat('Decision Point Started at: ', string));
    vr.odor2_end = 0;
    vr.decision_end = 1;
    vr.decision_runing = 1;
    vr.decision_label = 1;
    if vr.probetrials(vr.trialnum) == 0
        switch vr.conditions(vr.trialnum)
            case {3,4}
                vr.Deliver_water_left = 1;
            case {1,2}
                vr.Deliver_water_left = 0;
        end
    else
%         disp('Probe trial');
    end
    vr.counter = vr.timeElapsed+1; %1s delay once extended to make decision
end
vr.move = vr.timeElapsed - vr.counter;

% 
% Count the licks on the left and right lickports during the decision
% period, check for probe trials and act acordingly 
if (vr.decision_runing == 1)
    vr.decision_label = 1;
    if vr.probetrials(vr.trialnum) == 0
        switch vr.conditions(vr.trialnum)
            case {3,4}
                vr.Deliver_water_left = 1;
            case {1,2}
                vr.Deliver_water_left = 0;
        end
    end
   
    %vr.record_lick=inputSingleScan(vr.analogSess);
    if vr.record_lick(1) > 3
        vr.licks_left = vr.licks_left + 1;
    end
    
    if vr.record_lick(2) > 3
        vr.licks_right = vr.licks_right + 1;
    end
    vr.decision_time = [vr.decision_time vr.dt];
end
% 
% End of the decision point. Check for the presence of the punishment and
% act acordingly depending on the combination of odors being present within
% the experiment
if (round(vr.move,1) == 2.8) && (vr.decision_end == 1)
    vr.decision_label = 1;
    string1 = 'No. of licks: ';
%    string2 = 'No. of licks right: ';   
    vr.L(vr.trialnum,1)=vr.licks_left;
    vr.L(vr.trialnum,2)=vr.licks_right;
    disp(strcat(string1,num2str(vr.licks_left)));
%    disp(strcat(string2,num2str(vr.licks_right)));
%     disp(mean(vr.decision_time));
    vr.decision_time = [];
    
    switch vr.conditions(vr.trialnum)
        case {3,4}
           vr.go_licks = vr.go_licks + vr.licks_left;
           if vr.licks_left <= vr.licks_right
               vr.wrong_picked_right = 1;
%                disp('Wrong choice');
                vr.wrong = 1;
            if vr.punishment == true
                vr.wrong_comitment = 1; % put 1 here if wanting to retract
                % actuator as a punishment opposing to just not giving water
                % Also make sure to match the times with the arduino
            else
                vr.Deliver_water_left = 1;
                vr.correct = 1;
                disp('Delivering water to left lickport');
            end
           else 
%                   disp('Correct choice');
                if vr.punishment == true
                   vr.correct = 1;
                   vr.Deliver_water_left = 1;
                   disp('Delivering water to left lickport');
                else
                    vr.Deliver_water_left = 1;
                    vr.correct = 1;
                    disp('Delivering water to left lickport');
                end
            end
        case {1,2}
           vr.nogo_licks = vr.nogo_licks + vr.licks_left;
           if vr.licks_left > vr.licks_right
               vr.wrong_picked_left = 1;
%                disp('Wrong choice');
                vr.wrong = 1;
            if vr.puff_punishment == true
                vr.wrong_comitment = 1; % Same here
                vr.puffPun = 1;
            else
%                vr.Deliver_water_right = 1;
                vr.correct = 1;
                disp('Delivering water to right lickport');
            end
           else 
 %                  disp('Correct choice');
               if vr.punishment == true
                   vr.correct = 1;
 %                  vr.Deliver_water_right = 1;
                else
%                    vr.Deliver_water_right = 1;
                    vr.correct = 1;
%                    disp('Delivering water to right lickport');
                end
           end
     end 
    vr.decision_end = 0;
    vr.decision_runing = 0;
    vr.decision = 1;
    vr.reward_runing = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter;

% if (vr.reward_runing == 1)
%     vr.record_lick=inputSingleScan(vr.analogSess);
% end

% When the decision time has passed and the decision has been made being
% the intertrial perion (ITI), display that it has started and turn of the
% air. The duration of air being of is regulated by the arduino code which
% is set to 5000 ms - 02/05/2021
if ((vr.decision == 1) && (round(vr.move,1) == 4)) || (vr.wrong_comitment == 1)
    vr.AirON = 1;
    vr.decision_runing = 0;
    vr.Deliver_water_left = 0;
    vr.Deliver_water_right = 0;
    vr.licks_left = 0;
    vr.licks_right = 0;
    vr.Retract_actuator = 1;
    %vr.synchPulse = 1;
    
    vr.odor1_end = 0;
    string = num2str(round(vr.timeElapsed,1));
    disp(strcat('Decision Point Ended at: ', string));
    disp('ITI started');
    vr.decision = 0;
    vr.wrong_comitment = 0;
    vr.ITI_happened = 1;
    vr.counter = vr.timeElapsed;
end
vr.move = vr.timeElapsed - vr.counter;

% After 15 seconds ITI ends
if ((round(vr.move,1) == 5) && (vr.ITI_happened == 1)) %&& (vr.wrong_picked_left == 0) && (vr.wrong_picked_right == 0)
    disp('ITI ended');
    disp('------------------------------');
    if ((vr.wrong == 1))
%         disp('will repeat trial');
         vr.trialnum = vr.trialnum;
     else
%         disp('will move on');
         vr.trialnum = vr.trialnum + 1;
     end 
    vr.newtrial = 1;
    vr.ITI_happened = 0;
    vr.trialnum_true = vr.trialnum_true +1;
    vr.puffPun=0;
end
vr.move = vr.timeElapsed - vr.counter;


% % After 10 seconds have passed display that the intertrial period has ended
% % and begin new trial - if choice was incorrect
% if ((round(vr.move,1) == 5) && (vr.trialset == 1))
%     disp('ITI ended');
%     disp('------------------------------');
%     vr.L(vr.trialnum,3)=vr.conditions(vr.trialnum);
%     vr.L(vr.trialnum,4)=vr.delays(vr.trialnum);
%     vr.L(vr.trialnum,5)=vr.probetrials(vr.trialnum);
%     vr.newtrial = 1;
%     vr.trialnum_true = vr.trialnum_true + 1;
%     vr.trialset = 0;
%     vr.trialnum = vr.trialnum;
%     vr.wrong_picked_left = 0;
%     vr.wrong_picked_right = 0;
%     %vr.synchPulse = 0;
% end
%%%%%
% If all the trials have been executed end the experiment and display that
% the experiment has ended.
if (vr.trialnum_true == vr.trials + 1)
    vr.experimentEnded = true;
    disp('Experiment ended');
    string = num2str(round(vr.nogo_licks,1));
    disp(strcat('No-go licks: ', string));
    string = num2str(round(vr.go_licks,1));
    disp(strcat('Go licks: ', string));
end

% At the end of each iteration output a digital signal to NI-DAQ to trigger
% the stimuli which were defined to be triggered in the previous if
% statements by setting them to 1. Otherwise if 0 they will not be
% triggered
outputSingleScan(vr.DigitalSess,[vr.AirON vr.tone vr.Odor1_valve vr.Odor2_valve vr.Deliver_water_left vr.Deliver_water_right vr.Extend_actuator vr.Retract_actuator vr.synchPulse vr.puffPun]);

% Save all the relevant variables to a matrix called vr.M and save all of
% the licking data to the matrix called vr.D
vr.M(vr.iterations,1)=vr.timeElapsed;
vr.M(vr.iterations,2)=vr.AirON;
vr.M(vr.iterations,3)=vr.tone;
vr.M(vr.iterations,4)=vr.Odor1_valve;
vr.M(vr.iterations,5)=vr.Odor2_valve;
vr.M(vr.iterations,6)=vr.newtrial;
vr.M(vr.iterations,7)=vr.Extend_actuator;
vr.M(vr.iterations,8)=vr.Retract_actuator;
vr.M(vr.iterations,9)=vr.Deliver_water_left;
vr.M(vr.iterations,10)=vr.Deliver_water_right;
vr.M(vr.iterations,11)=vr.trialnum_true;
vr.M(vr.iterations,12)=vr.correct;
vr.M(vr.iterations,13)=vr.wrong;
vr.M(vr.iterations,14)=vr.decision_label;
vr.M(vr.iterations,15)= vr.record_lick(1); %licks left
vr.M(vr.iterations,16)= vr.record_lick(2); %licks right
vr.M(vr.iterations,17)= vr.synchPulse;
vr.M(vr.iterations,18)= vr.go_licks;
vr.M(vr.iterations,19)= vr.nogo_licks;
vr.M(vr.iterations,20)= vr.stim1_first;
vr.M(vr.iterations,21)= vr.stim2_first;
vr.D(vr.iterations,:)=vr.record_lick;


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)

% Stop all outputs to NI-DAQ by outputing 0 on every channel
outputSingleScan(vr.DigitalSess,[0 0 0 0 0 0 0 0 0 0]);


% Terminate both analog and digital sessions in NI-DAQ
vr.analogSess.stop
vr.DigitalSess.stop;

% Make a struct out of the matrix containing variables (vr.M) and the
% matrix containing data (vr.D) and save it to a separate file
vr.note = '_DNMS_gonogo';
vr.filename1 = [vr.finalPathname vr.finalname vr.note];
%exper = struct('M',vr.M,'Data',vr.D,'Licks',vr.L);
exper = struct('M',vr.M,'Data',vr.D);

if vr.punishment == true
    save([vr.filename1 '_Data' vr.str_punish '_refresher.mat'],'exper');
else
    save([vr.filename1 '_Data' '_refresher.mat'],'exper');
end

% Close everything
close all;

