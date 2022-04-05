function ok = wl_initialise(WL)
%WL_INITIALISE sets up the I\O for the the experiments
%   WL_INITIALISE(OBJ) where OBJ is a pointer to the main experiment class.
%   Note the passing of OBJ is handled inplictly by the ExperimentBase
%   class thus you do not need to wory about passing this parameter.
%
%   Calling WL.initialise inside your experiment class
%   equates to wl_initialise(WL) where WL is the
%   first argument of every method in the class also known as an internal
%   pointer to the WLect. this and WL are more commonly used.
%
%   Usage (within instance of ExperimenBase):
%       WL.initialise


%  random seed
WL.GW.rng = rng('shuffle');

%%%%%%%%%%%%%% GUI processsing %%%%%%%%%%%%%%
%Initialise GUI
ok = ~WL.GUI.process_params();
if ~ok % Abrupt termination(closing) of GUI
    fprintf('GUI Terminated abruptly\n');
    return;
end

% Read params from GUI
WL.GW = WL.GUI.pSet;
delete(WL.GUI.UIFigure);

%%%%%%%%%%%%%% Parse config filep %%%%%%%%%%%%%%

wl_cfg_defaults(WL);

feval(WL.GW.config_file, WL, WL.GW.config_type);
WL.TrialDataClean=WL.TrialData;


%%%%%%%%%%%%%% debugging diary %%%%%%%%%%%%%%

if WL.cfg.Debug    
    WL.GW.log_file = ['_diary_log_' strrep(num2str(GetSecs), '.', '_') '.txt']; %log file
    diary(WL.GW.log_file); 
else
    %ListenChar(2);
    diary off
end

%%%%%%%%%%%%%% Timer set up %%%%%%%%%%%%%%

WL.Timer.Paradigm.ExperimentTimer = wl_timer;  %total experiment time
WL.Timer.Paradigm.TrialTimer = wl_timer; % from the start of a trial
WL.Timer.Paradigm.InterTrialDelayTimer = wl_timer; %from the end of the trial
WL.State.Timer = wl_timer;
WL.Timer.System.TrialSave = wl_timer;

N=10000;
B=20; %burn in i.e. counts to skip
WL.Timer.Graphics.idle_loop = wl_timer(N,B); %timer of the idle function
WL.Timer.Graphics.idle_func_2_idle_func = wl_timer(N,B); %timer of the idle function
WL.Timer.Graphics.idle_func = wl_timer(N,B); %timer of the idle function
WL.Timer.Graphics.flip_check = wl_timer(N,B);

WL.Timer.Graphics.flip_2_flip=wl_timer(N,B);
WL.Timer.Graphics.flip_2_display_func=wl_timer(N,B);
WL.Timer.Graphics.flip_request_2_flip=wl_timer(N,B);
WL.Timer.Graphics.display_func = wl_timer(400,B);
WL.Timer.Graphics.flip_request = wl_timer(N,B);

WL.Timer.MovementDurationTimer = wl_timer;
WL.Timer.MovementReactionTimer = wl_timer;
WL.Timer.StimulusTime = wl_timer;

WL.FrameCounter.Stimulus = wl_frame_counter;

%%%%%%%%%%%%%% Parameter set up %%%%%%%%%%%%%%

WL.GW.error_msg = '';
WL.GW.MissTrial = 0;
WL.GW.TrialRunning = false;
WL.GW.ExitFlag = false;

%%%%%%%%%%%%%% Audio set up %%%%%%%%%%%%%%

InitializePsychSound(1); % reallyneedlowlatency = 1

dev = PsychPortAudio('GetDevices');

WL.GW.AudioDevice = nan;

% First, find first audio device that is an ouput.
for k=1:length(dev)
    if( dev(k).NrOutputChannels > 0 )
        WL.GW.AudioDevice = dev(k).DeviceIndex;
        WL.GW.AudioDeviceIndex = k;
        WL.GW.AudioDeviceName = dev(k).DeviceName;
        break;
    end
end

% Next, if an AISO audio output device is present, use that instead.
for k=1:length(dev)
    if( (length(strfind(dev(k).DeviceName,'ASIO')) == 1) && (dev(k).NrOutputChannels > 0) )
        WL.GW.AudioDevice = dev(k).DeviceIndex;
        WL.GW.AudioDeviceIndex = k;
        WL.GW.AudioDeviceName = dev(k).DeviceName;
        break;
    end
end

if ~isnan(WL.GW.AudioDevice ) && ispc
    fprintf('Audio device: Name=%s ID=%d Index=%d\n',WL.GW.AudioDeviceName,WL.GW.AudioDevice,WL.GW.AudioDeviceIndex);
    WL.GW.AudioHandle = PsychPortAudio('Open', WL.GW.AudioDevice, 1, 1, 48000,2);
else
    WL.GW.AudioHandle = PsychPortAudio('Open');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PsychPortAudio('Verbosity' ,WL.cfg.verbose);
Screen('Preference', 'Verbosity', WL.cfg.verbose); % usuallly 3

%find all keyboards and create a queue for each
[WL.Keyboard.key_id, ~, ~] = GetKeyboardIndices();

for k=1:length(WL.Keyboard.key_id)
    KbQueueCreate(WL.Keyboard.key_id(k));
    KbQueueStart(WL.Keyboard.key_id(k));
end

WL.Keyboard.KeyNames = KbName('KeyNames');

%%%%%%%%%%%%%% File names for saving %%%%%%%%%%%%%%

tmp=strrep(strrep(strrep(datestr(now),':','-'),' ','-'),'-','_');

WL.GW.table_file = ['_table_' tmp '.mat'];  %temp file for table
WL.GW.save_file = [WL.GW.save_file '.mat'];  %save file
WL.GW.log_file = ['_log_' strrep(num2str(GetSecs), '.', '_') '.txt'];

if( ispc )
    % Check if T: drive exists, ig not save path is current directory.
    if( exist('T:\\.','dir') == 0 )
        WL.GW.save_path = '.\';
    else
        WL.GW.save_path = 'T:\';
    end
    %WL.GW.save_path = '.\Data\'; % EAC Added
    
    WL.GW.save_file = [ WL.GW.save_path WL.GW.save_file ];
    WL.GW.table_file = [ WL.GW.save_path WL.GW.table_file ];
    WL.GW.log_file = [ WL.GW.save_path WL.GW.log_file ];
end

%check so as to not overwrite file
if exist(WL.GW.save_file,'file')
    warning OFF BACKTRACE % prevent line number being displayed
    warning(sprintf('File %s exists. Overwrite? y/n [n]:',WL.GW.save_file));
    pause(0.001);
    
    %ListenChar(2)
    [ ~,keyname ] = WL.keyboard_read(true); % true = wait
    %ListenChar(1)
    
    if( keyname == 'Y' )
        warning('Overwriting file')
        delete(WL.GW.save_file);
    else
        ok=false;
        warning('Aborting to prevent overwriting file')
        warning ON BACKTRACE
        return
    end
    warning ON BACKTRACE
end

%%%%%%%%%%%%%% choose whether to reload TrialData %%%%%%%%%%%%%%

if WL.GW.reload_table
    D = load('TrialData');
    WL.TrialData = D.TrialData;
    clear D;
else
    TrialData = WL.TrialData;
    save('TrialData','TrialData');
end

%%%%%%%%%%%%%% choose trials to run %%%%%%%%%%%%%%
if length(WL.GW.trials_to_run) == 1
    WL.TrialData=WL.TrialData(WL.GW.trials_to_run:end,:);
elseif length(WL.GW.trials_to_run) == 2
    WL.TrialData=WL.TrialData(WL.GW.trials_to_run(1):WL.GW.trials_to_run(2),:);
end

%%%%%%%%%%%%%% start graphics %%%%%%%%%%%%%%

wl_start_screen(WL);





