classdef DynamicsExample < wl_experiment_v1_6
    properties
        ObjectPosition = zeros(3,1);
        ObjectColor = [1 0 0 0.5];
        ObjectRadius = 0;
        
        HomeActiveColor = [0.8 0.8 0.8 0.7];
        HomeInactiveColor = [0.7 0.7 0.7 0.2];
        HomeColor = [0.8 0.8 0.8 0.7];
        
        MovementDurationTimeoutMessage = '';
        
        Explosion = [];
    end
    methods
        % must implement ALL abstract methods or matlab will complain.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function run(WL, varargin)
            try
                WL.GUI = wl_gui('DynamicsExample', 'test', 'cfg_DynamicsExample', 'OUT3ISO', varargin{:});
                % These are custom parameters specific for this experiment.
                WL.GUI.addParam('array', 'trials_to_run', []);
                WL.GUI.addParam('numeric', 'reload_table', 0);
                
                ok = WL.initialise();
                if ~ok
                    WL.printf('Initialisation aborted\n')
                    return;
                end
                WL.my_initialise();

                % wl_robot should check mouse flag so this would be a single line
                % what about paradigms that don't use robot (e.g. eye tracker)
                if ~WL.cfg.MouseFlag
                    WL.Robot = WL.robot(WL.cfg.RobotName);
                else
                    WL.Robot = WL.mouse(WL.cfg.RobotName);
                end
                
                % This new line required to use the Liberty.
                %WL.Liberty = wl_liberty('LIBERTY.CFG'); % 'LIBERTY.CFG' is a text configuration file.
                %WL.Hardware = wl_hardware(WL.Robot,WL.Liberty); % The Liberty object must be included here in the hardware list.
                
                WL.Hardware = wl_hardware(WL.Robot); 
                
                ok = WL.Hardware.Start();
                if( ok ) % This will eventually happen inside wl_robot
                    ok = WL.Robot.ForceMaxSet(WL.cfg.RobotForceMax);
                end
                
                if ok
                    %WL.test_timings(0.5);
                    WL.main_loop();
                end
                
                WL.Hardware.Stop();
                
                clear mex
                
            catch msg
                WL.close(msg); % Does everything we need to do before exiting.
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function my_initialise(WL, varargin)
            WL.state_init('INITIALIZE','SETUP','HOME','OBJECT','START','DELAY','GO','MOVEWAIT',...
                'MOVING','FINISH','NEXT','INTERTRIAL','EXIT','TIMEOUT','ERROR','REST');
            % Stimulus frame counter.
            WL.FrameCounter.Stimulus = wl_frame_counter();
            WL.Timer.Stimulus = wl_timer;
            WL.Timer.ObjectPulseTimer = wl_timer;
            WL.ObjectColor = WL.cfg.ObjectActiveColor;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function idle_func(WL)
            ok = WL.Hardware.GetLatest(WL);
            WL.state_process();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function keyboard_func(WL,keyname)
            %WL.printf('Key pressed: %s\n',keyname);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function display_func(WL)
            %fcount = WL.FrameCounter.Stimulus.GetCount();
            %dt = WL.Timer.dt.Toc;
            %WL.Timer.dt.Tic;
            
            Screen('BeginOpenGL', WL.Screen.window);
            
            if WL.State.Current >= WL.State.START && WL.State.Current <= WL.State.MOVING
                WL.ObjectPosition = WL.Robot.Position;
            end

            if WL.State.Current > WL.State.OBJECT && WL.State.Current <= WL.State.MOVING
                pT = min(1, WL.Timer.ObjectPulseTimer.GetTime/WL.cfg.ObjectPulseDuration);
                % Pulse the size of the object by -20% to signal pick-up
                WL.ObjectRadius = WL.Trial.ObjectRadius*(1+0.1*(cos(2*pi*pT)-1)); % NB (cos(x)-1) is in [0, -2]
            end
            
            % during Probe/channel trials cursor disappears at start of the trial and is not visible on return trials
            WL.GW.display_cursor = WL.State.Current >= WL.State.SETUP && WL.State.Current <= WL.State.INTERTRIAL;
            WL.GW.display_target = WL.State.Current >= WL.State.OBJECT;% && WL.State.Current <= WL.State.FINISH;
            WL.GW.display_object = WL.State.Current >= WL.State.OBJECT;% && WL.State.Current <= WL.State.FINISH;
            
            if WL.GW.display_cursor
                % Draw the cursor
                wl_draw_sphere(WL.Robot.Position, WL.cfg.CursorRadius, [0 1 1], 'Alpha', 1)
                % Draw the targets
                for ti = 1:size(WL.cfg.TargetPositions)
                    tmp = WL.cfg.TargetPositions(ti,:);
                    % Draw the current target brightly after the object has been grabbed
                    if WL.GW.display_target && all(tmp==WL.Trial.TargetPosition')
                        wl_draw_sphere(WL.Trial.TargetPosition, WL.cfg.TargetRadius, [1 1 0.0], 'Alpha', 1);
                    else % draw the inactive targets dimly
                        wl_draw_sphere(WL.cfg.TargetPositions(ti,:), WL.cfg.TargetRadius, [0.6 0.6 0.0], 'Alpha', 0.2);
                    end
                end
            end
            wl_draw_sphere(WL.cfg.HomePosition, WL.cfg.HomeRadius, WL.HomeColor(1:3), 'Alpha', WL.HomeColor(4))
            
            if WL.GW.display_object
                if ~isempty(WL.Explosion) && WL.State.Current > WL.State.MOVING
                    WL.Explosion.ExplodeProcess(WL);
                else
                    wl_draw_sphere(WL.ObjectPosition, WL.ObjectRadius, WL.ObjectColor(1:3), 'Alpha', WL.ObjectColor(4));
                end
            end
            
            Screen('EndOpenGL',WL.Screen.window)
            txt = sprintf('Trial = %i State = %s', WL.TrialNumber, WL.State.Name{WL.State.Current});
            
            WL.draw_text(txt, 'center', 100, 'col', [1 1 1]);
            WL.draw_text(txt, 20, WL.Screen.windowRect(4)-30, 'fontsize', 12, 'fliph', 0,'flipv', 0, 'col', [1 1 1]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function state_process(WL)
            if  WL.GW.TrialRunning && any(~WL.Robot.Active) % If robot is not active, abort current trial.
                WL.trial_abort('Handle Switch',WL.State.SETUP);
            end
            
            switch WL.State.Current % State processing.
                
                case WL.State.INITIALIZE % Initialization state.
                    WL.Timer.Paradigm.ExperimentTimer.Reset;
                    WL.state_next(WL.State.SETUP);
                    
                case WL.State.SETUP % Setup details of next trial, but only when robot stationary and active.
                    if all(WL.Robot.Active) %WL.robot_stationary() && all(WL.Robot.Active)
                        WL.trial_setup();
                        WL.HomeColor = WL.HomeActiveColor;
                        WL.state_next(WL.State.HOME);
                    end
                    
                case WL.State.HOME % Subject returns robot to home position (and stationary and active).
                    if (WL.robot_home() && all(WL.Robot.Active)) || (WL.Trial.FieldType == 3) %(WL.robot_stationary() &&  WL.robot_home() && all(WL.Robot.Active)) || (WL.Trial.FieldType == 3)
                        WL.HomeColor = WL.HomeInactiveColor;
                        %WL.ObjectColor = WL.cfg.ObjectInactiveColor;
                        % object is initialized below the cursor to make it easier to see
                        WL.ObjectPosition = WL.cfg.ObjectHomePosition + [0 0 -WL.Trial.ObjectRadius]';
                        WL.ObjectRadius = WL.Trial.ObjectRadius;
                        WL.FrameCounter.Stimulus.ResetCount();
                        WL.FrameCounter.Stimulus.StartCount();
                        WL.Timer.Stimulus.Tic();
                        WL.Timer.StimulusTime.Reset();
                        WL.state_next(WL.State.OBJECT);
                    end
                    
                case WL.State.OBJECT % Subject moves to the object under null field to "pick it up"
                    if (WL.robot_stationary() && WL.robot_at_object() && all(WL.Robot.Active)) || (WL.Trial.FieldType == 3)
                        WL.trial_start();
                        WL.play_sound(WL.cfg.pickbeep);
                        %WL.ObjectColor = WL.cfg.ObjectActiveColor;
                        WL.Timer.MovementReactionTimer.Reset();
                        WL.Timer.ObjectPulseTimer.Reset();
                        WL.state_next(WL.State.MOVEWAIT);
                        %WL.state_next(WL.State.START);
                    end
                    
%                 case WL.State.START % Start trial.
%                     WL.trial_start();
%                     if WL.Trial.FieldType == 3
%                         WL.state_next(WL.State.MOVING);
%                     else
%                         WL.state_next(WL.State.DELAY);
%                     end
%                     
%                 case WL.State.DELAY % Delay period before go signal.
%                     if WL.State.Timer.GetTime > WL.cfg.TrialDelay
%                         WL.FrameCounter.Stimulus.ResetCount();
%                         WL.FrameCounter.Stimulus.StartCount();
%                         WL.Timer.Stimulus.Tic();
%                         WL.state_next(WL.State.GO);
%                     elseif  WL.movement_started()
%                         WL.trial_abort('Moved Too Soon', WL.State.SETUP);
%                     end
%                     
%                 case WL.State.GO % Go signal to cue movement.
%                     WL.Timer.MovementReactionTimer.Reset();
%                     WL.play_sound(WL.cfg.highbeep);
%                     WL.Timer.StimulusTime.Reset();
%                     WL.state_next(WL.State.MOVEWAIT);
                    
                case WL.State.MOVEWAIT
                    if  WL.movement_started()
                        WL.Timer.MovementDurationTimer.Reset();
                        WL.Trial.MovementReactionTime = WL.Timer.MovementReactionTimer.GetTime;
                        WL.state_next(WL.State.MOVING);
                    elseif WL.Timer.MovementReactionTimer.GetTime>WL.cfg.MovementReactionTimeOut
                        WL.state_next(WL.State.TIMEOUT);
                    end
                    
                case WL.State.MOVING
                    if WL.movement_finished(WL.cfg.StopAtTarget)
                        ok = WL.Robot.RampDown();
                        mt = WL.Timer.MovementDurationTimer.GetTime;
                        if mt < WL.cfg.MovementTooFastTimeout
                            WL.MovementDurationTimeoutMessage = 'Too Fast';
                            WL.state_next(WL.State.TIMEOUT);
                            return
                            
                        elseif mt < WL.cfg.MovementTooFastWarning
                            WL.play_sound(WL.cfg.fastwarnbeep);
                        elseif mt > WL.cfg.MovementTooSlowWarning
                            WL.play_sound(WL.cfg.slowwarnbeep);
                        else
                            WL.play_sound(WL.cfg.placebeep);
                        end
                        
                        if ~isempty(WL.Explosion)
                            WL.Explosion.ExplodePop(WL.Robot.Position);
                        else
                            WL.ObjectPosition = WL.Robot.Position; % Update for the last time
                        end
                        WL.HomeColor = WL.HomeActiveColor;
                        %WL.ObjectColor = WL.cfg.ObjectInactiveColor;
                        WL.Trial.MovementDurationTime = WL.Timer.MovementDurationTimer.GetTime();
                        if( WL.Trial.FieldType ~= 3 ) % Not a PMove trial
                            WL.Trial.StimulusFrameCount = WL.FrameCounter.Stimulus.GetCount();
                            WL.Trial.StimulusTime = WL.Timer.Stimulus.Toc();
                            %WL.printf('Stimulus on for %d frames, %.3f ms (%.1f Hz).\n',WL.Trial.StimulusFrameCount,WL.Trial.StimulusTime,WL.Trial.StimulusFrameCount/WL.Trial.StimulusTime);
                        end
                        WL.state_next(WL.State.FINISH);
                    elseif  WL.Timer.MovementDurationTimer.GetTime > WL.cfg.MovementTooSlowTimeout
                        WL.MovementDurationTimeoutMessage = 'Too Slow';
                        WL.state_next(WL.State.TIMEOUT);
                    end
                    
                case WL.State.FINISH
                    if WL.State.Timer.GetTime > WL.cfg.FinishDelay % Trial has finished so stop trial.
                        %ok = WL.Robot.RampDown();
                        WL.trial_stop();
                        WL.Timer.Paradigm.InterTrialDelayTimer.Reset;
                        
                        if ~WL.trial_save()
                            WL.printf('Cannot save Trial %d.\n',WL.TrialNumber);
                            WL.state_next(WL.State.EXIT);
                        else
                            % non fatal error on too-slow trials.
                            %here we have changed FieldType ~=3  so that on return trials we do not get the 'too
                            %slow' messages
%                             if  WL.Trial.MovementDurationTime >= WL.cfg.MovementDurationTimeOut && (WL.Trial.FieldType ~= 3)
%                                 WL.error_state( 'Moved Too Slow',WL.State.NEXT);
%                             else
                            WL.state_next(WL.State.NEXT);
%                             end
                        end
                    end
                    
                case WL.State.NEXT
                    if WL.Trial.RestFlag==1
                        WL.state_next(WL.State.REST);
                    elseif  ~WL.trial_next()
                        WL.state_next(WL.State.EXIT);
                    else
                        WL.state_next(WL.State.INTERTRIAL);
                    end
                    
                case WL.State.INTERTRIAL % Wait for the intertrial delay to expire.
                    if WL. Timer.Paradigm.InterTrialDelayTimer.GetTime > WL.cfg.InterTrialDelay
                        WL.state_next(WL.State.SETUP);
                    end
                    
                case WL.State.EXIT
                    WL.GW.ExperimentSeconds = WL.Timer.Paradigm.ExperimentTimer.GetTime;
                    WL.GW.ExperimentMinutes = WL.GW.ExperimentSeconds / 60.0;
                    WL.printf('Game Over (%.1f minutes)',WL.GW.ExperimentMinutes);
                    WL.GW.ExitFlag = true;
                    
                case WL.State.TIMEOUT
                    switch WL.State.Last % Which state had the timeout?
                        case WL.State.MOVEWAIT
                            WL.trial_abort('Move After Beep',WL.State.SETUP);
                        case WL.State.MOVING
                            WL.trial_abort(WL.MovementDurationTimeoutMessage,WL.State.SETUP);
                        otherwise
                            WL.trial_abort(sprintf('%s TimeOut',WL.State.Name{WL.State.Last}),WL.State.SETUP);
                    end
                    
                case WL.State.ERROR
                    if  WL.State.Timer.GetTime > WL.cfg.ErrorWait
                        WL.error_resume();
                    end
                    
                case WL.State.REST
                    RestBreakRemainSeconds = (WL.cfg.RestBreakSeconds -  WL.State.Timer.GetTime);
                    WL.cfg.RestBreakRemainPercent = (RestBreakRemainSeconds / WL.cfg.RestBreakSeconds);
                    
                    if  RestBreakRemainSeconds < 0
                        WL.Trial.RestFlag = 0;
                        WL.state_next(WL.State.NEXT);
                    end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function trial_start(WL)
            WL.Timer.Paradigm.TrialTimer.Reset();
            
            WL.Trial.StimulusFrameCount = 0;
            WL.Trial.StimulusTime = 0;
            
            ok = WL.Hardware.DataStart(); % If error in datastart, throw error.
            if( ~ok )
                error('WL.Hardware.DataStart() Failed.');
            end
            
            switch( WL.Trial.FieldType )
                case 0 % Null field
                    ok = WL.Robot.FieldNull();
                case 1 % Exposure, viscous curl field
                    ViscousMatrix = WL.Trial.FieldConstants(1) * ...
                        [ 0 -1  0; ...
                          1  0  0; ...
                          0  0  0 ];
                    ok = WL.Robot.FieldViscous(ViscousMatrix);
                case 2 % Channel
                    ok = WL.Robot.FieldChannel(WL.Trial.TargetPosition,WL.Trial.FieldConstants(1),WL.Trial.FieldConstants(2));
                case 3 % Passive return movement
                    ok = WL.Robot.FieldPMove(WL.Trial.TargetPosition,0.5,0.2);
            end
            
            %WL.printf('TrialStart() Trial=%d, Field=%d\n',WL.TrialNumber,WL.Trial.FieldType);
            %WL.printf('RobotField=%d, Started=%d\n',WL.Trial.FieldType,ok);
            
            WL.Explosion = wl_explode(WL.Trial.ObjectRadius,WL.cfg.ObjectActiveColor,sqrt(WL.Trial.ObjectRadius^2/16));
            
            ok = WL.Robot.RampUp();
            WL.GW.TrialRunning = true;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function miss_trial(WL,MissTrialType)
         % this is only used if you want to change the order of trials after a miss-trial  
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out=robot_home(WL)
            err = norm(WL.Robot.Position(1:2) - WL.cfg.HomePosition(1:2));
            out = err<WL.cfg.HomeTolerance;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out=robot_at_object(WL)
            err = norm(WL.Robot.Position(1:2) - WL.cfg.ObjectHomePosition(1:2));
            out = err<WL.cfg.ObjectTolerance;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = movement_finished(WL,StopAtTarget)
            err = norm(WL.Robot.Position(1:2) - WL.Trial.TargetPosition(1:2));
            if (StopAtTarget)
                out = (err < WL.cfg.TargetTolerance) && WL.robot_stationary();
            else
                out = (err < WL.cfg.TargetTolerance);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = movement_started(WL)
            err = norm(WL.Robot.Position(1:2) - WL.cfg.ObjectHomePosition(1:2));
            out = err > 0.5;
            %out = ~WL.robot_at_object();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end





