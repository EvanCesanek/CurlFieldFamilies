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
        
        SkipToTest = false;
        
        BestStreak = 0;
        Streak = 0;
        StreakColor = [1 1 1];
        StreakPosnPx = [-100 -100];
        
        MPEBuffer = [];
        ScoreBuffer = [];
        Score = [];
        MissWarn = false;
        RobotPreviousPosition = [0 0];
    end
    methods
        % must implement ALL abstract methods or matlab will complain.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function run(WL, varargin)
            try
                WL.GUI = wl_gui('DynamicsExample', 'test', 'cfg_DynamicsExample', 'BASIC1', varargin{:});
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
                'MOVING','EXPLODE','FINISH','NEXT','INTERTRIAL','EXIT','TIMEOUT','ERROR','REST');
            % Stimulus frame counter.
            WL.FrameCounter.Stimulus = wl_frame_counter();
            WL.Timer.Stimulus = wl_timer;
            WL.Timer.ScoreFade = wl_timer;
            WL.Timer.ObjectPulseTimer = wl_timer;
            WL.ObjectColor = WL.cfg.ObjectActiveColor;
            WL.StreakPosnPx = WL.cm2pix(WL.cfg.HomePosition);
            
            WL.MPEBuffer = nan(1,WL.cfg.NumTrials);
            WL.ScoreBuffer = nan(1,WL.cfg.NumTrials);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function idle_func(WL)
            ok = WL.Hardware.GetLatest(WL);
            WL.state_process();
            WL.RobotPreviousPosition = WL.Robot.Position; %AFTER state_process
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function keyboard_func(WL,keyname)
            if strcmpi(keyname,'+')
                WL.SkipToTest = true;
            end
            %WL.printf('Key pressed: %s\n',keyname);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function display_func(WL)
            %fcount = WL.FrameCounter.Stimulus.GetCount();
            %dt = WL.Timer.dt.Toc;
            %WL.Timer.dt.Tic;
            
            Screen('BeginOpenGL', WL.Screen.window);
            
            if WL.State.Current >= WL.State.START && WL.State.Current <= WL.State.MOVING
                WL.ObjectPosition = WL.Robot.Position + WL.Trial.ObjectOffset;
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
                    elseif WL.cfg.DrawInactiveTargets % draw the inactive targets dimly
                        wl_draw_sphere(WL.cfg.TargetPositions(ti,:), WL.cfg.TargetRadius, [0.6 0.6 0.0], 'Alpha', 0.2);
                    end
                end
            end
            wl_draw_sphere(WL.cfg.HomePosition, WL.cfg.HomeRadius, WL.HomeColor(1:3), 'Alpha', WL.HomeColor(4))
            
            if WL.GW.display_object
                if ~isempty(WL.Explosion) && WL.State.Current > WL.State.MOVING
                    WL.Explosion.ExplodeProcess(WL);
                    if WL.MissWarn
                        wl_draw_sphere(WL.ObjectPosition, WL.ObjectRadius, WL.ObjectColor(1:3), 'Alpha', WL.ObjectColor(4));
                        wl_draw_rectangle_gl([WL.ObjectPosition(1:2) - WL.Trial.ObjectOffset(1:2); 0], -WL.Trial.ObjectOffset(1), WL.ObjectRadius, WL.Trial.TargetAngleDegreesCentered, WL.ObjectColor(1:3), [1 0.5], 'Alpha', 1);
                    end
                else
                    wl_draw_sphere(WL.ObjectPosition, WL.ObjectRadius, WL.ObjectColor(1:3), 'Alpha', WL.ObjectColor(4));
                    if strcmpi(WL.cfg.Feature,'location')
                        if WL.State.Current == WL.State.OBJECT
                            wl_draw_sphere(WL.cfg.ObjectHomePosition, WL.cfg.CursorRadius, [1 1 1], 'Alpha', 0.7)
                            wl_draw_rectangle_gl(WL.cfg.ObjectHomePosition + [0 0 -WL.Trial.ObjectRadius]', -WL.cfg.Offsets(WL.Trial.ObjectId), WL.ObjectRadius, WL.Trial.TargetAngleDegreesCentered, WL.ObjectColor(1:3), [1 0.5], 'Alpha', 1);
                        else
                            wl_draw_rectangle_gl(WL.Robot.Position + [0 0 -WL.Trial.ObjectRadius]', -WL.cfg.Offsets(WL.Trial.ObjectId), WL.ObjectRadius, WL.Trial.TargetAngleDegreesCentered, WL.ObjectColor(1:3), [1 0.5], 'Alpha', 1);
                        end
                    end
                end
            end
            
            if ~isempty(WL.Score)
                WL.draw_text(sprintf('+%i',WL.Score), WL.StreakPosnPx(1), WL.StreakPosnPx(2), 'col', [1 1 1 max(0, 1-WL.Timer.ScoreFade.GetTime())]);
            end
            WL.draw_text(sprintf('Score: %i', nansum(WL.ScoreBuffer)), WL.StreakPosnPx(1), WL.StreakPosnPx(2)+180, 'fontsize',16, 'col', WL.StreakColor);
            
            Screen('EndOpenGL',WL.Screen.window)
            txt = sprintf('Trial = %i State = %s', WL.TrialNumber, WL.State.Name{WL.State.Current});
            
            %WL.draw_text(txt, 'center', 100, 'col', [1 1 1]);
            WL.draw_text(txt, 20, WL.Screen.windowRect(4)-30, 'fontsize', 12, 'fliph', 0,'flipv', 0, 'col', [1 1 1]);
            %WL.draw_text(sprintf('%i',WL.Streak), WL.StreakPosnPx(1), WL.StreakPosnPx(2), 'col', WL.StreakColor);
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
                        %WL.ObjectPosition = WL.cfg.ObjectHomePosition + [0 0 -WL.Trial.ObjectRadius]';
                        WL.ObjectPosition = WL.cfg.ObjectHomePosition + WL.Trial.ObjectOffset;% + [0 0 -WL.Trial.ObjectRadius]';
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
%                     elseif WL.State.Timer.GetTime > WL.cfg.ObjectPickupTimeOut
%                         WL.state_next(WL.State.TIMEOUT);
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
                            WL.update_streak(1);
                        elseif mt > WL.cfg.MovementTooSlowWarning
                            WL.play_sound(WL.cfg.slowwarnbeep);
                            WL.update_streak(0);
                        else
                            WL.update_streak(1);
                            WL.play_sound(WL.cfg.placebeep);
                        end
                        
                        WL.ObjectPosition = WL.Robot.Position + WL.Trial.ObjectOffset; % Update for the last time
                        if ~isempty(WL.Explosion)
                            WL.Explosion.ExplodePop(WL.ObjectPosition);
                        end
                        WL.Trial.MovementDurationTime = WL.Timer.MovementDurationTimer.GetTime();
                        if( WL.Trial.FieldType ~= 3 ) % Not a PMove trial
                            WL.Trial.StimulusFrameCount = WL.FrameCounter.Stimulus.GetCount();
                            WL.Trial.StimulusTime = WL.Timer.Stimulus.Toc();
                            %WL.printf('Stimulus on for %d frames, %.3f ms (%.1f Hz).\n',WL.Trial.StimulusFrameCount,WL.Trial.StimulusTime,WL.Trial.StimulusFrameCount/WL.Trial.StimulusTime);
                        end
                        WL.state_next(WL.State.EXPLODE);
                    elseif WL.Timer.MovementDurationTimer.GetTime > 2.5 || (WL.Trial.Timed && WL.Timer.MovementDurationTimer.GetTime > WL.cfg.MovementTooSlowTimeout)
                        WL.MovementDurationTimeoutMessage = 'Too Slow';
                        WL.update_streak(0);
                        WL.state_next(WL.State.TIMEOUT);
                    elseif ~WL.cfg.StopAtTarget && WL.cfg.RequireSlice && norm(WL.Robot.Position-WL.cfg.ObjectHomePosition) > WL.cfg.TargetDistance
                        WL.ObjectPosition = WL.Robot.Position + WL.Trial.ObjectOffset; % Update for the last time
                        ok = WL.Robot.RampDown();
                        WL.MissWarn = true;
                        WL.play_sound(WL.cfg.tooslowbeep);
                        WL.Trial.MovementDurationTime = WL.Timer.MovementDurationTimer.GetTime();
                        if( WL.Trial.FieldType ~= 3 ) % Not a PMove trial
                            WL.Trial.StimulusFrameCount = WL.FrameCounter.Stimulus.GetCount();
                            WL.Trial.StimulusTime = WL.Timer.Stimulus.Toc();
                        end
                        WL.state_next(WL.State.FINISH);
                    end
                    
                case WL.State.EXPLODE
                    if WL.Explosion.ExplodeState==4
                        WL.state_next(WL.State.FINISH);
                    end
                    
                case WL.State.FINISH
                    %[ok, WL.GW.FrameData, names] = WL.Hardware.DataGet(); %doesn't work (or work fast enough?)
                    if WL.State.Timer.GetTime > WL.cfg.FinishDelay % Trial has finished so stop trial.
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

                            WL.MPEBuffer(WL.Trial.TrialNumber) = WL.compute_mpe(WL.Trial.TargetAngle);
                            WL.Score = round(max(0,100-100*WL.MPEBuffer(WL.Trial.TrialNumber)));
                            if WL.MissWarn
                                WL.Score = 0;
                            end
                            WL.ScoreBuffer(WL.Trial.TrialNumber) = WL.Score;
                            WL.Timer.ScoreFade.Reset();

                            WL.state_next(WL.State.NEXT);
%                             end
                        end
                    end
                    
                case WL.State.NEXT
                    if WL.MissWarn
                        if WL.State.Timer.GetTime < 2
                            WL.cfg.ClearColor(1) = abs(cos(WL.State.Timer.GetTime*4*pi)-1)/6;
                            return;
                        else
                            WL.cfg.ClearColor(1) = 0;
                        end
                    end
                    if WL.Trial.RestFlag==1
                        WL.state_next(WL.State.REST);
                    elseif  ~WL.trial_next()
                        WL.state_next(WL.State.EXIT);
                    else
                        WL.state_next(WL.State.INTERTRIAL);
                    end
                    
                case WL.State.INTERTRIAL % Wait for the intertrial delay to expire.
                    if WL.Timer.Paradigm.InterTrialDelayTimer.GetTime > WL.cfg.InterTrialDelay
                        if WL.SkipToTest
                            WL.TrialNumber = WL.cfg.FirstTestTrial;
                            WL.SkipToTest = false;
                        end
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
                    ok = WL.Robot.FieldChannel(WL.Trial.TargetPosition - WL.Trial.ObjectOffset,WL.Trial.FieldConstants(1),WL.Trial.FieldConstants(2));
                case 3 % Passive return movement
                    ok = WL.Robot.FieldPMove(WL.Trial.TargetPosition,0.5,0.2);
            end
            
            %WL.printf('TrialStart() Trial=%d, Field=%d\n',WL.TrialNumber,WL.Trial.FieldType);
            %WL.printf('RobotField=%d, Started=%d\n',WL.Trial.FieldType,ok);
            
            WL.Explosion = wl_explode(WL.Trial.ObjectRadius,WL.cfg.ObjectActiveColor,sqrt(WL.Trial.ObjectRadius^2/16));
            WL.MissWarn = false;
            
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
            if (StopAtTarget)
                err = norm(WL.ObjectPosition(1:2) - WL.Trial.TargetPosition(1:2));
                out = (err < WL.cfg.CursorRadius+WL.cfg.TargetRadius) &&  WL.robot_stationary();
            else
                if (norm(WL.Robot.Position-WL.cfg.ObjectHomePosition) > WL.cfg.TargetDistance)
                    d = WL.Robot.Position(1:2) - WL.RobotPreviousPosition(1:2);
                    f = WL.RobotPreviousPosition(1:2) - WL.cfg.ObjectHomePosition(1:2);
                    a = d' * d;
                    b = 2 * (f' * d);
                    c = f' * f - (WL.cfg.TargetDistance^2);
                    discriminant = sqrt(b*b-4*a*c);
                    t = (-b + discriminant) / (2 * a);
                    WL.Trial.Cross = WL.RobotPreviousPosition(1:2) + t*d;
                    err = norm(WL.Trial.Cross + WL.Trial.ObjectOffset(1:2) - WL.Trial.TargetPosition(1:2));
                    %err = norm(WL.ObjectPosition(1:2) - WL.Trial.TargetPosition(1:2));
                    out = (err < WL.cfg.CursorRadius+WL.cfg.TargetRadius);
                else
                    out = false;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = movement_started(WL)
            err = norm(WL.Robot.Position(1:2) - WL.cfg.ObjectHomePosition(1:2));
            out = err > 0.5;
            %out = ~WL.robot_at_object();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = update_streak(WL,in)
            if in==0
                WL.Streak = 0;
                WL.StreakColor = [1 1 1];
            else
                WL.Streak = WL.Streak + 1;
                if WL.Streak > WL.BestStreak
                    WL.BestStreak = WL.Streak;
                    %WL.StreakColor = [1 0 1];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function mpe = compute_mpe(WL, target_angle)
            rotationangle = -(target_angle - pi/2)*180/pi;
            if contains(version,'2017')
                rotmat = EulerRotationMatrix('z',rotationangle,'D',0);
            else
                rotmat = rotz(rotationangle);
            end
            xyz_r = rotmat * WL.GW.FrameData{1}.RobotPosition(:,WL.GW.FrameData{1}.State==WL.State.MOVING);
            mpe = max(abs(xyz_r(1,:)));
        end
    end
end





