function cfg_ControlPoints(WL,cfg_name)

%%%%%%%%%%%%%%% general parameters

%WL.cfg.RobotName = 'ROBOT_vBOT-Left'; % VIOLET, BEGONIA rig
WL.cfg.RobotName = 'ROBOT_Stiff2Bot-Left'; % HIBISCUS rig
%WL.cfg.RobotName = 'ROBOT_Stiff2Bot-Right'; % SUNFLOWER rig
%WL.cfg.RobotName = 'ROBOT_WristBOT-Left'; % BUTTERCUP rig
%WL.cfg.RobotName = 'ROBOT_3BOT-01'; % PANSY rig
%WL.cfg.RobotName = 'ROBOT_3BOT-11'; % LAVENDER rig
%WL.cfg.RobotName = 'ROBOT_3BOT-10'; % GERANIUM rig
%WL.cfg.RobotName = 'ROBOT_Small-Left'; % CROCUS rig
%WL.cfg.RobotForceMax = 10; % IMPORTANT! When testing, limit force to 10N

%WL.cfg.OculusRift = true;
WL.cfg.RotateScreen = true; % true on VIOLET & BEGONIA rigs
WL.cfg.SkipSyncTests = true;
WL.cfg.ScreenIndex = 0;
%WL.cfg.MouseFlag = true;
%WL.cfg.SmallScreen = false;
%WL.cfg.ScreenSize = 2*[100 100 640 640*16.8/29.8]; 
%WL.cfg.ClearColor = [ 0 0 0 ];
%WL.cfg.trial_save = true;
%WL.cfg.verbose = 0;
%WL.cfg.vol = 0.5;
%WL.cfg.Debug = true;

%to rerun previous table 
%WL.GW.reload_table=false
%WL.GW.trials_to_run=[3 9]

% possible values for the next two found by running ResolutionTest
% WL.cfg.ScreenFrameRate=75; 
% WL.cfg.ScreenResolution=[1280 1024]; 
% WL.cfg.ScreenIndex =2

%%%%%%%%%%%%%%% experiment speciific parameters
WL.cfg.CursorRadius = 0.25;
WL.cfg.HomeRadius = 0.33;
WL.cfg.TargetRadius = 0.33;

WL.cfg.RectWidth = 8;
WL.cfg.RectHeight = 1;

WL.cfg.HomeTolerance = WL.cfg.HomeRadius;
WL.cfg.TargetTolerance	= WL.cfg.TargetRadius;
WL.cfg.StationarySpeed = 5; % cm/s
WL.cfg.StationaryTime = 0.1; % s

WL.cfg.MovementReactionTimeOut = 1.0;
WL.cfg.MovementDurationTimeOut = 5.0;
WL.cfg.InterTrialDelay = 0.25;
WL.cfg.RestBreakSeconds = 45;
WL.cfg.TrialDelay = 0.15;
WL.cfg.FinishDelay = 0.0;
WL.cfg.ErrorWait = 1.5;

WL.cfg.TargetDistance = 10;
WL.cfg.HomePosition = [0 0 0]';

WL.cfg.NumObjects = 5;
% Object radius should be spaced via sqrt, not linear, because perceived size ~ circle area 
%WL.cfg.ObjectRadius = linspace(1.5*WL.cfg.CursorRadius, 5*WL.cfg.CursorRadius, WL.cfg.NumObjects);
WL.cfg.ObjectRadius = 0.5;
% Linear family of viscous gains across the objects (NB: centered on baseline gain, see below)
WL.cfg.ViscousGainModifier = 0.019 * ((1:WL.cfg.NumObjects) - median(1:WL.cfg.NumObjects));
% Outlier object
WL.cfg.OutlierId = 3;
% Set the modifier for the outlier
WL.cfg.ViscousGainModifier(WL.cfg.OutlierId) = WL.cfg.ViscousGainModifier(WL.cfg.NumObjects); % HEAVY
%WL.cfg.ViscousGainModifier(WL.cfg.OutlierId) = 2*WL.cfg.ViscousGainModifier(WL.cfg.NumObjects); % EXTRA HEAVY 

% Default object color/shape -- override below if desired
WL.cfg.ObjectInactiveColor = [1 0 1 0.45];
WL.cfg.ObjectActiveColor = [1 0 1 0.75];

WL.cfg.HomeActiveColor = [0.8 0.8 0.8 0.7];
WL.cfg.HomeInactiveColor = [0.7 0.7 0.7 0.2];

try
    WL.cfg.highbeep = WL.load_beeps(500,0.05);
    WL.cfg.pickbeep = WL.load_beeps([330 660],[0.03 0.02]);
    WL.cfg.placebeep = WL.load_beeps([440 660],[0.05 0.12]);
    WL.cfg.lowbeep = WL.load_beeps([250 150],[0.5 0.5]);
    WL.overide_cfg_defaults();
catch
    
end

% Null field.
Field{1}.FieldType = 0;
Field{1}.FieldConstants	= [0 0];
% Field{1}.FieldAngle	= 0.0;
% Field{1}.Rotation	= 0;
% Field{1}.FieldMartrix =	eye(3);

% Viscous curl field
Field{2}.FieldType = 1;
Field{2}.FieldConstants	= [0.15 0]; % NB: this is the *baseline* viscous gain (see note above)
% Field{2}.FieldAngle	= 90.0;
% Field{2}.Rotation	= 0;
% Field{2}.FieldMartrix =	eye(3);

% Channel trial
Field{3}.FieldType = 2;
Field{3}.FieldConstants	= [-30.000  -0.05];
% Field{3}.FieldAngle	= 0.0;
% Field{3}.Rotation	= 0;
% Field{3}.FieldMartrix =	eye(3);

% FieldType = 3 are passive return trials this is defined in wl_movement_return
% Perhaps there should be a comment somewhere explicitly stating that
% FirldType 3 is reserved for this in case students attempt to define it as
% else?

WL.cfg.Field = Field;

WL.cfg.ObjectId = num2cell(1:WL.cfg.NumObjects);

WL.cfg.TargetAngle = num2cell([pi/4 pi/2 3*pi/4]);

switch upper(cfg_name) % Specify parameters unique to each experiment (VMR/FF)
    
    case 'OUT1' % Outlier, single target
        ExposureFam.Trial.Index.ObjectId = [1 2 4 5];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2];
        ExposureFam.Trial.Index.Field = [2 2 2 2];
        ExposureFam.Permute = true;
        W0 = WL.parse_trials(ExposureFam);
        
        ExposureFamTest.Trial.Index.ObjectId = [1 2 4 5 1 2 4 5];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 2 2 2 2];
        ExposureFamTest.Trial.Index.Field = [2 2 2 2 3 3 3 3];
        ExposureFamTest.Permute = true;
        bl = length(ExposureFamTest.Trial.Index.Field); % block length
        fc = find(ExposureFamTest.Trial.Index.Field==3,1); % first channel trial index
        ExposureFamTest.Location = sparse(bl,bl);
        ExposureFamTest.Location(fc:end, 1:3:end) = 1; % no more than two channels in a row
        ExposureFamTest.Adjacency = sparse(bl,bl);
        for ti = 1:(fc-1)
            ExposureFamTest.Adjacency(ti,ti+fc) = 1; % channels should not immediately follow fields of same object
        end
        W1 = WL.parse_trials(ExposureFamTest);
        
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5];
        ExposureAll.Trial.Index.TargetAngle = [2 2 2 2 2];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2];
        ExposureAll.Permute = true;
        W2 = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 1 2 3 4 5];
        ExposureTest.Trial.Index.TargetAngle = [2 2 2 2 2 2 2 2 2 2];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 3 3 3 3 3];
        ExposureTest.Permute = true;
        bl = length(ExposureTest.Trial.Index.Field); % block length
        fc = find(ExposureTest.Trial.Index.Field==3,1); % first channel trial index
        ExposureTest.Location = sparse(bl,bl);
        ExposureTest.Location(fc:end, 1:3:end) = 1; % no more than two channels in a row
        ExposureTest.Adjacency = sparse(bl,bl);
        for ti = 1:(fc-1)
            ExposureTest.Adjacency(ti,ti+fc) = 1; % channels should not immediately follow fields of same object
        end
        W3 = WL.parse_trials(ExposureTest);
        %T = parse_tree(50*W0 + 10*W1 + 50*W0 + 10*W1 + 50*W0 + 10*W1 + 50*W0 + 10*W1); % 1120 = 960 field (240/object) + 160 channel (40/object)
        T = parse_tree(25*W0 + 50*W2 + 10*W3 + 50*W2 + 10*W3 + 50*W2 + 10*W3 + 50*W2 + 10*W3); % 1500 = 100 fam field (25/object) + 1200 field (240/object) + 200 channel (40/object)
        

    case 'OUT3ISO' % Outlier, single target for all objects + two single-object targets
        WL.cfg.ISOtarget = [3 5];
        
        ExposureFam.Trial.Index.ObjectId = [1 2 4 5 WL.cfg.ISOtarget([2 2])];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 3 3];
        ExposureFam.Trial.Index.Field = [2 2 2 2 2 2];
        ExposureFam.Permute = true;
        W0 = WL.parse_trials(ExposureFam);
        
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2])];
        ExposureAll.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2];
        ExposureAll.Permute = true;
        WA = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget(1) WL.cfg.ISOtarget(2)];
        ExposureTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3 2 2 2 2 2 1 3];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3];
        ExposureTest.Permute = true;
        ExposureTest.Location = sparse(length(ExposureTest.Trial.Index.Field),length(ExposureTest.Trial.Index.Field));
        ExposureTest.Location(find(ExposureTest.Trial.Index.Field==3,1):end, 1:3:end) = 1;
        WB = WL.parse_trials(ExposureTest);
        
        GeneralizeTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 1 2 3 4 5];
        GeneralizeTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3 1 1 1 1 1 3 3 3 3 3];
        GeneralizeTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
        GeneralizeTest.Permute = true;
        GeneralizeTest.Location = sparse(length(GeneralizeTest.Trial.Index.Field),length(GeneralizeTest.Trial.Index.Field));
        GeneralizeTest.Location(find(GeneralizeTest.Trial.Index.Field==3,1):end, 1:3:end) = 1;
        WC = WL.parse_trials(GeneralizeTest);
        
        %T = parse_tree(30*W0 + 50*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        T = parse_tree(8*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        
    case 'OUT3ISO_TEST' % Outlier, single target for all objects + two single-object targets
        WL.cfg.ISOtarget = [3 5];
        
%         ExposureFam.Trial.Index.ObjectId = [1 2 4 5 WL.cfg.ISOtarget([2 2])];
%         ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 3 3];
%         ExposureFam.Trial.Index.Field = [2 2 2 2 2 2];
%         ExposureFam.Permute = true;
%         W0 = WL.parse_trials(ExposureFam);
        
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2])];
        ExposureAll.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2];
        ExposureAll.Permute = true;
        WA = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget(1) WL.cfg.ISOtarget(2)];
        ExposureTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3 2 2 2 2 2 1 1 1 3 3 3 2 2 2 2 2 1 3];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3];
        ExposureTest.Permute = true;
        ExposureTest.Location = sparse(length(ExposureTest.Trial.Index.Field),length(ExposureTest.Trial.Index.Field));
        ExposureTest.Location(find(ExposureTest.Trial.Index.Field==3,1):end, 1:2:end) = 1;
        WB = WL.parse_trials(ExposureTest);
        
        GeneralizeTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 1 2 3 4 5];
        GeneralizeTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 3 3 3 2 2 2 2 2 1 1 1 3 3 3 1 1 1 1 1 3 3 3 3 3];
        GeneralizeTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
        GeneralizeTest.Permute = true;
        GeneralizeTest.Location = sparse(length(GeneralizeTest.Trial.Index.Field),length(GeneralizeTest.Trial.Index.Field));
        GeneralizeTest.Location(find(GeneralizeTest.Trial.Index.Field==3,1):end, 1:2:end) = 1;
        WC = WL.parse_trials(GeneralizeTest);
        
        %T = parse_tree(30*W0 + 50*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        T = parse_tree(8*WA + 5*WB + 5*WC + 5*WB + 5*WC);

    otherwise
        error('cfg name invalid')
end

z = zeros(rows(T),1);

%create more table parameters
xoffsets = linspace(-1, -WL.cfg.RectWidth+1, WL.cfg.NumObjects);
T.ObjectOffset = xoffsets(T.ObjectId)' .* [cos(T.TargetAngle-pi/2) sin(T.TargetAngle-pi/2) zeros(size(T,1),1)];
T.TargetPosition = WL.cfg.HomePosition' + T.ObjectOffset + WL.cfg.TargetDistance * [cos(T.TargetAngle) sin(T.TargetAngle) zeros(size(T,1),1)];
T.TargetAngleDegreesCentered = 180/pi * (T.TargetAngle-pi/2);

T.MovementReactionTime = z;
T.MovementDurationTime = z;

% Update field constants according to family modifiers
T.FieldConstants(T.FieldType==1, 1) = T.FieldConstants(T.FieldType==1, 1) + WL.cfg.ViscousGainModifier(T.ObjectId(T.FieldType==1))';
%T.ObjectRadius = WL.cfg.ObjectRadius(T.ObjectId)';

T.TrialNumber = (1:size(T,1))';

%add passive return movements
%T = wl_movement_return(WL.cfg, T, 'passive');

%T.RestFlag([10 40])=true


WL.TrialData = T;
