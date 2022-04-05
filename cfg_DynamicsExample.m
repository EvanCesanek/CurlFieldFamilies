function cfg_DynamicsExample(WL,cfg_name)

%%%%%%%%%%%%%%% general parameters

WL.cfg.RobotName = 'ROBOT_vBOT-Left'; % VIOLET, BEGONIA rig
%WL.cfg.RobotName = 'ROBOT_Stiff2Bot-Left'; % HIBISCUS rig
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
%WL.cfg.trial_save = false;
%WL.cfg.verbose = 0;
%WL.cfg.vol = 0.5;
%WL.cfg.Debug = true;

%to rerun previous table 
%WL.GW.reload_table=false
%WL.GW.trials_to_run=[3 9]

% possible values for the next two found by running ResolutionTest
% WL.cfg.ScreenFrameRate=75; 
% WL.cfg.ScfeenResolution=[1280 1024]; 
% WL.cfg.ScreenIndex =2

%%%%%%%%%%%%%%% experiment speciific parameters
WL.cfg.CursorRadius = 0.25;
WL.cfg.HomeRadius = 0.5;
WL.cfg.TargetRadius = 0.33;

WL.cfg.HomeTolerance = WL.cfg.HomeRadius;
WL.cfg.TargetTolerance	= WL.cfg.TargetRadius;
WL.cfg.ObjectTolerance	= 0.375; % this is the minimum object radius
WL.cfg.StationarySpeed = 5; % cm/s
WL.cfg.StationaryTime = 0.1; % s

% Require movement to end stationary within target?
WL.cfg.StopAtTarget = false;

WL.cfg.MovementReactionTimeOut = 1.0;
WL.cfg.MovementTooSlowTimeout = 1.0;
WL.cfg.MovementTooSlowWarning = 0.7;
WL.cfg.MovementTooFastWarning = 0.3;
WL.cfg.MovementTooFastTimeout = 0.15;
WL.cfg.InterTrialDelay = 0.0;
WL.cfg.RestBreakSeconds = 45;
WL.cfg.TrialDelay = 0.15;
WL.cfg.FinishDelay = 0.33;
WL.cfg.ErrorWait = 1.5;

WL.cfg.ObjectPulseDuration = 0.15;

WL.cfg.TargetDistance = 10;

WL.cfg.NumTargets = 5;
targetspacing = 2*pi/WL.cfg.NumTargets;

WL.cfg.NumObjects = 5;
% Object radius should be spaced via sqrt, not linear, because perceived size ~ circle area 
%WL.cfg.ObjectRadius = linspace(1.5*WL.cfg.CursorRadius, 5*WL.cfg.CursorRadius, WL.cfg.NumObjects);
WL.cfg.ObjectRadius = WL.cfg.CursorRadius*sqrt(linspace(2.25,15,WL.cfg.NumObjects));
% Linear family of viscous gains across the objects (NB: centered on baseline gain, see below)
WL.cfg.ViscousGainModifier = 0.019 * ((1:WL.cfg.NumObjects) - median(1:WL.cfg.NumObjects));
% Outlier object
WL.cfg.OutlierId = 3;

WL.cfg.ObjectHomePosition = [0 0 0]';
WL.cfg.HomePosition = WL.cfg.ObjectHomePosition + [cos(-pi/4) sin(-pi/4) 0]'*WL.cfg.TargetDistance/3;

% Set the modifier for the outlier, e.g. equal to the greatest gain in the family
WL.cfg.ViscousGainModifier(WL.cfg.OutlierId) = WL.cfg.ViscousGainModifier(WL.cfg.NumObjects);
%WL.cfg.ViscousGainModifier(WL.cfg.OutlierId) = 2*WL.cfg.ViscousGainModifier(WL.cfg.NumObjects);

% Default object color/shape -- override below if desired
WL.cfg.ObjectInactiveColor = [1 0 1 0.45];
WL.cfg.ObjectActiveColor = [1 0 1 0.75];
WL.cfg.ObjectShape = 'circle';

try
    WL.cfg.highbeep = WL.load_beeps(500,0.05);
    WL.cfg.pickbeep = WL.load_beeps([330 660],[0.03 0.02]);
    WL.cfg.placebeep =  WL.load_beeps([660 1 880 1 880],[0.05 0.02 0.05 0.05 0.07]);
    WL.cfg.slowwarnbeep = WL.load_beeps([250],[0.25]);
    WL.cfg.fastwarnbeep = WL.load_beeps([4000 1 4000],[0.05 0.05 0.05]);
    WL.cfg.tooslowbeep = WL.load_beeps([250 150],[0.5 0.5]);
    WL.cfg.toofastbeep = WL.load_beeps([4000],[0.2]);
    WL.overide_cfg_defaults();
catch
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%items that can change with each trial - either specify indexed structure or indexed cell
%array - need to be made global

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch upper(cfg_name)
    
    case 'OUT1ISO'
        lowerlim = pi/5;
        upperlim = 3*pi/2;
        WL.cfg.TargetAngleRandomShift = lowerlim + rand*(upperlim - lowerlim - targetspacing);
    
    case {'OUT3ISO', 'OUT2ISO','OUT2ISOL'} % fixed 100 - 160 - 220 // 135 - 180 - 225
        targetspacing = pi/4; %pi/3;
        WL.cfg.TargetAngleRandomShift = 3*pi/4; %pi/2+10*pi/180;
        WL.cfg.NumTargets = 8;
    
    otherwise
        %WL.cfg.TargetAngleRandomShift = 0; % "target 1" is to the right
        WL.cfg.TargetAngleRandomShift = pi/2; % "target 1" is straight ahead
        %WL.cfg.TargetAngleRandomShift = rand*2*pi/WL.cfg.NumTargets; % "target 1" direction is randomized     
end

TargetAngle = WL.cfg.TargetAngleRandomShift:targetspacing:(4*pi);
TargetAngle = TargetAngle(1:WL.cfg.NumTargets);
TargetAngle = num2cell(TargetAngle);

ObjectId = num2cell(1:WL.cfg.NumObjects);

WL.cfg.Field = Field;
WL.cfg.TargetAngle = TargetAngle;
WL.cfg.ObjectId = ObjectId;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial Setup % parameters that are used in both visuomotor rotation and force field adaptation experiments (VMR/FF)

%PreExposure.Trial.Index.TargetAngle = [ 1 2 3 4 5 ]; % indexes TargetAngle = [0] [1.5708] [3.1416] [4.7124] [6.2832] which are radians
%PreExposure.Permute = true; % whether to permute within block
%PreExposure.Location = sparse(5,5);
%PreExposure.Adjacency = sparse(5,5);
%PreExposure.Location(5,1) = 1; % forbidden channel trial first or last % index 5 cannot be in position 1 and 2?
%PreExposure.Location(5,2) = 1;
%PreExposure.Adjacency(3,5) = 1; % forbidden movement to same trial % index 3 cannot be next to 5
%PreExposure.Adjacency(5,3) = 1; % forbidden channel trial first or last

%Exposure = PreExposure;
%PostExposure = PreExposure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch upper(cfg_name) % Specify parameters unique to each experiment (VMR/FF)
    
    case 'PRE' % Pre-training 1 target
        WL.cfg.ObjectInactiveColor = [0 1 0 0.45]; % Green object for familiarization
        WL.cfg.ObjectActiveColor = [0 1 0 0.75];
        WL.cfg.ObjectShape = 'circle';
        Pretraining.Trial.Index.ObjectId = 5;
        Pretraining.Trial.Index.TargetAngle = 1;
        Pretraining.Trial.Index.Field = 1; % Null field
        A = WL.parse_trials(Pretraining);
        T = parse_tree(30*A);
    
    case 'OUT1' % Outlier, single target
    case 'OUT1NT' % Outlier, single target, no training phase
        % Training phase
        % In each block, one field trial with each of four "family" objects
        %   plus one channel trial with one of the "family" objects
        % So to get one channel trial with each object it takes 4 blocks
        % Channel trial can't be first
        % Channel trial can't come before or after field trial with same object
        ji = 0;
        for j = [1 2 4 5]
            ji = ji+1;
            ExposureFam.Trial.Index.ObjectId = [1 2 4 5 j];
            ExposureFam.Trial.Index.TargetAngle = [1 1 1 1 1];
            ExposureFam.Trial.Index.Field = [2 2 2 2 3];
            ExposureFam.Permute = true;
            ExposureFam.Location = sparse(5,5);
            ExposureFam.Location(5,1) = 1;
            ExposureFam.Adjacency = sparse(5,5);
            ExposureFam.Adjacency(5,ji) = 1;
            ExposureFam.Adjacency(ji,5) = 1;
            A{ji} = WL.parse_trials(ExposureFam);
        end
        % Causes parse_tree to sample randomly without replacement
        WA = wl_node.permute([A{1} A{2} A{3} A{4}]);
        
        % Test phase
        % Same as training except we add in the outlier
        ji = 0;
        for j = [1 2 3 4 5]
            ji = ji+1;
            ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 j];
            ExposureAll.Trial.Index.TargetAngle = [1 1 1 1 1 1];
            ExposureAll.Trial.Index.Field = [2 2 2 2 2 3];
            ExposureAll.Permute = true;
            ExposureAll.Location = sparse(6,6);
            ExposureAll.Location(6,1) = 1;
            ExposureAll.Adjacency = sparse(6,6);
            ExposureAll.Adjacency(6,ji) = 1;
            ExposureAll.Adjacency(ji,6) = 1;
            B{ji} = WL.parse_trials(ExposureAll);
        end
        WB = wl_node.permute([B{1} B{2} B{3} B{4} B{5}]);
        
        % Alternating channel/field phase
        ChannelAll.Trial.Index.ObjectId    = [1 1 2 2 3 3 4 4 5 5];
        ChannelAll.Trial.Index.TargetAngle = [1 1 1 1 1 1 1 1 1 1];
        ChannelAll.Trial.Index.Field       = [2 3 2 3 2 3 2 3 2 3];
        ChannelAll.Permute = true;
        ChannelAll.Adjacency = sparse(10,10);
        ChannelAll.Location = sparse(10,10);
        for j = 2:2:10
            ChannelAll.Location(j,1) = 1; % channels can't go first (they will be the even trials)
            ChannelAll.Adjacency(j-1,j)  = 1; % don't repeat object in field-channel pair
            ChannelAll.Adjacency(j,j-1)  = 1; % don't repeat object in field-channel pair
            for jj = (j+2):2:10
                ChannelAll.Adjacency(j,jj)  = 1; % don't allow consecutive channels
                ChannelAll.Adjacency(jj,j)  = 1; % don't allow consecutive channels
            end
        end
        WC = WL.parse_trials(ChannelAll);
        
        if strcmpi(cfg_name,'OUT1')
            T = parse_tree(20*WA + 25*WB + 10*WC);
        elseif strcmpi(cfg_name,'OUT1NT')
            T = parse_tree(40*WB + 10*WC);
        end
        
        
    case 'OUT5'
        % Train phase: All five objects, separate targets
        WL.cfg.TOMap = 1+mod((1:5) + randi(5), 5); % random shift (so weakest is not always straight ahead)
        ji = 0;
        for j = [1 2 3 4 5]
            ji = ji+1;
            ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 j];
            ExposureAll.Trial.Index.TargetAngle = [WL.cfg.TOMap WL.cfg.TOMap(j)];
            ExposureAll.Trial.Index.Field = [2 2 2 2 2 3];
            ExposureAll.Permute = true;
            ExposureAll.Location = sparse(6,6);
            ExposureAll.Location(6,1) = 1;
            ExposureAll.Adjacency = sparse(6,6);
            ExposureAll.Adjacency(6,ji) = 1;
            ExposureAll.Adjacency(ji,6) = 1;
            B{ji} = WL.parse_trials(ExposureAll);
        end
        WB = wl_node.permute([B{1} B{2} B{3} B{4} B{5}]);
        
        % Generalize phase: Add channel trials to other locations
        % 5 sub-blocks of 5 field trials + 5 channels
        % Ensures diff types of field trials are more evenly spaced
        % Also avoids possibility of diff objects to same test location without recent field trial at that location
        clear tmp2
        for ti = 1:5
            tmp(ti,:) = Shuffle(WL.cfg.TOMap);
            % Ensure that all objects are tested at diff locations in each sub-block
            %  -- sorts columns and checks for any repeats (i.e. diff==0)
            while any(max(diff(sort(tmp,1),1,1)==0)) 
                tmp(ti,:) = Shuffle(tmp(ti,:));
            end
        end
        % Now each column of tmp lists the channel target for each object in each sub-block
        for bi = 1:5
            GeneralizeTarget.Trial.Index.ObjectId = [1:5 1:5];
            GeneralizeTarget.Trial.Index.TargetAngle = [WL.cfg.TOMap tmp(:,bi)'];
            GeneralizeTarget.Trial.Index.Field = [2 2 2 2 2 3 3 3 3 3];
            GeneralizeTarget.Permute = true;
            GeneralizeTarget.Location = sparse(10,10);
            GeneralizeTarget.Location(6:end,[1 4 7]) = 1;
%             GeneralizeTarget.Adjacency = sparse(10,10);
%             GeneralizeTarget.Adjacency(6:end,6:end) = 1;
            D{bi} = WL.parse_trials(GeneralizeTarget);
        end
        WD = wl_node.permute([D{1} D{2} D{3} D{4} D{5}]);
        
        %T = parse_tree(40*WB + 20*WD);
        
        T = parse_tree(5*WB + 5*WD);
        %T = parse_tree(1000*WD);
        
    case 'OUT1ISO' % Outlier, single target + isolated outlier target, no training phase
        WL.cfg.ISOtarget = 5;
        % Early phase - all fields
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget WL.cfg.ISOtarget WL.cfg.ISOtarget];
        ExposureAll.Trial.Index.TargetAngle = [1 1 1 1 1 2 2 2];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2 2 2 2];
        ExposureAll.Permute = true;
        WA = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget WL.cfg.ISOtarget WL.cfg.ISOtarget 1 2 3 4 5 WL.cfg.ISOtarget];
        ExposureTest.Trial.Index.TargetAngle = [1 1 1 1 1 2 2 2 1 1 1 1 1 2];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 3 3 3 3 3 3];
        ExposureTest.Permute = true;
        ExposureTest.Location = sparse(14,14);
        ExposureTest.Location(9:end,[1 4 7 10 13]) = 1;
        WB = WL.parse_trials(ExposureTest);
        
        GeneralizeTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget WL.cfg.ISOtarget WL.cfg.ISOtarget 1 2 3 4 5 WL.cfg.ISOtarget];
        GeneralizeTest.Trial.Index.TargetAngle = [1 1 1 1 1 2 2 2 2 2 2 2 2 1];
        GeneralizeTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 3 3 3 3 3 3];
        GeneralizeTest.Permute = true;
        GeneralizeTest.Location = sparse(14,14);
        GeneralizeTest.Location(9:end,[1 4 7 10 13]) = 1;
        WC = WL.parse_trials(GeneralizeTest);
        
        %T = parse_tree(2*WA + 1*WB + 1*WC);
        T = parse_tree(25*WA + 10*WB + 2*WA + 10*WC);
        
    case 'OUT2ISO' % Outlier, single target + isolated outlier target, no training phase
        WL.cfg.ISOtarget = [3 5];
        
        ExposureFam.Trial.Index.ObjectId = [1 2 4 5 WL.cfg.ISOtarget([2 2])];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 3 3];
        ExposureFam.Trial.Index.Field = [2 2 2 2 2 2];
        ExposureFam.Permute = true;
        W0 = WL.parse_trials(ExposureFam);
        
        % Early phase - all fields
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
        
        %T = parse_tree(2*WA + 1*WB + 1*WC);
        %T = parse_tree(30*W0 + 50*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        T = parse_tree(15*W0 + 25*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        
    case 'OUT2ISOL' % Outlier, single target + isolated outlier target, no training phase
        WL.cfg.ISOtarget = [3 1];
        
        ExposureFam.Trial.Index.ObjectId = [1 2 4 5 WL.cfg.ISOtarget([2 2])];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 7 7];
        ExposureFam.Trial.Index.Field = [2 2 2 2 2 2];
        ExposureFam.Permute = true;
        W0 = WL.parse_trials(ExposureFam);
        
        % Early phase - all fields
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2])];
        ExposureAll.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 7 7 7];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2];
        ExposureAll.Permute = true;
        WA = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget(1) WL.cfg.ISOtarget(2)];
        ExposureTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 7 7 7 2 2 2 2 2 1 7];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3];
        ExposureTest.Permute = true;
        ExposureTest.Location = sparse(length(ExposureTest.Trial.Index.Field),length(ExposureTest.Trial.Index.Field));
        ExposureTest.Location(find(ExposureTest.Trial.Index.Field==3,1):end, 1:3:end) = 1;
        WB = WL.parse_trials(ExposureTest);
        
        GeneralizeTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1]) WL.cfg.ISOtarget([2 2 2]) 1 2 3 4 5 1 2 3 4 5];
        GeneralizeTest.Trial.Index.TargetAngle = [2 2 2 2 2 1 1 1 7 7 7 1 1 1 1 1 7 7 7 7 7];
        GeneralizeTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
        GeneralizeTest.Permute = true;
        GeneralizeTest.Location = sparse(length(GeneralizeTest.Trial.Index.Field),length(GeneralizeTest.Trial.Index.Field));
        GeneralizeTest.Location(find(GeneralizeTest.Trial.Index.Field==3,1):end, 1:3:end) = 1;
        WC = WL.parse_trials(GeneralizeTest);
        
        %T = parse_tree(2*WA + 1*WB + 1*WC);
        %T = parse_tree(30*W0 + 50*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        T = parse_tree(25*WA + 5*WB + 5*WC + 5*WB + 5*WC);

    case 'OUT3ISO'
        WL.cfg.ISOtarget = [3 5]; %[3 5] or [3 1]
        
        if all(WL.cfg.ISOtarget==[3 5])
            WL.cfg.ISOposn = [1 3];
        elseif all(WL.cfg.ISOtarget==[3 1])
            WL.cfg.ISOposn = [1 7];
        end
        
        ExposureFam.Trial.Index.ObjectId = [1 2 4 5 WL.cfg.ISOtarget([2 2])];
        ExposureFam.Trial.Index.TargetAngle = [2 2 2 2 WL.cfg.ISOposn([2 2])];
        ExposureFam.Trial.Index.Field = [2 2 2 2 2 2];
        ExposureFam.Permute = true;
        W0 = WL.parse_trials(ExposureFam);
        
        ExposureAll.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1 2 2 2])];
        ExposureAll.Trial.Index.TargetAngle = [2 2 2 2 2 WL.cfg.ISOposn([1 1 1 2 2 2])];
        ExposureAll.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2];
        ExposureAll.Permute = true;
        WA = WL.parse_trials(ExposureAll);
        
        ExposureTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1 2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget([1 1 1 2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget([1 2])];
        ExposureTest.Trial.Index.TargetAngle = [2 2 2 2 2 WL.cfg.ISOposn([1 1 1 2 2 2]) 2 2 2 2 2 WL.cfg.ISOposn([1 1 1 2 2 2]) 2 2 2 2 2 WL.cfg.ISOposn([1 2])];
        ExposureTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3];
        ExposureTest.Permute = true;
        ExposureTest.Location = sparse(length(ExposureTest.Trial.Index.Field),length(ExposureTest.Trial.Index.Field));
        ExposureTest.Location(find(ExposureTest.Trial.Index.Field==3,1):end, 1:2:end) = 1;
        WB = WL.parse_trials(ExposureTest);
        
        GeneralizeTest.Trial.Index.ObjectId = [1 2 3 4 5 WL.cfg.ISOtarget([1 1 1 2 2 2]) 1 2 3 4 5 WL.cfg.ISOtarget([1 1 1 2 2 2]) 1 2 3 4 5 1 2 3 4 5];
        GeneralizeTest.Trial.Index.TargetAngle = [2 2 2 2 2 WL.cfg.ISOposn([1 1 1 2 2 2]) 2 2 2 2 2 WL.cfg.ISOposn([1 1 1 2 2 2]) WL.cfg.ISOposn([1 1 1 1 1 2 2 2 2 2])];
        GeneralizeTest.Trial.Index.Field = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3];
        GeneralizeTest.Permute = true;
        GeneralizeTest.Location = sparse(length(GeneralizeTest.Trial.Index.Field),length(GeneralizeTest.Trial.Index.Field));
        GeneralizeTest.Location(find(GeneralizeTest.Trial.Index.Field==3,1):end, 1:2:end) = 1;
        WC = WL.parse_trials(GeneralizeTest);
        
        %T = parse_tree(20*WA + 5*WB + 5*WC + 5*WB + 5*WC);
        T = parse_tree(30*W0 + 50*WA + 5*WB + 5*WC + 5*WB + 5*WC);

    otherwise
        error('cfg name invalid')
end

z = zeros(rows(T),1);

%create more table parameters
R = [cos(T.TargetAngle) sin(T.TargetAngle) z];
T.TargetPosition = bsxfun(@plus, bsxfun(@times, WL.cfg.TargetDistance,R), WL.cfg.ObjectHomePosition');
WL.cfg.TargetPositions = unique(T.TargetPosition, 'rows');

%T.ObjectReactionTime = z;
T.MovementReactionTime = z;
T.MovementDurationTime = z;

% Update field constants according to family modifiers
T.FieldConstants(T.FieldType==1, 1) = T.FieldConstants(T.FieldType==1, 1) + WL.cfg.ViscousGainModifier(T.ObjectId(T.FieldType==1))';
T.ObjectRadius = WL.cfg.ObjectRadius(T.ObjectId)';

T.TrialNumber = (1:size(T,1))';
%T.ObjectRotate = randi(360,size(T,1),1);

%add passive return movements
%T = wl_movement_return(WL.cfg, T, 'passive');

%T.RestFlag([10 40])=true


WL.TrialData = T;