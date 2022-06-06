function TrialData = analyze_subject(subi,sesi,fn,varargin)
% ANALYZE_SUBJECT  Perform WL subject data pre-processing for vBOT curl field adaptation experiments
%   TrialData = ANALYZE_SUBJECT(subi, sesi, fn, [fn2][, datadir])
%
%   Returns an expanded TrialData, including a variety of additional pre-processed features/data.
%   (e.g., aligned and interpolated trajectories for averaging, MPE, CumForce, Adaptation Index, etc.) 
%
%   Arguments:
%   'subi' - integer label for the Subject field of TrialData
%   'sesi' - integer label for the Session field of TrialData
%   'fn'   - matfile name, output of WL experiment
%
%   Optional:
%   'fn2'     - a second matfile, if WL experiment had to be restarted
%   'datadir' - directory where fn and fn2 are located (default: './Data')

    %% Set parameters
    fn2 = [];
    datadir =  './Data/';

    if nargin>=4 && ~isempty(varargin{1})
        fn2 = varargin{1};
    end

    if nargin>=5 && ~isempty(varargin{2})
        datadir = varargin{2};
    end

    %% Load input files
    load([datadir fn], 'TrialData', 'TimeStamp', 'RobotPosition', 'RobotVelocity', 'RobotForces');
    if fn2
        TrialData_pre = TrialData;
        TimeStamp_pre = TimeStamp;
        sizepre = size(TimeStamp_pre);
        RobotPosition_pre = RobotPosition;
        RobotVelocity_pre = RobotVelocity;
        RobotForces_pre = RobotForces;
    
        load([datadir fn2], 'TrialData', 'TimeStamp', 'RobotPosition', 'RobotVelocity', 'RobotForces');
        sizepost = size(TimeStamp);
        if sizepre(2) > sizepost(2)
            nanpad = nan(sizepost(1),3,sizepre(2)-sizepost(2));
            TimeStamp = [TimeStamp squeeze(nanpad(:,1,:))];
            RobotPosition = cat(3, RobotPosition, nanpad);
            RobotVelocity = cat(3, RobotVelocity, nanpad);
            RobotForces = cat(3, RobotForces, nanpad);
        elseif sizepre(2) < sizepost(2)
            nanpad = nan(sizepre(1),3,sizepost(2)-sizepre(2));
            TimeStamp_pre = [TimeStamp_pre squeeze(nanpad(:,1,:))];
            RobotPosition_pre = cat(3, RobotPosition_pre, nanpad);
            RobotVelocity_pre = cat(3, RobotVelocity_pre, nanpad);
            RobotForces_pre = cat(3, RobotForces_pre, nanpad);
        end
        TrialData.TrialNumber = TrialData.TrialNumber + max(TrialData_pre.TrialNumber);
        TrialData.block_count = TrialData.block_count + max(TrialData_pre.block_count);
        TrialData = [TrialData_pre; TrialData];
        TimeStamp = [TimeStamp_pre; TimeStamp];
        RobotPosition = [RobotPosition_pre; RobotPosition];
        RobotVelocity = [RobotVelocity_pre; RobotVelocity];
        RobotForces = [RobotForces_pre; RobotForces];
    end
    
    TrialData.MPE = nan(size(TrialData,1),1);
    TrialData.TPE = TrialData.MPE;
    TrialData.PEPV = TrialData.MPE;
    TrialData.CumForce = TrialData.MPE;
    TrialData.IdealCumForce = TrialData.MPE;
    TrialData.IdealInterpCumForce = TrialData.MPE;
    TrialData.ForceAtPeakVel = TrialData.MPE;
    TrialData.AdaptIndex = TrialData.MPE;

    TrialData.Session = repelem(sesi,size(TrialData,1),1);
    TrialData.Subject = repelem(subi,size(TrialData,1),1);
    TrialData.ChannelTrial = TrialData.FieldType==2;
    TrialData.FieldDirection = sign(TrialData.FieldConstants(:,1));
    
    %% Preprocess
    angles = unique(TrialData.TargetAngle)';
    objects = unique(TrialData.ObjectId)';
    numobjects = max(objects);
    interplength = 100; % How many interpolated data points?

    TrialData.X = nan(size(TrialData,1),interplength);
    TrialData.Y = nan(size(TrialData,1),interplength);
    TrialData.Force = nan(size(TrialData,1),interplength);
    TrialData.ForceIdeal = nan(size(TrialData,1),interplength);
    TrialData.ForceIdealInterp = nan(size(TrialData,1),interplength);
    
    objectviscousgain = nan(numobjects,1);
    tmp = unique(TrialData.FieldDirection(~TrialData.ChannelTrial & TrialData.ObjectId~=3));
    if length(tmp)~=1
        warning('ERROR: Family objects do not all share same field direction.');
    end
    outlierinterpviscousgain = 0.15*tmp;
    
    for ai = 1:length(angles)
        target = TrialData.TargetAngle==angles(ai); %abs(TrialData.TargetAngle-angles(ai))<0.0001;
    
        for oi = objects
            object = TrialData.ObjectId==oi;
            % viscous gain associated with the object
            objectviscousgain(oi) = unique(TrialData.FieldConstants(~TrialData.ChannelTrial & object,1));

            testtrialsel = ~TrialData.MissTrial & TrialData.ChannelTrial & object & target;
            expotrialsel = ~TrialData.MissTrial & ~TrialData.ChannelTrial & object & target;
            trialselnames = {'channel trials', 'field trials'};

            tsi = 0;
            for ts = {testtrialsel,expotrialsel}
                tsi = tsi+1;
                trialsel = ts{1};
            
                % keep track of the number of trials included in trialsel
                fprintf('Found %i %s with target angle = %.2f, object ID = %i.\n',sum(trialsel),trialselnames{tsi},angles(ai),oi);
        
                if ~any(trialsel)
                    fprintf('  No trials found, skipping this angle-object combination.\n');
                    continue % skip this object-angle combo if we didn't use it
                end
        
                % **** Preprocessing ****
                
                % positions
                x = squeeze(RobotPosition(trialsel, 1, :));
                y = squeeze(RobotPosition(trialsel, 2, :));
                % velocities
                xv = squeeze(RobotVelocity(trialsel, 1, :));
                yv = squeeze(RobotVelocity(trialsel, 2, :));
                % forces
                xf = squeeze(RobotForces(trialsel, 1, :));
                yf = squeeze(RobotForces(trialsel, 2, :));
                if size(x,2)==1
                    x = x';
                    y = y';
                    xv = xv';
                    yv = yv';
                    xf = xf';
                    yf = yf';
                end
                % target angles
                thetas = TrialData.TargetAngle(trialsel);
                % times
                t = TimeStamp(trialsel,:);
            
                % rotate
                for ti = 1:length(thetas)
                    z = zeros(1,size(x,2));
                    rotationangle = -(thetas(ti) - pi/2)*180/pi;
                    if exist('rotz','builtin')
                        rotmat = rotz(rotationangle);
                    else
                        rotmat = EulerRotationMatrix('z',rotationangle);
                    end
                    xyz_r = rotmat * [x(ti,:); y(ti,:); z];
                    xyzv_r = rotmat * [xv(ti,:); yv(ti,:); z];
                    xyzf_r = rotmat * [xf(ti,:); yf(ti,:); z];
                    x(ti,:) = xyz_r(1,:);
                    y(ti,:) = xyz_r(2,:);
                    xv(ti,:) = xyzv_r(1,:);
                    yv(ti,:) = xyzv_r(2,:);
                    xf(ti,:) = xyzf_r(1,:);
                    yf(ti,:) = xyzf_r(2,:);
                end
            
                % Manual cropping using custom criteria (requires trajectory rotation to pi/2)  
                % (note that cropping based on WL State == MOVING doesn't work well) 
                ystartmove = 0.5;
                yendmove = 9;
                ystart = y > ystartmove;
                ystartcell = mat2cell(ystart,ones(size(ystart,1),1));
                startframe = cellfun(@(x) find(x,1), ystartcell);
                ynear = y > yendmove;
                ynearcell = mat2cell(ynear,ones(size(ynear,1),1));
                endframe = cellfun(@(x) find(x,1), ynearcell);
                for j = 1:size(x,1) %numtrialsincluded?
                    % tail
                    t(j,endframe(j):end) = nan;
                    x(j,endframe(j):end) = nan;
                    y(j,endframe(j):end) = nan;
                    xf(j,endframe(j):end) = nan;
                    yv(j,endframe(j):end) = nan;
                    % head
                    t(j,1:startframe(j)) = nan;
                    x(j,1:startframe(j)) = nan;
                    y(j,1:startframe(j)) = nan;
                    xf(j,1:startframe(j)) = nan;
                    yv(j,1:startframe(j)) = nan;
                end
                % Compute ideal force profile, based on gain and y-velocity profiles
                xfideal = -objectviscousgain(oi)*yv;
                % set t=0 at first "ismoving" frame
                tcell = mat2cell(t,ones(size(t,1),1));
                t = cell2mat(cellfun(@(x) x-min(x), tcell, 'UniformOutput', false));
    
                % **** Analysis ****
    
                % Compute integral of force wrt time
                xfcum = sum( (xf+[diff(xf,[],2)/2 nan(size(xf,1),1)]) .* [diff(t,[],2) nan(size(xf,1),1)], 2, 'omitnan');
                xfidealcum = -unique(TrialData.FieldDirection(trialsel))*sum( (xfideal+[diff(xfideal,[],2)/2 nan(size(xf,1),1)]) .* [diff(t,[],2) nan(size(xf,1),1)], 2, 'omitnan');
                xfinterpideal = -objectviscousgain(oi)*yv; % Compute ideal force profile, based on ASSUMED INTERPOLATED GAIN and y-velocity profiles
                if oi==3
                    xfinterpideal = -outlierinterpviscousgain*yv; % interpolated viscous gain for the outlier
                end
                xfidealinterpcum = sum( (xfinterpideal+[diff(xfinterpideal,[],2)/2 nan(size(xf,1),1)]) .* [diff(t,[],2) nan(size(xf,1),1)], 2, 'omitnan');
                TrialData.CumForce(trialsel) = xfcum-xfidealcum+mean(xfidealcum); % Compute it as a within-trial difference, and add the mean back on (so it's positive linear, not error clustered around 0)
                TrialData.IdealCumForce(trialsel) = mean(xfidealcum);
                TrialData.IdealInterpCumForce(trialsel) = mean(xfidealinterpcum);
    
                % Force at peak velocity
                [~, peakvelindices] = max(yv,[],2);
                peakvelindices = sub2ind(size(yv),1:size(yv,1),peakvelindices');
                TrialData.ForceAtPeakVel(trialsel) = xf(peakvelindices);
    
                % Perpendicular errors (terminal, maximum, and at peak velocity) 
                TrialData.TPE(trialsel) = x(sub2ind(size(t),1:size(t,1),cellfun(@(x) find(~isnan(x),1,'last'), num2cell(t',1))));
                if unique(TrialData.FieldDirection(trialsel)==1)
                    TrialData.MPE(trialsel) = min(x,[],2);
                elseif unique(TrialData.FieldDirection(trialsel)==-1)
                    TrialData.MPE(trialsel) = max(x,[],2);
                end
                TrialData.PEPV(trialsel) = x(peakvelindices);
                
                % Adaptation index
                adaptindex = nan(size(x,1),1);
                for ti = 1:size(x,1)
                    mdl = fitlm(array2table([xfideal(ti,:)' xf(ti,:)']),'Var2 ~ Var1 - 1');
                    adaptindex(ti) = mdl.Coefficients.Estimate;
                end
                TrialData.AdaptIndex(trialsel) = adaptindex;
            
                % Interpolating trajectories to fixed length for averaging
                xout = linspace(ystartmove,yendmove,interplength);
                xfi = nan(size(x,1),length(xout));
                yvi = nan(size(x,1),length(xout));
                xi = nan(size(x,1),length(xout));
                yi = nan(size(x,1),length(xout));
                for ti = 1:size(x,1)
                    yk = y(ti,:); % predictor ("x") variable is the y position
                    [~,idx] = unique(yk); % predictor values must be unique
                    idx = intersect(idx,find(~isnan(yk))); % ... and they can't be nans
                    xfi(ti,:) = interp1(y(ti,idx)', xf(ti,idx)', xout);
                    yvi(ti,:) = interp1(y(ti,idx)', yv(ti,idx)', xout);
                    xi(ti,:)  = interp1(y(ti,idx)', x(ti,idx)',  xout);
                    yi(ti,:)  = interp1(y(ti,idx)', y(ti,idx)',  xout);
                end
    
                % compute ideal interpolated force profile
                xfideali = -objectviscousgain(oi)*yvi;
    
                TrialData.X(trialsel,:) = xi;
                TrialData.Y(trialsel,:) = repmat(xout,size(t,1),1);
                TrialData.Force(trialsel,:) = xfi;
                TrialData.ForceIdeal(trialsel,:) = xfideali;
                if oi==3
                    TrialData.ForceIdealInterp(trialsel,:) = -outlierinterpviscousgain*yvi;
                end
            end

        end
    end
    return
end