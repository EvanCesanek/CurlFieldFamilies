clearvars
clc

if ~exist('/Volumes/eac2257/experiments/CurlFieldFamilies/Data','dir')
    [~, ssid] = system("/Sy*/L*/Priv*/Apple8*/V*/C*/R*/airport -I | grep -w SSID | awk '{print $2}'");
    if ~contains(ssid,'Columbia') % if we're not on columbia network open vpn app
        system('open -a "Cisco AnyConnect Secure Mobility Client"')
    end
end

%% once logged into vpn app, mount Engram
if ~exist('/Volumes/eac2257/experiments/CurlFieldFamilies/Data','dir')
    system('open smb://locker-smb.engram.rc.zi.columbia.edu/wolpert-locker/users/eac2257')
end

% if ~exist('/Volumes/wolpert-locker/users','dir')
%     system('open smb://locker-smb.engram.rc.zi.columbia.edu/wolpert-locker')
% end

%%
cd ~/Documents/wolpertlab/CurlFieldFamilies
system('rsync -av /Volumes/eac2257/experiments/CurlFieldFamilies/Data/ ./Data')

%% go to the data directory
%cd /Volumes/eac2257/experiments/CurlFieldFamilies/Data

%%
datadir = '~/Documents/wolpertlab/CurlFieldFamilies/Data/';

%fn = 'ec1.mat';
%fn = 'ec_out5.mat';
%fn = 'dw_out5.mat';
analyzetrain = true; % if analyzing out5, analyze early = true, else false
%fn = 'ec1_outiso.mat';
%fn = 'ec1_isof.mat';
%fn = 'ec1_isof5.mat';
%fn = 'ec2_iso5.mat';
%fn = 'ec4_iso5.mat';
%fn = 'ec2_iso3.mat';
%fn = 'ec2_iso3xl.mat';

%fn = 'zz_out2iso.mat';
%fn = 'jh_out2iso.mat';

%fn = 'ec_out2iso.mat';
%fn = 'ec2_out2iso.mat';
%fn = 'ec3_out2iso.mat';
%fn = 'ec4_out2iso.mat';
%fn = 'ec6_out2iso.mat';
%fn = 'ec_out2isol.mat';
%fn = 'ec3_out2isol.mat';
%fn = 'ec7_out2iso.mat';
%fn = 'ec8_out3iso.mat';
%fn = 'ec9_out3iso.mat';
%fn = 'ec10_out3iso.mat';
%fn = 'ec11_out3iso.mat';

%fn = 'ec_cptrain2.mat';
%fn = 'ec_cptest.mat';
%fn = 'ec1_cpout3iso.mat';
%fn = 'ec2_cpout3iso.mat';

%fn = 'ec_s1.mat';

for sub = {'atf','ck','vp'}
    for ses = {'_s1.mat','_s2.mat','_s3.mat'}
        fn = [sub ses];
        clearvars -except fn
        if strcmpi(fn,'atf_s1.mat')
            fn2 = 'atf_s1b.mat';
        elseif strcmpi(fn,'vp_s2.mat')
            fn2 = 'vp_s2b.mat';
        end

    %fn = 'atf_s1.mat';
    %fn2 = 'atf_s1b.mat';
    %fn = 'atf_s2.mat';
    %fn = 'atf_s3.mat';
    %fn = 'ck_s1.mat';
    %fn = 'ck_s2.mat';
    %fn = 'ck_s3.mat';
    %fn = 'vp_s1.mat';
    %fn = 'vp_s2.mat';
    %fn2 = 'vp_s2b.mat';
    %fn = 'vp_s3.mat';

load([datadir fn]);
if exist('fn2','var')
    TrialData_pre = TrialData;
    TimeStamp_pre = TimeStamp;
    sizepre = size(TimeStamp_pre);
    RobotPosition_pre = RobotPosition;
    RobotVelocity_pre = RobotVelocity;
    RobotForces_pre = RobotForces;

    load([datadir fn2]);
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
    TrialData = [TrialData_pre; TrialData];
    TimeStamp = [TimeStamp_pre; TimeStamp];
    RobotPosition = [RobotPosition_pre; RobotPosition];
    RobotVelocity = [RobotVelocity_pre; RobotVelocity];
    RobotForces = [RobotForces_pre; RobotForces];
end

%%
notmiss = ~TrialData.MissTrial;
channel = TrialData.FieldType==2;
expofam = strcmp(TrialData.block_name,'ExposureFam');
expoall = strcmp(TrialData.block_name,'ExposureAll');
chanall = strcmp(TrialData.block_name,'ChannelAll');

gentarg = strcmp(TrialData.block_name,'GeneralizeTarget');

famtest = strcmp(TrialData.block_name,'TestFam');
expotest = strcmp(TrialData.block_name,'ExposureTest');
gentest = strcmp(TrialData.block_name,'GeneralizeTest');

if ~any(strcmp('TargetAngle',TrialData.Properties.VariableNames))
    TrialData.TargetAngle = (pi/2)*ones(size(TrialData,1),1);
end

figh = figure(100);
clf;
colors = parula(6); colors = colors(1:5,:); colormap(colors);
%scatter(TrialData(notmiss,:).TrialNumber,TrialData.TargetAngle(notmiss,:),20*(1+str2num(num2str(TrialData(notmiss,:).FieldType==2))),TrialData(notmiss,:).ObjectId,'filled');
scatter(TrialData(notmiss & ~channel,:).TrialNumber,TrialData.TargetAngle(notmiss & ~channel,:),20*(1+str2num(num2str(TrialData(notmiss & ~channel,:).FieldType==2))),TrialData(notmiss & ~channel,:).ObjectId,'Marker','|');
hold on
scatter(TrialData(notmiss & channel,:).TrialNumber,TrialData.TargetAngle(notmiss & channel,:),20*(1+str2num(num2str(TrialData(notmiss & channel,:).FieldType==2))),TrialData(notmiss & channel,:).ObjectId,'Marker','|','LineWidth',2);
ylim([0,pi]);
exportgraphics(figh,[fn '_design.pdf']);

switch fn
    case 'ec1.mat'
        firstlateblock = 54;
        lateblock = TrialData.block_count >= firstlateblock;

    case 'ec_out5.mat'
        expolate = ismember(TrialData.block_count, 21:40);

    case 'dw_out5.mat'
        expolate = ismember(TrialData.block_count, 1:5);

end

angles = (pi/4):(pi/4):(3*pi/4); %unique(TrialData.TargetAngle)';
numangles = length(angles);
objects = unique(TrialData.ObjectId)';
numobjects = max(objects);

numdataframes = size(RobotPosition,3);
interplength = 100; % How many interpolated data points?
earlyinterpframes = 1:100; % Compute cumulative force over these frames
thresh = 5; % Outlier threshold (# MADs)

% initialize arrays
objectviscousgain = nan(numobjects,1);
targetviscousgain = nan(numangles,1);
outlierinterpviscousgain = 0.15;

multipleanglesperobject = false;

for ai = 1:length(angles)
    target = TrialData.TargetAngle==angles(ai); %abs(TrialData.TargetAngle-angles(ai))<0.0001;
    % viscous gain associated with the target during training (by virtue of one-to-one object-target pairings)
    % only relevant for OUT5 condition (should be identical to objectviscous gain in OUT1 condition)
    targetviscousgain(ai) = mean(unique(TrialData.FieldConstants(~channel & target,1)));

    for oi = objects
        object = TrialData.ObjectId==oi;
        % viscous gain associated with the object during training
        objectviscousgain(oi) = unique(TrialData.FieldConstants(~channel & object,1));

        switch fn
            case {'ec_cptrain2.mat'}
                ph = 3;
                trialsel = notmiss & channel & object & target & strcmpi(TrialData.block_short, 'W1') & TrialData.block_index_count>(ph*10) & TrialData.block_index_count<=(ph+1)*10;

            case {'ec_cptest.mat'}
                ph = 3;
                trialsel = notmiss & channel & object & target & strcmpi(TrialData.block_short, 'W3') & TrialData.block_index_count>(ph*10) & TrialData.block_index_count<=(ph+1)*10;

            case {'ec_outiso.mat','ec1_outiso.mat','ec1_isof.mat','ec1_isof5.mat', 'ec2_iso5.mat', 'ec4_iso5.mat', 'ec2_iso3.mat', 'ec2_iso3xl.mat'}
                leniso = 10; % number of trials in expotest phase
                oiso = regexp(fn,'\d*','Match');
                if isempty(oiso)
                    if any(strcmpi(fn,{'ec_outiso.mat','ec1_outiso.mat'}))
                        oiso = 3;
                    elseif strcmpi(fn,{'ec1_isof.mat'})
                        oiso = 2;
                    end
                else
                    oiso = str2double(oiso{numel(oiso)});
                end
                trialsel = notmiss & channel & object & target & (expotest | gentest);% & TrialData.block_count <= 35;
                multipleanglesperobject = true;
            
            case 'ec1.mat'
                trialsel = notmiss & channel & object & target & (expoall | chanall) & lateblock;
        
            case {'ec_out5.mat', 'dw_out5.mat'}
                if (analyzetrain)
                    trialsel = notmiss & channel & object & target & expolate;
                else
                    trialsel = notmiss & channel & object & target & gentarg;
                    multipleanglesperobject = true;
                end

            %case {'ec2_cpout3iso.mat', 'ec1_cpout3iso.mat', 'ec8_out3iso.mat', 'ec7_out2iso.mat', 'ec3_out2isol.mat', 'ec_out2isol.mat', 'ec6_out2iso.mat', 'ec4_out2iso.mat', 'ec3_out2iso.mat', 'ec2_out2iso.mat','ec_out2iso.mat','jh_out2iso.mat','zz_out2iso.mat'}
            otherwise
                leniso = 10; % number of trials in expotest phase
                %oiso = [3 5];
                trialsel = notmiss & channel & object & target & (famtest | expotest | gentest);% & TrialData.block_count <= 85;
                multipleanglesperobject = true;
                analyzetrajectories = false;
                if analyzetrajectories
                    %trialsel = notmiss & object & target & ~channel & TrialData.block_count >= 81; %& ismember(TrialData.block_count,21:25);
                    trialsel = notmiss & object & target & ~channel;
                end
        end
        
        if ~any(trialsel)
            continue % skip this angle if we didn't use it... shouldn't really happen except in pilot 'ec_out5.mat'
        end

        % keep track of the number of trials included in trialsel
        tmp = sum(trialsel);
        if ~exist('numtrialsincluded','var')
            numtrialsincluded = tmp;
            widearraysize = [150,numdataframes,numobjects,numangles]; % start with size(dim1) = 150 so we have plenty of room
            wideinterparraysize = [150,interplength,numobjects,numangles];
            % initialize arrays
            adaptindex = nan([150,numobjects,numangles]);
            
            t_wide = nan(widearraysize);
            x_wide = nan(widearraysize);
            y_wide = nan(widearraysize);
            yv_wide = nan(widearraysize);
            xf_wide = nan(widearraysize);
            xfideal_wide = nan(widearraysize);
            xfspatial_wide = nan(widearraysize);
            
            xi_wide = nan(wideinterparraysize);
            xfi_wide = nan(wideinterparraysize);
            xfideali_wide = nan(wideinterparraysize);
            xfspatiali_wide = nan(wideinterparraysize);
        end
        if tmp ~= numtrialsincluded
            warning('trialsel found %i trials for object %i angle %i, expected %i',tmp,oi,ai,numtrialsincluded);
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
        if size(x,2)==1 %numtrialsincluded==1
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
            rotmat = rotz(rotationangle);
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
        yendmove = 9.5;
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
        % set t=0 at first "ismoving" frame
        tcell = mat2cell(t,ones(size(t,1),1));
        t = cell2mat(cellfun(@(x) x-min(x), tcell, 'UniformOutput', false));

        t_wide(1:size(t,1),:,oi,ai) = t;
        x_wide(1:size(t,1),:,oi,ai) = x;
        y_wide(1:size(t,1),:,oi,ai) = y;
        yv_wide(1:size(t,1),:,oi,ai) = yv;
        xf_wide(1:size(t,1),:,oi,ai) = xf;

        % **** Analysis ****
    
        % Compute ideal force profile, based on gain and y-velocity profiles
        xfideal_wide(1:size(t,1),:,oi,ai) = -objectviscousgain(oi)*yv;
        for ti = 1:size(x,1)
            % compute adaptation index in each trial
            adaptindex(ti,oi,ai) = fitlm(array2table([xfideal_wide(ti,:,oi,ai)' xf(ti,:)']),'Var2 ~ Var1 - 1').Coefficients.Estimate;
        end
    
        % interpolate to fixed positions so we can average over profiles
        xout = linspace(ystartmove,yendmove,interplength);
        xfi = nan(size(x,1),length(xout));
        yvi = nan(size(x,1),length(xout));
        xi = nan(size(x,1),length(xout));
        yi = nan(size(x,1),length(xout));
        for ti = 1:size(x,1)
            yk = y(ti,:); % predictor ("x") variable is the y position
            [~,idx] = unique(yk); % predictor values must be unique
            idx = intersect(idx,find(~isnan(yk))); % ... and they can't be nans
            xfk = xf(ti,:); % forces against channel
            xfi(ti,:) = interp1(y(ti,idx)',xf(ti,idx)',xout);
            yvk = xf(ti,:); % forward velocities
            yvi(ti,:) = interp1(y(ti,idx)',yv(ti,idx)',xout);
            xk = x(ti,:);
            xi(ti,:) = interp1(y(ti,idx)',x(ti,idx)',xout);
        end
        % compute ideal interpolated force profile
        xfideali = -objectviscousgain(oi)*yvi;

        % store in wide format for easier averaging: trials x ypos x object 
        xfideali_wide(1:size(t,1),:,oi,ai) = xfideali; 
        xfi_wide(1:size(t,1),:,oi,ai) = xfi;
        xi_wide(1:size(t,1),:,oi,ai) = xi;

        if oi==3
            xfideali = -outlierinterpviscousgain*yvi;
            xfideali_wide(1:size(t,1),:,numobjects+1,ai) = xfideali;
        end
    end
end

%figure(999);clf;plot3(-xfideali_wide(:,:,3,2)',-xfi_wide(:,:,3,2)',xout');hold on;plot3([0 10],[0 10],[0 10],':k');hold off

%% Plot adaptation index
% for oi = 1:5
%     adaptindexiallmean(oi) = fitlm(array2table([mean(xfideali_wide(:,2:99,oi),1)' mean(xfidealiall(:,2:99),1)']),'Var2 ~ Var1 - 1').Coefficients.Estimate;
% end
% adaptmean = mean(adaptindex,1);
% for oi = 1:5
%     ci = bootci(10000,@mean,adaptindex(:,oi));
%     adaptci(oi,1) = ci(1)-adaptmean(oi);
%     adaptci(oi,2) = ci(2)-adaptmean(oi);
% end
% figh = figure(3);
% clf
% hold on
% for bi = 1:5
%     bar(bi,adaptmean(bi),'FaceColor',colors(bi,:));
%     errorbar(bi,adaptmean(bi),adaptci(bi,1),adaptci(bi,2),'.','LineWidth',2,'Color',colors(bi,:)*0.75);
%     scatter(bi,adaptindexiallmean(bi),'Marker','o','MarkerFaceColor',min(colors(bi,:)*1.33,[1 1 1]),'MarkerEdgeColor','k','LineWidth',1.2);
% end
% yline(1.00,'--','LineWidth',2);
% xlabel('Object ID');
% ylabel('Adaptation');


%% Plot cumulative force by movement direction (angle
figh = figure(4+10*multipleanglesperobject);
clf
clear axes
tlo = tiledlayout('flow','TileSpacing','compact');
for ai = 1:numangles
    if isempty(tlo.Children) || multipleanglesperobject
        % only need subplots if there are multiple angles per object combos
        axes(ai) = nexttile([10,1]);
    end
    hold on
    for oi = 1:numobjects
        if true || analyzetrajectories || ~ismember(oi,oiso) || (contains(fn,'out2iso') && ((ai==1 && oi==5) || (ai==3 && oi==3) || ai==2))
            xfi = xfi_wide(:,:,oi,ai);
            xfideali = xfideali_wide(:,:,oi,ai);
            %firstinterpframe = find(~isnan(xfi(1,:)),1);
            rowstokeep = find(~all(isnan(xfi),2));
            xfi = xfi(rowstokeep,:);
            xfideali = xfideali(rowstokeep,:);
            if isempty(xfi)
                continue
            end
    
            xficum{ai,oi} = sum(-xfi(:,earlyinterpframes),2,'omitnan');
            xficum_outidx{ai,oi} = find(~isoutlier(xficum{ai,oi},'ThresholdFactor',thresh));
            xfidealicum{ai,oi} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
            %xfidealicum(ai,oi) = cellfun(@(x,y) y(~ismaxormin(x)),xficum(ai,oi),xfidealicum(ai,oi),'UniformOutput',false);
            %xficum(ai,oi) = cellfun(@(x) x(~ismaxormin(x)),xficum(ai,oi),'UniformOutput',false);
            xfidealicum(ai,oi) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),xfidealicum(ai,oi),'UniformOutput',false);
            if oi==3
                xfideali = xfideali_wide(:,:,numobjects+1,ai);
                xfideali = xfideali(rowstokeep,:);
                xfidealicum{ai,numobjects+1} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
                xfidealicum(ai,numobjects+1) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),xfidealicum(ai,numobjects+1),'UniformOutput',false);
            end
            xficum(ai,oi) = cellfun(@(x) x(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),'UniformOutput',false);
            xficummean(ai,oi) = mean(xficum{ai,oi});
            xfidealicummean(ai,oi) = mean(xfidealicum{ai,oi});
            bar(oi,xficummean(ai,oi),0.5,'FaceColor',colors(oi,:));
            if oi==3
                xfidealicummean(ai,numobjects+1) = mean(xfidealicum{ai,numobjects+1});
            end
            if numtrialsincluded>1
                xficumci(ai,oi,:) = bootci(10000,@mean,xficum{ai,oi});
                xficumci(ai,oi,:) = xficumci(ai,oi,:)-xficummean(ai,oi);
                xfidealicumci(ai,oi,:) = bootci(10000,@mean,xfidealicum{ai,oi});
                xfidealicumci(ai,oi,:) = xfidealicumci(ai,oi,:)-xfidealicummean(ai,oi);
                errorbar(oi,xficummean(ai,oi),xficumci(ai,oi,1),xficumci(ai,oi,2),'.','LineWidth',2,'Color',colors(oi,:)*0.75);
                errorbar(oi,xfidealicummean(ai,oi),xfidealicumci(ai,oi,1),xfidealicumci(ai,oi,2),'Marker','none','MarkerFaceColor','black','MarkerEdgeColor','none','LineWidth',1.2,'Color',min(colors(oi,:)*1.33,[1 1 1]));
                if oi==3
                    xfidealicumci(ai,numobjects+1,:) = bootci(10000,@mean,xfidealicum{ai,numobjects+1});
                    xfidealicumci(ai,numobjects+1,:) = xfidealicumci(ai,numobjects+1,:)-xfidealicummean(ai,numobjects+1);
                    errorbar(oi,xfidealicummean(ai,numobjects+1),xfidealicumci(ai,numobjects+1,1),xfidealicumci(ai,numobjects+1,2),'Marker','none','MarkerFaceColor','black','MarkerEdgeColor','none','LineWidth',1.2,'Color',[0.5 0.5 0.5]);
                end
            end

        else
            xfi = xfi_wide(1:leniso,:,oi,ai);
            xfideali = xfideali_wide(1:leniso,:,oi,ai);
            %firstinterpframe = find(~isnan(xfi(1,:)),1);
            %rowstokeep = ~isnan(xfi(:,firstinterpframe));
            rowstokeep = find(~all(isnan(xfi),2));
            xfi = xfi(rowstokeep,:);
            xfideali = xfideali(rowstokeep,:);
            
            xficum{ai,numobjects+1} = sum(-xfi(:,earlyinterpframes),2,'omitnan');
            xfidealicum{ai,numobjects+1} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
            xfidealicum(ai,numobjects+1) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,numobjects+1),xfidealicum(ai,numobjects+1),'UniformOutput',false);
            xficum(ai,numobjects+1) = cellfun(@(x) x(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,numobjects+1),'UniformOutput',false);
            xficummean(ai,numobjects+1) = mean(xficum{ai,numobjects+1});
            xfidealicummean(ai,numobjects+1) = mean(xfidealicum{ai,numobjects+1});
            bar(oi-0.25,xficummean(ai,numobjects+1),0.5,'FaceColor',colors(oi,:));
            if numtrialsincluded>1
                xficumci(ai,numobjects+1,:) = bootci(10000,@mean,xficum{ai,numobjects+1});
                xficumci(ai,numobjects+1,:) = xficumci(ai,numobjects+1,:)-xficummean(ai,numobjects+1);
                xfidealicumci(ai,numobjects+1,:) = bootci(10000,@mean,xfidealicum{ai,numobjects+1});
                xfidealicumci(ai,numobjects+1,:) = xfidealicumci(ai,numobjects+1,:)-xfidealicummean(ai,numobjects+1);
                errorbar(oi-0.25,xficummean(ai,numobjects+1),xficumci(ai,numobjects+1,1),xficumci(ai,numobjects+1,2),'.','LineWidth',2,'Color',colors(oi,:)*0.75);
                errorbar(oi-0.25,xfidealicummean(ai,numobjects+1),xfidealicumci(ai,numobjects+1,1),xfidealicumci(ai,numobjects+1,2),'Marker','none','MarkerFaceColor','black','MarkerEdgeColor','none','LineWidth',1.2,'Color',min(colors(oi,:)*1.33,[1 1 1]));
            end

            xfi = xfi_wide(((leniso+1):end),:,oi,ai);
            xfideali = xfideali_wide(((leniso+1):end),:,oi,ai);
            %firstinterpframe = find(~isnan(xfi(1,:)),1);
            %rowstokeep = ~isnan(xfi(:,firstinterpframe));
            rowstokeep = find(~all(isnan(xfi),2));
            xfi = xfi(rowstokeep,:);
            xfideali = xfideali(rowstokeep,:);
            if isempty(xfi)
                continue
            end
    
            xficum{ai,oi} = sum(-xfi(:,earlyinterpframes),2,'omitnan');
            xfidealicum{ai,oi} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
            xfidealicum(ai,oi) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),xfidealicum(ai,oi),'UniformOutput',false);
            xficum(ai,oi) = cellfun(@(x) x(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),'UniformOutput',false);
            xficummean(ai,oi) = mean(xficum{ai,oi});
            xfidealicummean(ai,oi) = mean(xfidealicum{ai,oi});
            bar(oi+0.25,xficummean(ai,oi),0.5,'FaceColor',colors(oi,:));
            if numtrialsincluded>1
                xficumci(ai,oi,:) = bootci(10000,@mean,xficum{ai,oi});
                xficumci(ai,oi,:) = xficumci(ai,oi,:)-xficummean(ai,oi);
                xfidealicumci(ai,oi,:) = bootci(10000,@mean,xfidealicum{ai,oi});
                xfidealicumci(ai,oi,:) = xfidealicumci(ai,oi,:)-xfidealicummean(ai,oi);
                errorbar(oi+0.25,xficummean(ai,oi),xficumci(ai,oi,1),xficumci(ai,oi,2),'.','LineWidth',2,'Color',colors(oi,:)*0.75);
                errorbar(oi+0.25,xfidealicummean(ai,oi),xfidealicumci(ai,oi,1),xfidealicumci(ai,oi,2),'Marker','none','MarkerFaceColor','black','MarkerEdgeColor','none','LineWidth',1.2,'Color',min(colors(oi,:)*1.33,[1 1 1]));
            end

            
        end
    end
end
if multipleanglesperobject
    linkaxes(axes)
end
xlabel(tlo,'Object ID');
ylabel(tlo,'Cum. Lat. Force (N)');
if contains(fn,{'ec_cptrain','ec_cptest'})
    ylim([0 800]);
end
annotation('rectangle',[0 0 1 1],'Color','w');
exportgraphics(figh,[fn '_bar1.pdf']);


% Test for linear effect of size
try
    for li = 1:3
        lol = xficum(li,1:end);
        y = cell2mat(lol')';
        ym = cellfun(@mean,lol);
        Xm = 1:size(lol,2);
        X = [];
        for oi = 1:size(lol,2)
            X = [X; oi*ones(size(lol{oi}))];
        end
        res = fitlm(X,y);
        res = fitlm(Xm,ym);
        p(li) = res.Coefficients.pValue(2);
    end
catch

end

%% Simulation: outlier performance only affects intercept, not slope
% X = repelem(1:5,10)';
% y = X+randn(size(X));
% fitlm(X,y).Coefficients.Estimate'
% for i = 1:100000
%     y(X==3) = y(X==3)+randn(size(X(X==3)));
% end
% fitlm(X,y).Coefficients.Estimate'

%% Plot cumulative force by movement direction (angle
figh = figure(4+20*multipleanglesperobject);
clf
clear axes
tlo = tiledlayout('flow','TileSpacing','compact');
for ai = 1:numangles
    if isempty(tlo.Children) || multipleanglesperobject
        % only need subplots if there are multiple angles per object combos
        axes(ai) = nexttile([10,1]);
    end
    hold on
    for oi = 1:numobjects
        xfi = xfi_wide(:,:,oi,ai);
        xfideali = xfideali_wide(:,:,oi,ai);
        %firstinterpframe = find(~isnan(xfi(1,:)),1);
        rowstokeep = find(~all(isnan(xfi),2));
        xfi = xfi(rowstokeep,:);
        xfideali = xfideali(rowstokeep,:);
        if isempty(xfi)
            continue
        end

        xficum{ai,oi} = sum(-xfi(:,earlyinterpframes),2,'omitnan');
        xfidealicum{ai,oi} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
        %xfidealicum(ai,oi) = cellfun(@(x,y) y(~ismaxormin(x)),xficum(ai,oi),xfidealicum(ai,oi),'UniformOutput',false);
        %xficum(ai,oi) = cellfun(@(x) x(~ismaxormin(x)),xficum(ai,oi),'UniformOutput',false);
        xfidealicum(ai,oi) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),xfidealicum(ai,oi),'UniformOutput',false);
        if oi==3
            xfideali = xfideali_wide(:,:,numobjects+1,ai);
            xfideali = xfideali(rowstokeep,:);
            xfidealicum{ai,numobjects+1} = sum(-xfideali(:,earlyinterpframes),2,'omitnan');
            xfidealicum(ai,numobjects+1) = cellfun(@(x,y) y(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),xfidealicum(ai,numobjects+1),'UniformOutput',false);
        end
        xficum(ai,oi) = cellfun(@(x) x(~isoutlier(x,'ThresholdFactor',thresh)),xficum(ai,oi),'UniformOutput',false);
        xferricum{ai,oi} = xficum{ai,oi}-xfidealicum{ai,oi};
        xferricummean(ai,oi) = mean(xferricum{ai,oi});
        bar(oi,xferricummean(ai,oi),0.5,'FaceColor',colors(oi,:));
        if oi==3
            xferricum{ai,numobjects+1} = xficum{ai,oi}-xfidealicum{ai,numobjects+1};
            xferricummean(ai,numobjects+1) = mean(xferricum{ai,numobjects+1});
            bar(oi,xferricummean(ai,numobjects+1),0.5,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);%colors(oi,:));
        end
        if numtrialsincluded>1
            xferricumci(ai,oi,:) = bootci(10000,@mean,xferricum{ai,oi});
            xferricumci(ai,oi,:) = xferricumci(ai,oi,:)-xferricummean(ai,oi);
            errorbar(oi,xferricummean(ai,oi),xferricumci(ai,oi,1),xferricumci(ai,oi,2),'.','LineWidth',2,'Color',colors(oi,:)*0.75);
            if oi==3
                xferricumci(ai,numobjects+1,:) = bootci(10000,@mean,xferricum{ai,numobjects+1});
                xferricumci(ai,numobjects+1,:) = xferricumci(ai,numobjects+1,:)-xferricummean(ai,numobjects+1);
                errorbar(oi,xferricummean(ai,numobjects+1),xferricumci(ai,numobjects+1,1),xferricumci(ai,numobjects+1,2),'.','LineWidth',2,'Color',[0.5 0.5 0.5]);%colors(oi,:)*0.75);
            end
        end
    end
end
if multipleanglesperobject
    linkaxes(axes)
end
xlabel(tlo,'Object ID');
ylabel(tlo,'Cum Lat. Force Error (N)');
if contains(fn,{'ec_cptrain','ec_cptest'})
    ylim([0 800]);
end
annotation('rectangle',[0 0 1 1],'Color','w');
exportgraphics(figh,[fn '_bar2.pdf']);


%%
if false && multipleanglesperobject % only need this plot if there are multiple angles-per-object combos
    figure(5+10*multipleanglesperobject)
    clf
    clear axes
    tlo = tiledlayout('flow');
    for oi = 1:numobjects
        axes(oi) = nexttile([10,1]);
        hold on
        for ai = 1:numangles
            bar(ai,xficummean(ai,oi),'FaceColor',colors(oi,:));
            if numtrialsincluded>1
                errorbar(ai,xficummean(ai,oi),xficumci(ai,oi,1),xficumci(ai,oi,2),'.','LineWidth',2,'Color',colors(oi,:)*0.75);
                errorbar(ai,xfidealicummean(ai,oi),xfidealicumci(ai,oi,1),xfidealicumci(ai,oi,2),'Marker','o','MarkerFaceColor','black','MarkerEdgeColor','none','LineWidth',1.2,'Color',min(colors(oi,:)*1.33,[1 1 1]));
            end
        end
    end
    xlabel(tlo,'Target ID');
    ylabel(tlo,'Cum. Lat. Force (N)');
    linkaxes(axes);
end

%% Plot individual trajectories
figh = figure(1+10*multipleanglesperobject);
clf
tlo = tiledlayout(numobjects, numangles);
for oi = 1:numobjects
    for ai = 1:numangles
        axes(sub2ind([numobjects,numangles],oi,ai)) = nexttile;
        hold on
        
        % to remove outlier trials
        xfi = xfi_wide(:,:,oi,ai);
        rowstokeep = find(~all(isnan(xfi),2));
        rowstokeep = rowstokeep(ismember(rowstokeep,xficum_outidx{ai,oi}));
        
        plot(xfideal_wide(rowstokeep,:,oi,ai)', y_wide(rowstokeep,:,oi,ai)', '--');
        plot(xf_wide(rowstokeep,:,oi,ai)', y_wide(rowstokeep,:,oi,ai)', '-')
        xlim([-5,5]);
        axis auto
    end
end
linkaxes(axes);
xlabel(tlo,sprintf('locations'),'FontSize',16);
ylabel(tlo,sprintf('objects'),'FontSize',16);
corner = nexttile(sub2ind([numangles numobjects],1,numobjects));
xlabel(corner,'force','FontSize',8);
ylabel(corner,'y','FontSize',8)

exportgraphics(figh,[fn '_traj.pdf']);
%% Plot average trajectories
figh = figure(2+10*multipleanglesperobject);
clf
clear axes
tlo = tiledlayout('flow');
for oi = 1:numobjects
    for ai = 1:numangles
        xfi = xfi_wide(:,:,oi,ai);
        xfideali = xfideali_wide(:,:,oi,ai);
        %firstinterpframe = find(~isnan(xfi(1,:)),1);
        %rowstokeep = ~isnan(xfi(:,firstinterpframe));
        rowstokeep = find(~all(isnan(xfi),2));
        
        % to remove outlier trials
        rowstokeep = rowstokeep(ismember(rowstokeep,xficum_outidx{ai,oi}));
        
        xfi = xfi(rowstokeep,:);
        xfideali = xfideali(rowstokeep,:);
        if isempty(tlo.Children) || multipleanglesperobject
            % only need subplots if there are multiple angles per object combos
            axes(sub2ind([numobjects,numangles],oi,ai)) = nexttile;
        end
        plot(mean(xfideali,1,'omitnan'),xout,'LineWidth',1.5,'LineStyle','--','Color',colors(oi,:));
        hold on
        se = std(xfi,0,1)/sqrt(size(xfi,1));
        m = mean(xfi,1,'omitnan');
        xout_omitnan = xout(~isnan(m));
        se = se(~isnan(m));
        m = m(~isnan(m));
        fill([m-se fliplr(m+se)], [xout_omitnan fliplr(xout_omitnan)], colors(oi,:), 'LineStyle','none', 'FaceAlpha',0.4);
        plot(m,xout_omitnan,'LineWidth',2,'Color',colors(oi,:));
        %axis equal
        %xfidealiall = reshape(permute(xfideali_wide,[2 1 3]),size(xfideali_wide,2),[])';
        %plot(mean(xfidealiall(:,2:99),1),xout,'LineWidth',1.75,'LineStyle','-.','Color',[0 0 0]);
    end
end
xlabel(tlo,'Lateral Force (N)');
ylabel(tlo,'Y Position');
if multipleanglesperobject
    linkaxes(axes);
end
exportgraphics(figh,[fn '_trajmean.pdf']);

%% Plot average trajectories
if analyzetrajectories
    figure(3+10*multipleanglesperobject);
    clf
    clear axes
    tlo = tiledlayout('flow');
    for oi = 1:numobjects
        for ai = 1:numangles
            xfi = xi_wide(:,:,oi,ai);
            %firstinterpframe = find(~isnan(xfi(1,:)),1);
            %rowstokeep = ~isnan(xfi(:,firstinterpframe));
            rowstokeep = find(~all(isnan(xfi),2));
            xfi = xfi(rowstokeep,:);
            if isempty(tlo.Children) || multipleanglesperobject
                % only need subplots if there are multiple angles per object combos 
                axes(sub2ind([numobjects,numangles],oi,ai)) = nexttile;
            end
            if isempty(xfi)
                continue
            end
            hold on
            se = std(xfi,0,1)/sqrt(size(xfi,1));
            m = mean(xfi,1,'omitnan');
            xout_omitnan = xout(~isnan(m));
            se = se(~isnan(m));
            m = m(~isnan(m));
            plot(xfi',xout,'LineWidth',1,'Color',[colors(oi,:) 0.3]);
            %fill([m-se fliplr(m+se)], [xout_omitnan fliplr(xout_omitnan)], colors(oi,:), 'LineStyle','none', 'FaceAlpha',0.4);
            plot(m,xout_omitnan,'LineWidth',2,'Color',colors(oi,:));
            if contains(fn,{'out3iso','_s'})
                plot([-0.33 0.33],[10 10],'k-','LineWidth',3)
            else
                plot([-0.5 0.5],[10 10],'k-','LineWidth',3)
            end
            %axis equal
            %xfidealiall = reshape(permute(xfideali_wide,[2 1 3]),size(xfideali_wide,2),[])';
            %plot(mean(xfidealiall(:,2:99),1),xout,'LineWidth',1.75,'LineStyle','-.','Color',[0 0 0]);
            xlim([-2 2]);
        end
    end
    xlabel(tlo,'X Position');
    ylabel(tlo,'Y Position');
    if multipleanglesperobject
        linkaxes(axes);
    end
    exportgraphics(figh,[fn '_trajfield.pdf']);
end


%%
% close all
% figh = figure(99);
% clf
% 
% object = TrialData.ObjectId==3;
% 
% trialsel = notmiss & channel & object & chanall;
% 
% x = squeeze(RobotPosition(trialsel, 1, :));
% y = squeeze(RobotPosition(trialsel, 2, :));
% x(State(trialsel,:)~=9) = nan;
% y(State(trialsel,:)~=9) = nan;
% 
% xf = squeeze(RobotForces(trialsel, 1, :));
% yf = squeeze(RobotForces(trialsel, 2, :));
% 
% yv = squeeze(RobotVelocity(trialsel, 1, :));
% 
% t = TimeStamp(trialsel, :);
% 
% ynear = y > 9;
% ynearcell = mat2cell(ynear,ones(size(ynear,1),1));
% endframe = cellfun(@(x) find(x,1), ynearcell);
% for i = 1:size(y,1)
%     x(i,endframe(i):end) = nan;
%     y(i,endframe(i):end) = nan;
%     xf(i,endframe(i):end) = nan;
%     yf(i,endframe(i):end) = nan;
%     yv(i,endframe(i):end) = nan;
%     t(i,endframe(i):end) = nan;
% end
% 
% %plot(x',y','Color','blue');
% %axis equal
% plot(t',yv','Color','blue');

function out = ismaxormin(x)
    out = zeros(size(x));
    [~,maxi] = max(x);
    [~,mini] = min(x);
    out(maxi) = 1;
    out(mini) = 1;
end