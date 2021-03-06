clearvars
clc
colors = parula(6); colors = colors(1:5,:);
if ispc
    addpath('U:/gramm')
end

if contains(version(),'2017')
    timeline_marker = '+';
    correct_marker = '+';
else
    timeline_marker = '|';
    correct_marker = '_';
end

%% Single subject
%TrialData = analyze_subject(1,1,'eac_locstop2');
%TrialData = analyze_subject(1,1,'eac_locstopXL');
%TrialData = analyze_subject(1,1,'eac_locsliceTRI');
TrialData = analyze_subject(1,1,'eac_locstopTRI');

%TrialData = analyze_subject(1,1,'zz_locstop');

% subi = 0;
% for sub = {'eac_locstop2','eac_locstopTRI'}
%     subi = subi+1;
%     sesi = 1;
%     fn = sub{1};
%     fprintf('\nProcessing %s\n',fn);
%     td = analyze_subject(subi,sesi,fn);
%     if subi==1 && sesi==1
%         TrialData = td;
%     else
%         TrialData = [TrialData; td];
%     end
% end
%TrialData = analyze_subject(1,1,'eac_locstopbasic');
TrialData = analyze_subject(1,1,'zz_locstop');
%TrialData = analyze_subject(1,1,'eac_locstopCW2');
%TrialData = analyze_subject(1,1,'eac_pilotNEG2.mat');

%% Batch
% subi = 0;
% for sub = {'atf','ck','vp','ik','dg'}
%     subi = subi+1;
%     sesi = 0;
%     for ses = {'_s1.mat','_s2.mat','_s3.mat'}
%         sesi = sesi+1;
%         fn = [sub{1} ses{1}];
% 
%         fn2 = [];
%         if strcmpi(fn,'atf_s1.mat')
%             fn2 = 'atf_s1b.mat';
%         elseif strcmpi(fn,'vp_s2.mat')
%             fn2 = 'vp_s2b.mat';
%         end
% 
%         fprintf('\nProcessing %s\n',fn);
% 
%         td = analyze_subject(subi,sesi,fn,fn2);
%         if subi==1 && sesi==1
%             TrialData = td;
%             TrialData.Timed = ones(size(TrialData,1),1);
%         else
%             TrialData = [TrialData; td];
%         end
%     end
% end

subi = 0;
for sub = {'lr','ta'}
    subi = subi+1;
    sesi = 0;
    for ses = {'1.mat','2.mat'}
        sesi = sesi+1;
        fn = [sub{1} ses{1}];

        fn2 = [];

        fprintf('\nProcessing %s\n',fn);

        td = analyze_subject(subi,sesi,fn,fn2);
        if subi==1 && sesi==1
            TrialData = td;
        else
            TrialData = [TrialData; td];
        end
    end
end

%writetable(TrialData,'dat.csv');
%save('TrialDataN5.mat',"TrialData")
%clearvars
%load('TrialDataN5.mat')

%% Timeline
figure(100);
clf;
subset = ~TrialData.MissTrial & TrialData.Subject==1;
g1 = gramm('x',TrialData.TrialNumber(subset),'y',TrialData.TargetAngle(subset),'size',TrialData.ChannelTrial(subset),...
    'color',TrialData.ObjectId(subset),'row',TrialData.Session(subset));
g1.axe_property('YLim',[0 9*pi/8]);
g1.no_legend();
g1.set_point_options('markers',{timeline_marker});
g1.set_color_options('map',colors);
g1.geom_point();
g1.set_names('x','Block (WL)','y','Target Angle (rad.)','row','Session');
%g1.geom_abline('intercept',0,'slope',0,'style','k:');
g1.draw();

for mi = 1:length(get([g1.results.geom_point_handle],'MarkerFaceColor'))
    set(g1.results.geom_point_handle(mi),'MarkerEdgeColor',get(g1.results.geom_point_handle(mi),'MarkerFaceColor'),...
        'LineWidth',get(g1.results.geom_point_handle(mi),'MarkerSize')/10);
end

td = TrialData(subset & TrialData.block_trial==1, :);
g1.update('x',td.TrialNumber,'y',repelem(0,size(td,1)),'color',[],'size',[],'row',td.Session,'label',td.block_count);
g1.no_legend();
g1.geom_point();
g1.draw();

set([g1.results.geom_point_handle],'MarkerEdgeColor',[0 0 0],'MarkerSize',8,'LineWidth',1.5);

td = td(mod(td.block_count,5)==1, :);
g1.update('x',td.TrialNumber,'y',repelem(pi/8,size(td,1)),'color',[],'size',[],'row',td.Session,'label',td.block_count);
g1.no_legend();
g1.set_color_options('map',[0 0 0]);
g1.geom_label('HorizontalAlignment','center','VerticalAlignment','middle');
g1.draw();

%g1.export('file_name','MPE_N5.pdf','file_type','pdf');

%% MPE
figure(1);
clf;
% Set TrialData.TargetAngle~=pi/2 for MPE on field trials to peripheral targets

% Batch - before channel phase:
% subset = ~TrialData.MissTrial & ~TrialData.ChannelTrial & TrialData.TargetAngle==pi/2 & ...
%     (   (TrialData.Session==1 & ismember(TrialData.block_count,81:100)) | ...
%         (TrialData.Session==2 & ismember(TrialData.block_count,36:50)) | ...
%         (TrialData.Session==3 & ismember(TrialData.block_count,26:40))    );

% Batch - during channel phase:
% subset = ~TrialData.MissTrial & ~TrialData.ChannelTrial & TrialData.TargetAngle==pi/2 & ...
%     (   (TrialData.Session==1 & TrialData.block_count>100) | ...
%         (TrialData.Session==2 & TrialData.block_count>50) | ...
%         (TrialData.Session==3 & TrialData.block_count>40)    );

% Second batch (length objects)
subset = ~TrialData.MissTrial & ~TrialData.ChannelTrial & TrialData.TargetAngle==pi/2 & ... 
    (TrialData.Session==1 & ismember(TrialData.block_count,71:80));
    %(TrialData.Session==2 & ismember(TrialData.block_count,41:52));
    %(TrialData.Session==1 & ismember(TrialData.block_count,51:60));
%     (TrialData.Session==2 & ismember(TrialData.block_count,1:16));
% Additional - other periods:
subset = ~TrialData.MissTrial & ~TrialData.ChannelTrial & ...%TrialData.TargetAngle==pi/2 & ... 
    (TrialData.Session==1 & ismember(TrialData.block_count,61:80));

%     (TrialData.Session==1 & ismember(TrialData.block_count,1:16));
%     (TrialData.Session==1 & ismember(TrialData.block_count,46:60));
%     (TrialData.Session==2 & ismember(TrialData.block_count,41:52));

g1 = gramm('x',TrialData.ObjectId(subset),'y',TrialData.MPE(subset),'color',TrialData.ObjectId(subset),...
    'row',TrialData.Subject(subset),'column',TrialData.Session(subset));
g1.axe_property('YLim',[-2.5 5.5],'XLim',[0.5 5.5]);
g1.no_legend();
g1.set_color_options('map',colors);
g1.stat_summary('geom','bar','width',3);
g1.set_names('x','Object','y','MPE (cm)','row','Subject','column','Session');
g1.geom_abline('intercept',0,'slope',0,'style','k:');
g1.draw();

g1.update();
g1.set_color_options('map',colors*0.7);
g1.no_legend();
g1.stat_summary('geom','errorbar','width',1,'type','sem');
g1.draw();

%g1.export('file_name','MPE_PreChannel_N5.pdf','file_type','pdf');

%% Channel Trials - exclude outliers (grouped subject x object x target x session)
subset = ~TrialData.MissTrial & TrialData.ChannelTrial;
[groups,s,o,t,ss] = findgroups(TrialData.Subject(subset),TrialData.ObjectId(subset),TrialData.TargetAngle(subset), TrialData.Session(subset));
a = splitapply(@(x1){isoutlier(x1,'median','ThresholdFactor',5)},TrialData.CumForce(subset),groups);
for k = 1:length(a)
    tmp = TrialData.CumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k));
    tmp(a{k}) = NaN;
    TrialData.CumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k)) = tmp;

    tmp = TrialData.IdealCumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k));
    tmp(a{k}) = NaN;
    TrialData.IdealCumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k)) = tmp;

    tmp = TrialData.IdealInterpCumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k));
    tmp(a{k}) = NaN;
    TrialData.IdealInterpCumForce(subset & TrialData.Subject==s(k) & TrialData.ObjectId==o(k) & TrialData.TargetAngle==t(k) & TrialData.Session==ss(k)) = tmp;
end
fprintf('\n%i AF outlier trials excluded (%.2f%s of total).\n',sum(isnan(TrialData.CumForce(subset))),100*sum(isnan(TrialData.CumForce(subset)))/size(TrialData(subset,:),1),'%');

%% Channel Trial Integrated Forces
loc = 2; % 1 = outlier (NE), 2 = straight, all (N), 3 = largest (NW)
figure(2+loc);
clf;
g1 = gramm('x',TrialData.ObjectId(subset),'y',TrialData.FieldDirection(subset).*TrialData.CumForce(subset),'color',TrialData.ObjectId(subset),...
    'row',TrialData.Subject(subset),'column',TrialData.Session(subset),'subset',TrialData.TargetAngle(subset)==loc*pi/4);
g1.axe_property('YLim',[-2.5 2.5],'XLim',[0.5 5.5]);
g1.no_legend();
g1.set_color_options('map',colors);
g1.stat_summary('geom','bar','width',3);
g1.set_names('x','Object','y','Cum. Force (N)','row','Subject','column','Session');
g1.geom_abline('intercept',0,'slope',0,'style','k:');
g1.draw();

g1.update();
g1.set_color_options('map',colors*0.7);
g1.no_legend();
g1.stat_summary('geom','errorbar','width',1,'type','sem');
g1.draw();

b = splitapply(@nanmean, -TrialData.IdealInterpCumForce(subset), groups);

g1.update('x',o,'y',b,'color',[],'row',s,'column',ss,'subset',t==loc*pi/4);
g1.set_color_options('map',colors*0.75);
g1.no_legend();
g1.set_point_options('markers',{correct_marker});
g1.geom_point();
g1.geom_line();
g1.draw();

set([g1.results.geom_line_handle],'Color',[.6 .6 .6],'LineWidth',1,'LineStyle',':');
set([g1.results.geom_point_handle],'MarkerEdgeColor',[.6 .6 .6],'LineWidth',2);

c = splitapply(@nanmean, -TrialData.IdealCumForce(subset), groups);
sel = t==loc*pi/4 & o==3;
ndat = sum(sel);
nses = length(unique(ss(sel)));
nsub = length(unique(s(sel)));

if max(TrialData.Session)==3
    % Getting the outlier cum force on the right panels is tricky
    if loc==1 % the outlier location is easy because we don't need to skip any panels
        x = o;
        y = c;
        row = s;
        column = ss;
        subs = sel;
    elseif loc==2 % middle location, we need to skip the first column  (outlier is not in Session 1)
        x = repelem(3,ndat+1)';
        y = [nan; c(sel)];
        row = [1; repelem(1:nsub,2)'];
        column = [1; repmat([2 3]',nsub,1)];
        subs = [];
    elseif loc==3 % right location, we need to skip the first two columns (outlier is not in Sessions 1-2)
        x = repelem(3,ndat+2)';
        y = [nan; nan; c(sel)];
        row = [1; 1; repelem(1:nsub,1)'];
        column = [1; 2; repelem(3,nsub)'];
        subs = [];
    end

    g1.update('x',x,'y',y,'row',row,'column',column,'subset',subs);
    
else
    % Very simple when it's one subject one session, just one data point
    g1.update('x',repelem(3,nsub),'y',c(o==3),'row',s(o==3));
end

g1.set_color_options('map',colors*0);
g1.no_legend();
g1.set_point_options('markers',{correct_marker});
g1.geom_point();
g1.draw();

set([g1.results.geom_point_handle],'MarkerEdgeColor',colors(3,:)*0.5,'LineWidth',2);

%g1.export('file_name','CLF_N5.pdf','file_type','pdf');
