clearvars
clc

datadir = '~/Documents/wolpertlab/CurlFieldFamilies/Data/';
df = [];

subi = 0;
%for sub = {'atf','ck','vp','ik','dg'}
for sub = {'ik','dg'}
    subi = subi+1;
    sesi = 0;
    for ses = {'_s1.mat','_s2.mat','_s3.mat'}
        sesi = sesi+1;
        fn = [sub{1} ses{1}];
        clearvars -except sub ses fn datadir df subi sesi
        fn2 = [];
        if strcmpi(fn,'atf_s1.mat')
            fn2 = 'atf_s1b.mat';
        elseif strcmpi(fn,'vp_s2.mat')
            fn2 = 'vp_s2b.mat';
        end

        % filename, [fn2], datadir, [analyzetrajectories?], [thresh], [cuminterpframes], [barplotylims; ...] 
        %df = [df; analyze_subject(subi,sesi,fn,fn2,datadir,false,4,1:100,[100 700])];
        df = [df; analyze_subject(subi,sesi,fn,fn2,datadir,false,4,1:100)];
        %df = [df; analyze_subject(subi,sesi,fn,fn2,datadir,false,3,1:50,[0 750; 100 350])];
    end
end

%writetable(df,'dat.csv');

df = [df; analyze_subject(0,0,'dg_s3.mat',[],'~/Documents/wolpertlab/CurlFieldFamilies/Data/',true,Inf,1:100)];
