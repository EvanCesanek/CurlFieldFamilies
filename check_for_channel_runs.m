
clearvars
WL = DynamicsExample;
cfg_name = 'OUT5';

%run

lol = T.block_original_trial >= 6;
cum = 0;
count = 0;
for di = lol'
    if di==0
        cum = 0;
    else
        cum = cum+1;
    end
    if cum==3
        count = count+1;
    end
end