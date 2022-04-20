function [] = connect_and_sync_data()

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
    
    %%
    cd ~/Documents/wolpertlab/CurlFieldFamilies
    system('rsync -av /Volumes/eac2257/experiments/CurlFieldFamilies/Data/ ./Data')

end
