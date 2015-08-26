%calculate power for different specific trials
out = ent_behavior(dirs,exper);
ana = mm_ft_elecGroups(ana);

conds = {'flckr0','flckr6','flckr10','flckr20'};

%freqs = {[5.9 6.1] [9.9 10.1] [19.9 20.1]};
freqs = {[4 8] [8 12] [12 30]};
freqstr = {'rand','theta','alpha','beta'};

cfg = [];
%define time, freq, and channels of interest
cfg.avgoverfreq = 'yes';

cfg.latency = [0 1];
cfg.avgovertime = 'yes';

chanstr = 'PS2';
cfg.channel = ana.elecGroups{ismember(ana.elecGroupsStr,chanstr)};
cfg.avgoverchan = 'yes';

%get correct trials
cfg.trials = 'all';
cfg.avgoverrpt = 'no';
avgdata = [];
for isub = 1:length(exper.subjects)
    
    fprintf('\n%s\n------------\n',exper.subjects{isub});
    
    for icond = 1:length(conds)
        
        subdata = data_pow.ses1.(conds{icond}).sub(isub).data;
        
        for ifreq = 1:length(freqs)
            
            cfg.frequency = freqs{ifreq};
            
            tmp = ft_selectdata(cfg, subdata);
            ntrials = size(tmp.trialinfo,1);
            tmp = [tmp.trialinfo(:,[1 2 5 7 9]) repmat(tmp.time,ntrials,1) repmat(round(tmp.freq),ntrials,1) tmp.powspctrm];

            avgdata = cat(1,avgdata,tmp);
            
        end
        
    end
    
    fprintf('-----------\n');
    
end
vars = {'subno', 'blk', 'ordr', 'pfrq', 'crrct', 'avgtm', 'avgfrq', 'pow'};
ds = mat2dataset(avgdata,'VarNames',vars);


%%
% compare means

grpout = grpstats(ds,{'subno','pfrq','avgfrq','crrct'},{'mean','sem'},'datavars',{'pow'});
gdind = grpout.pfrq == grpout.avgfrq | grpout.pfrq==0;
grpout = grpout(gdind,:);


x1 = grpout.crrct==11 & grpout.avgfrq==6 & grpout.pfrq==6;
dsx1 = grpout(x1,:);
x1 = dsx1.mean_pow;
x1(end+1) = 0;
x0 = grpout.crrct==11 & grpout.avgfrq==6 & grpout.pfrq==0;
dsx0 = grpout(x0,:);
x0 = dsx0.mean_pow;
deltaX = x1-x0;

z1 = grpout.crrct==10 & grpout.avgfrq==6 & grpout.pfrq==6;
dsz1 = grpout(z1,:);
z1=dsz1.mean_pow;
z0 = grpout.crrct==10 & grpout.avgfrq==6 & grpout.pfrq==0;
dsz0 = grpout(z0,:);
z0=dsz0.mean_pow;
deltaZ = z1-z0;

