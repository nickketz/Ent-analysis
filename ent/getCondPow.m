%adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/test/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/pilot/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
load(adFile)

%%
%calculate power for different specific trials
out = ent_behavior(dirs,exper);
ana = mm_ft_elecGroups(ana);

conds = {'flckr0','flckr6','flckr10','flckr20'};
subnum = regexp(exper.subjects,'_([0-9]+)$','tokens');
subnum = cellfun(@(x) (x{1}),subnum);
subnum = cellfun(@(x) (str2double(x)),subnum);


freqs = {[5.9 6.1] [9.9 10.1] [19.9 20.1]}; frqbandstr = 'ent';
%freqs = {[4 8] [8 12] [12 30]}; frqbandstr = 'full';
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

%data loading details
cfg_load = [];
cfg_load.loadMethod = 'seg';
cfg_load.latency = 'all';
cfg_load.frequency = 'all';
cfg_load.keeptrials = 'yes';
cfg_load.equatetrials = 'no';
cfg_load.rmPreviousCfg = true;
cfg_load.ftype = 'pow';
cfg_load.output = 'pow';
cfg_load.transform = '';
cfg_load.norm_trials = 'single';
cfg_load.baseline_type = 'zscore';
%cfg_load.baseline_time = []; bslnstr = 'nobsln';
cfg_load.baseline_time = [-0.3 -0.1]; bslnstr = ['bsln' num2str(cfg_load.baseline_time(1),'%.01f') 'to' num2str(cfg_load.baseline_time(2),'%.01f')];
cfg_load.baseline_data = 'pow';
cfg_load.saveFile = false;
cfg_load.rmevoked = 'no';
cfg_load.rmevokedfourier = 'no';
cfg_load.rmevokedpow = 'no';


for isub = 1:length(exper.subjects)
    
    fprintf('\n%s\n------------\n',exper.subjects{isub});
    
    for icond = 1:length(conds)
        if exper.nTrials.(conds{icond})(isub)==0
            fprintf('no trials for %s %s, skipping\n',exper.subjects{isub},conds{icond});
            continue
        end
                
        
        %have to check for subsequent correct column(subcorr)
        myvars = {'subn','blkn','ordr','pfrq','subcorr'};
        colinds = [];
        for ivar = 1:length(myvars)
            colinds = cat(2,colinds,find(ismember(ana.trl_order.(conds{1}),myvars{ivar})));
        end


        %cull exper struct to specific subjects
        tmpexper = ent_rmSubs(exper,~ismember(exper.subjects,exper.subjects{isub}));
        %define trials of interest
        tmpana = ana;
        tmpana.eventValues = {conds(icond)};

        
        %load data
        [subdata,tmpexper] =  mm_ft_loadData_multiSes(cfg_load,tmpexper,dirs,tmpana);
        subdata = subdata.ses1.(conds{icond}).sub.data;
        if size(subdata.trialinfo,2)~=length(ana.trl_order.(conds{icond})) || mean(subdata.trialinfo(:,ismember(to.(conds{icond}),'subn'))) ~= subnum(isub)
            error('mismatch of trialinfo columns for %s %i\nprobably haven''t added subcorr to trialinfo, see getCorinfo.m',exper.subjects{isub},conds{icond});
        end
        
        for ifreq = 1:length(freqs)
            
            cfg.frequency = freqs{ifreq};
            
            tmp = ft_selectdata(cfg, subdata);
            ntrials = size(tmp.trialinfo,1);
            
            tmp = [tmp.trialinfo(:,colinds) repmat(tmp.time,ntrials,1) repmat(round(tmp.freq),ntrials,1) tmp.powspctrm];

            avgdata = cat(1,avgdata,tmp);
            
        end
        
    end
    
    fprintf('-----------\n');
    
end
vars = {'subno', 'blk', 'ordr', 'pfrq', 'crrct', 'avgtm', 'avgfrq', 'pow'};
ds = mat2dataset(avgdata,'VarNames',vars);

%change avgfreq to match pfrq in beta 
ind = ds.avgfrq==21;
ds.avgfrq(ind) = 20;

%%
%save data

fname = fullfile('/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/Ranalysis',sprintf('avgpow_n%i_%sband_%s_%s.txt',length(exper.subjects),frqbandstr,chanstr,bslnstr));
export(ds,'File',fname);

