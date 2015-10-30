%adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/test/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/pilot/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
load(adFile)

%%
%calculate power for different specific trials
out = ent_behavior(dirs,exper);
ana = mm_ft_elecGroups(ana);

conds = {'flckr6','flckr10','flckr20'};
subnum = regexp(exper.subjects,'_([0-9]+)$','tokens');
subnum = cellfun(@(x) (x{1}),subnum);
subnum = cellfun(@(x) (str2double(x)),subnum);

clustfreqs = {[5.8 6.2] [9.8 10.2] [19.8 20.2]}; 
%freqs = {[5.8 6.2] [9.8 10.2] [19.8 20.2]}; frqbandstr = 'ent';
freqs = {[4 8] [8 12] [12 30]}; frqbandstr = 'full';
freqstr = {'theta','alpha','beta'};

cfg = [];
%define time, freq, and channels of interest
cfg.avgoverfreq = 'yes';

cfg.latency = [0 1];
cfg.avgovertime = 'no';

%chanstr = 'PS2';
%cfg.channel = ana.elecGroups{ismember(ana.elecGroupsStr,chanstr)};
cfg.avgoverchan = 'no';

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
cfg_load.baseline_time = []; bslnstr = 'nobsln';
%cfg_load.baseline_time = [-0.3 -0.1]; bslnstr = ['bsln' num2str(cfg_load.baseline_time(1),'%.01f') 'to' num2str(cfg_load.baseline_time(2),'%.01f')];
cfg_load.baseline_data = 'pow';
cfg_load.saveFile = false;
cfg_load.rmevoked = 'no';
cfg_load.rmevokedfourier = 'no';
cfg_load.rmevokedpow = 'no';

cfg_clust.latency = [0 1];

for isub = 1:length(exper.subjects)
    
    fprintf('\n%s\n------------\n',exper.subjects{isub});
    
    for icond = 1:length(conds)
        if exper.nTrials.(conds{icond})(isub)==0
            fprintf('no trials for %s %s, skipping\n',exper.subjects{isub},conds{icond});
            continue
        end
        
        %load cluster
        % set the directory to load the file from
        dirClusStat = fullfile(dirs.saveDirProc,sprintf('tfr_stat_clus_%d_%d',round(cfg_clust.latency(1)*1000),round(cfg_clust.latency(2)*1000)));
        vs_str = sprintf('%svs%s',conds{icond},'flckr0');
        fprintf('\nLoading new cluster %s, %d--%d ms, %.1f--%.1f Hz\n',vs_str,round(cfg_clust.latency(1)*1000),round(cfg_clust.latency(2)*1000),clustfreqs{icond}(1),clustfreqs{icond}(2));
        
        savedFile = fullfile(dirClusStat,sprintf('tfr_stat_clus_%s_%.1f_%.1f_%d_%d.mat',vs_str,...
            clustfreqs{icond}(1),clustfreqs{icond}(2),...
            round(cfg_clust.latency(1)*1000),round(cfg_clust.latency(2)*1000)));
        
        if exist(savedFile,'file')
            fprintf('Loading %s\n',savedFile);
            load(savedFile);
        else
            error('cluster file missing');
        end
        clustind = stat_clus.(vs_str).posclusterslabelmat==1;
        
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
        
        
        
        for ibsln = 1:2
            if ibsln == 1
                fprintf('\nExtracting condition %s average power from cluster\n',conds{icond});
                tmpana.eventValues = {conds(icond)};
                %load data
                [subdata,tmpexper] =  mm_ft_loadData_multiSes(cfg_load,tmpexper,dirs,tmpana);
                subdata = subdata.ses1.(conds{icond}).sub.data;
                
                if size(subdata.trialinfo,2)~=length(ana.trl_order.(conds{icond})) || mean(subdata.trialinfo(:,ismember(to.(conds{icond}),'subn'))) ~= subnum(isub)
                    error('mismatch of trialinfo columns for %s %i\nprobably haven''t added subcorr to trialinfo, see getCorinfo.m',exper.subjects{isub},conds{icond});
                end
            else
                fprintf('\nExtracting condition %s average power from cluster\n','flckr0');
                tmpana.eventValues = {{'flckr0'}};
                [subdata,tmpexper] =  mm_ft_loadData_multiSes(cfg_load,tmpexper,dirs,tmpana);
                subdata = subdata.ses1.flckr0.sub.data;

                if size(subdata.trialinfo,2)~=length(ana.trl_order.(conds{icond})) || mean(subdata.trialinfo(:,ismember(to.(conds{icond}),'subn'))) ~= subnum(isub)
                    error('mismatch of trialinfo columns for %s %i\nprobably haven''t added subcorr to trialinfo, see getCorinfo.m',exper.subjects{isub},conds{icond});
                end
            end
            
            
            for ifreq = 1:length(freqs)
                
                cfg.frequency = freqs{ifreq};
                
                tmp = ft_selectdata(cfg, subdata);
                tmpclust = nan(size(tmp.powspctrm));
                for i = 1:size(tmp.powspctrm,1), tmpclust(i,:,:,:) = clustind; end
                psclust = tmp.powspctrm(logical(tmpclust));
                psclust = reshape(psclust,size(tmp.trialinfo,1),sum(reshape(clustind,numel(clustind),1)));
                tmp.powspctrm = nanmean(psclust,2);
                ntrials = size(tmp.trialinfo,1);
                
                tmp = [tmp.trialinfo(:,colinds) repmat(mean(cfg.latency),ntrials,1) repmat(round(mean(freqs{ifreq})),ntrials,1) tmp.powspctrm];
                
                avgdata = cat(1,avgdata,tmp);
                
            end
            
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

fname = fullfile('/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/Ranalysis',sprintf('avgpow_n%i_%sband_%s_%s.txt',length(exper.subjects),frqbandstr,'clust',bslnstr));
export(ds,'File',fname);

