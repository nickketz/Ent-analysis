%% cluster statistics

files.figPrintFormat = 'fig';
files.saveFigs = 0;

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverfreq = 'yes';
%cfg_ft.avgoverfreq = 'no';

cfg_ft.parameter = 'powspctrm';

% debugging
%cfg_ft.numrandomization = 5000;

cfg_ft.numrandomization = 500;
cfg_ft.clusteralpha = .05;
cfg_ft.alpha = .1;
cfg_ft.neighbors = ft_prepare_neighbours(struct('elecfile',files.elecfile,'method','distance'));
%cfg_ft.minnbchan = 2;

cfg_ana = [];
cfg_ana.conditions = {'flckr20','flckr0'};

%cfg_ana.frequencies = [3 8; 8 12; 12 30; 30 50];
%cfg_ana.frequencies = [5.8 6.2; 9.8 10.2; 19.8 20.2];
cfg_ana.frequencies = [19.8 20.2];
cfg_ana.latencies = [0 1];
%cfg_ana.latencies = [0 0.5; 0.5 1.0];

for lat = 1:size(cfg_ana.latencies,1)
    cfg_ft.latency = cfg_ana.latencies(lat,:);
    for fr = 1:size(cfg_ana.frequencies,1)
        cfg_ft.frequency = cfg_ana.frequencies(fr,:);
        
        %[stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_pow);  
%         mydir = fullfile(dirs.saveDirProc,['tfr_stat_clus_' num2str(cfg_ft.latency(1)*1000) '_' num2str(cfg_ft.latency(2)*1000)]);
%         myfile = ['tfr_stat_clus_' cfg_ana.conditions{1} 'vs' cfg_ana.conditions{2} '_' sprintf('%.01f_%.01f',cfg_ft.frequency) '_' num2str(cfg_ft.latency(1)*1000) '_' num2str(cfg_ft.latency(2)*1000) '.mat'];
%         load(fullfile(mydir,myfile));
        
    end
end
cfg_plot.frequency = cfg_ft.frequency;
cfg_plot.latency = cfg_ft.latency;

%% plot the cluster statistics

%files.saveFigs = 0;

%cfg_ft = [];
%cfg_ft.alpha = .12;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;

% not averaging over frequencies - only works with ft_multiplotTFR
%files.saveFigs = 0;
%cfg_ft.avgoverfreq = 'yes';
%cfg_ft.avgoverfreq = 'no';
cfg_ft.interactive = 'yes';
cfg_plot.mask = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.maskstyle = 'saturation';
cfg_ft.maskalpha = 0.1;
cfg_plot.ftFxn = 'ft_multiplotTFR';
% http://mailman.science.ru.nl/pipermail/fieldtrip/2009-July/002288.html
% http://mailman.science.ru.nl/pipermail/fieldtrip/2010-November/003312.html

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  for fr = 1:size(cfg_plot.frequencies,1)
    cfg_ft.frequency = cfg_plot.frequencies(fr,:);
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs);

  end
end

%% make average plot for significant elecs

cfg_ft = [];
cfg_ft.alpha = .1;

cfg_ft.avgoverfreq = 'yes';
cfg_ft.avgovertime = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.mask = 'yes';
cfg_ft.highlightsizeseries  = repmat(15,6,1);

cfg_ft.maskstyle = 'opacity';
cfg_ft.transp = 1;  
%cfg_ft.maskstyle = 'saturation';
cfg_ft.maskalpha = 0.1;
cfg_ft.layout = ft_prepare_layout([],ana);
%cfg_ft.highlightseries ={'numbers','numbers','numbers','numbers','numbers','numbers'};
cfg_ft.highlightcolorpos = [1 1 1];
vsstr = fieldnames(stat_clus);

for iclus = 1:1%length(stat_clus.(vsstr{1}).posclusters)
    if stat_clus.(vsstr{1}).posclusters(iclus).prob < cfg_ft.alpha
        cfg_ft.clusnum = iclus;
        outdata = nk_ft_avgclustplot(stat_clus,cfg_plot,cfg_ft,dirs,files,1);        
        % plot cluster average power over time
        cfg_ft.time = [-1 1.5];       
        cfg_ft.conds = cat(2,cfg_ana.conditions);
        cfg_ft.conds = ana.eventValues{1};
        for icond = 1:length(cfg_ft.conds)
            cond1 = tokenize(cfg_ft.conds{icond},'_');
            cfg_plot.colors(icond,:) = rgb(cond1{2});
        end        
        outdata = nk_ft_avgpowerbytime(data_freq,stat_clus,cfg_plot,cfg_ft,dirs,files,0);
        %outdatafreq = nk_ft_avgpowerbyfreq(data_freq,stat_clus,cfg_plot,cfg_ft,dirs,files,1);
    else
        continue;
    end
end


%% descriptive statistics: ttest

%cfg_ana = [];
vsstr = fieldnames(stat_clus);
vsstr = vsstr{1};
ind = find(stat_clus.(vsstr).posclusterslabelmat==1);
[x,y,z] = ind2sub(size(stat_clus.(vsstr).posclusterslabelmat),ind);
% define which regions to average across for the test
%cfg_ana.rois = {{'PS'},{'FS'},{'LPS','RPS'},{'PS'},{'PS'}};
cfg_ana.rois = {stat_clus.(vsstr).label(unique(x))};
% define the times that correspond to each set of ROIs
%cfg_ana.latencies = [0.2 0.4; 0.6 1.0; 0.5 0.8; 0.5 1.0; 0.5 1.0];
cfg_ana.latencies = stat_clus.(vsstr).time(unique(z));
% define the frequencies that correspond to each set of ROIs
%cfg_ana.frequencies = [3 8; 3 8; 3 8; 3 8; 8 12];
cfg_ana.frequencies = stat_clus.(vsstr).cfg.frequency;

%cfg_plot.conditions = {{'THR','THF'},{'NTR','NTF'},{'THR','NTR'},{'THF','NTF'}};
cfg_ana.conditions = {'Word_Redbgt_Comb_Succ','Word_Green_Comb_Succ'};
%cfg_ana.conditions = {{'TH','NT'},{'TH','B'},{'NT','B'}};
%cfg_ana.conditions = {'all'};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.parameter = 'powspctrm';
cfg_ft.correctm = 'fdr';

cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
% line plot parameters
%cfg_plot.ylims = repmat([-1 1],size(cfg_ana.rois'));
cfg_plot.ylims = repmat([-100 100],size(cfg_ana.rois'));

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_freq);
end


