%% plot avg pow over time window 

cfg.latency = [0 1];
cfg.avgovertime = 'yes';

cfg_ana.frequency = {[3 8],[8 12],[12 30]};
cfg.avgoverfreq = 'yes';

roi = 'noEyeABH';
roiind = ismember(ana.elecGroupsStr,roi);
cfg.channel = ana.elecGroups{roiind};
cfg.avgoverchan = 'yes';
avgdata = [];
for ifreq = 1:length(cfg_ana.frequency)
    
    cfg.frequency = cfg_ana.frequency{ifreq};
    conds = {'flckr0','flckr6','flckr10','flckr20'};
    %get avg data for each condition
    bardata = [];
    for icond = 1:length(conds)
        for isub = 1:length(exper.subjects)
            avgdata(isub).(conds{icond}) =  ft_selectdata(cfg, data_pow.ses1.(conds{icond}).sub(isub).data);
            bardata(isub,icond) = avgdata(isub).(conds{icond}).powspctrm;
        end
    end
    figure('color','white');
    hold on    
    bar(mean(bardata));
    plot(bardata','.','markersize',30);
    box off
    set(gca,'Xtick',1:4);
    set(gca,'Xticklabel',conds);
    ylabel(sprintf('%i-%i avg power',cfg.frequency(1), cfg.frequency(2)));
    set(gca,'fontsize',24);
end



