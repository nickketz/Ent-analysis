function clusplots(stat_clus,data_pow,cfg_plot,cfg_ana,dirs,exper,files,ana)

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
cfg_ft.maskalpha = 0.1;
cfg_ft.layout = ft_prepare_layout([],ana);
cfg_ft.highlightcolorpos = [1 1 1];
cfg_ft.highlightcolorneg = [1 1 1];

%plot topo 
nk_ft_avgclustplot(stat_clus,cfg_plot,cfg_ft,dirs,files,1);

if ~isfield(cfg_ana,'conds') 
    cfg_ft.conds = cat(2,cfg_ana.conditions); 
else
    cfg_ft.conds = cfg_ana.conds;
    %cfg_plot.colors = distinguishable_colors(length(cfg_ft.conds));
end

%plot avg over time
nk_ft_avgclus(data_pow.ses1,stat_clus,cfg_plot,cfg_ft,dirs,exper.badSub,files,1);


