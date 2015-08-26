function out = ent_behavior(dirs,exper,cfg)
%Behavioral analysis of Entrain experiment
%
% input
%   cfg: config structure
%    -filefiltstr: regexp filter on file listing, default is 'EntAssoc_*_Sub[0-9]{1,2}\.txt'
%    -blocks: array of which blocks to break down output by, output will be a
%    strucutred array with stats calculated for each of the blocks in the
%    array
%    -precision: precision of Freq to adhere to when assigning stims
%    to condition bins, i.e. 10 rounds observed pctPeriod values to the
%    nearest 10th of a percent.
%    -doplots: bool to do plots or not.
%    -badsubs: subject numbers to exclude from analysis
%
%
% output
%   out: output structure
%    -data: data structure from entReadLogs
%    -results: results of analysis
%

%read logs
if ~exist('cfg','var')
    cfg.filefiltstr = 'EntAssocWordFreq_12blks_1reps_16trials_Sub[0-9]{1,2}\.txt';
else
    if ~isfield(cfg,'filefiltstr')
        cfg.filefiltstr = 'EntAssocWordFreq_12blks_1reps_16trials_Sub[0-9]{1,2}\.txt';
    end
end


logdata = entReadlogs(dirs,cfg);

%verify subject number matches logdata
snums = logdata.data.Subj(1,:);
lnums = regexp(exper.subjects,'_([0-9]+)$','tokens');
for isub = 1:length(snums)
    if snums(isub) ~= str2double(lnums{isub}{1})
        error('subject number logs don''t match exper data');
    end
end

if ~isfield(cfg,'doplots')
    cfg.doplots = 0;
end

%block indicies
nblks = max(max(logdata.data.Block));
indblk = logdata.data.Block(:,1) == 1;
for iblk = 2:nblks
    indblk = cat(2,indblk, logdata.data.Block(:,1) == iblk);
end

if isfield(cfg,'block')
    if strcmp(cfg.block,'all') & ~iscell(cfg.block)
        cfg.block = {cfg.block};
    end
    origdata = logdata;
    blkdata = getBlkData(logdata);
else
    cfg.block = {'all'};
end

for iblk = cfg.block
    
    if ~strcmp(iblk,'all')
        logdata.data = blkdata{iblk};
    else
        iblk = 1;
        if isfield(cfg,'badblocks')
            tmp = rmBlkData(logdata,cfg.badblocks);
            logdata.data = tmp;
        end
    end
    
    %phase indicies
    iphase2 = logdata.data.Phase(:,1) == 2;
    iphase3 = logdata.data.Phase(:,1) == 3;
    iphase1 = logdata.data.Phase(:,1) == 1;
    
    p1acc = mean(logdata.data.Correct(iphase1,:));
    p2acc = mean(logdata.data.Correct(iphase2,:));
    
    %create indicies
    %resp don't know to old/new question
    rdknow = ~cellfun(@isempty,regexp(logdata.data.Answer,'^D'));
    %resp don't know to pair question
    rpdknow = ~cellfun(@isempty,regexp(logdata.data.Answer,'OD'));    
    %resp old
    rold = ~cellfun(@isempty,regexp(logdata.data.Answer,'^O[.]*'));
    %old trial
    iold = ~cellfun(@isempty,regexp(logdata.data.StimType,'Old'));
    iold = iold & ~rdknow;
    %old correct
    cold = logdata.data.Correct>1;
    %old correct with pair correct
    cpold = logdata.data.Correct>10;
    
    %resp new
    rnew = ~cellfun(@isempty,regexp(logdata.data.Answer,'^N'));
    %new trial
    inew = ~cellfun(@isempty,regexp(logdata.data.StimType,'New'));
    inew = inew & ~rdknow;
    %new correct
    cnew = logdata.data.Correct == 1 & repmat(iphase3,1,size(logdata.data.Correct,2));
    

    
    
    %determine condition indicies from data
    if isfield(cfg,'precision')
        precision = cfg.precision;
    else
        precision = 10;
    end
    conds = unique(round(precision*logdata.data.Freq(iphase2,:)))/precision;
    conds = conds(~isnan(conds));
    indConds = {};
    for icond = 1:length(conds)
        indConds{icond} = logdata.data.Freq==conds(icond);
    end

    
    %dprime
    hit = sum(rold & iold)./sum(iold);
    hit(hit==1) = .99;
    miss = sum(rnew & iold)./sum(iold);
    cr = sum(rnew & inew)./sum(inew);
    fa = sum(rold & inew)./sum(inew);
    fa(fa==0) = 1/(2*sum(inew(:,1)));
    dprime = norminv(hit)-norminv(fa);
    
    %pair hit rate
    phit = sum(rold & cpold)./sum(rold & iold & ~rpdknow);
    
    %dprime by cond
    chit = cell(size(conds)); cmiss = chit; ccr = chit; cfa = chit; cdprime = chit; cphit = chit;
    for icond = 1:length(conds)
        chit{icond} = sum(rold & iold & indConds{icond}) ./sum(iold & indConds{icond});
        chit{icond}(chit{icond}==1) = .99; %hack to avoid infinite dprime
        cmiss{icond} = sum(rnew & iold & indConds{icond})./sum(iold & indConds{icond});
        ccr{icond} = sum(rnew & inew) ./sum(inew);
        cfa{icond} = sum(rold & inew) ./sum(inew);
        cfa{icond}(cfa{icond}==0) = 1/(2*sum(inew(:,1)));%hack to avoid infinte dprime
        cdprime{icond} = norminv(chit{icond})-norminv(cfa{icond});
        %pair acc (norm based on responded 'old')
        cphit{icond} = sum(cpold & indConds{icond} & rold)./sum(rold & ~rpdknow & indConds{icond}); %discards 'don't know' trials
%        cphit{icond} = sum(cpold & indConds{icond} & rold)./sum(rold & indConds{icond}); %considers 'don't know' trials as incorrect
    end
    
    
    %response bias
    tmp = repmat(iphase3,1,size(logdata.data.Freq,2));
    pctold = sum(rold)./sum(tmp);
    pctnew = sum(rnew)./sum(tmp);
    pctdknow = sum(rdknow&tmp)./sum(tmp);
    pctpdknow = sum(rpdknow)./sum(rold);
    
    
    if cfg.doplots
        %dprime across conditions
        figure('color','white');
        mycolors = get(gca,'ColorOrder');
        plot([hit',fa'], 's', 'markersize', 20);
        hold on
        plot(dprime,'k.','markersize',30);
        %plot(0:length(dprime),repmat(mean(dprime),length(dprime)+1),'--k','linewidth',2);
        shadedErrorBar(0:length(dprime)+1,repmat(mean(dprime),1,length(dprime)+2),repmat(ste(dprime'),1,length(dprime)+2),{'--k','linewidth',2},1);
        shadedErrorBar(0:length(hit)+1,repmat(mean(hit),1,length(hit)+2),repmat(ste(hit'),1,length(hit)+2),{'--','linewidth',2,'color',mycolors(1,:),'markerfacecolor',mycolors(1,:)},1);
        shadedErrorBar(0:length(fa)+1,repmat(mean(fa),1,length(fa)+2),repmat(ste(fa'),1,length(fa)+2),{'--','linewidth',2,'color',mycolors(2,:),'markerfacecolor',mycolors(2,:)},1);
        xlabel('subject number','fontsize',18);
        legend({'hit-rate','fa-rate','dprime'},'fontsize',18,'location','best');
        set(gca,'fontsize',18);
        box off
        
        
        %dprime by condition, works with multiple conditions
        figure('color','white');
        hold on
        tmp = cell2mat(cdprime)';
        plot(tmp,'.','markersize',30);
        ylabel('dprime','fontsize',18);
        xlabel('subject number','fontsize',18);
        mycolors = get(gca,'ColorOrder');
        myconds = conds;
        for icond = 1:length(myconds)
            shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
        end
        legend({num2str(myconds)},'fontsize',18,'location','best');
        set(gca,'fontsize',18);
        box off
        hold off
        
        
        %pair-hit across and within conditions
        %tmp = cphit;
        %pairhit = reshape(cell2mat(tmp),[size(tmp{1},2) size(tmp,1)]);
        %tmp = pairhit;
        tmp = cell2mat(cphit)';
        figure('color','white');
        hold on
        plot(tmp,'.','markersize',30);
%        plot(phit','k.','markersize',40);
        ylabel('pair hit-rate','fontsize',18);
        xlabel('subject number','fontsize',18);
        mycolors = get(gca,'ColorOrder');
        %mycolors = cat(1,[0 0 0],mycolors);
        myconds = conds;
        %shadedErrorBar(0:size(phit',1)+1,repmat(mean(phit'),[1 size(phit',1)+2]),repmat(ste(phit'),[1 size(phit',1)+2]),{'--','linewidth',2,'color',[0 0 0],'markerfacecolor',[0 0 0]},1);
        for icond = 1:length(myconds)
            shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
        end
        legend({num2str(myconds)},'fontsize',18,'location','best');
        set(gca,'fontsize',18);
        box off
        hold off

        
        %resp bias
        figure('color','white');
        b = bar([pctold' pctnew' pctdknow' pctpdknow']);
        mycolors = nan(4,3);
        cmap = colormap;
        mycolors(1,:) = cmap(1,:);
        mycolors(2,:) = cmap(round(size(cmap,1)/3),:);
        mycolors(3,:) = cmap(round(size(cmap,1)*2/3),:);
        mycolors(4,:) = cmap(end,:);
        hold on;
        shadedErrorBar(0:length(pctold)+1,repmat(mean(pctold),1,length(pctold)+2),repmat(ste(pctold'),1,length(pctold)+2),{'--','linewidth',2,'color',mycolors(1,:),'markerfacecolor',mycolors(1,:)},1);
        shadedErrorBar(0:length(pctnew)+1,repmat(mean(pctnew),1,length(pctnew)+2),repmat(ste(pctnew'),1,length(pctnew)+2),{'--','linewidth',2,'color',mycolors(2,:),'markerfacecolor',mycolors(2,:)},1);
        shadedErrorBar(0:length(pctdknow)+1,repmat(mean(pctdknow),1,length(pctdknow)+2),repmat(ste(pctdknow'),1,length(pctdknow)+2),{'--','linewidth',2,'color',mycolors(3,:),'markerfacecolor',mycolors(3,:)},1);
        shadedErrorBar(0:length(pctpdknow)+1,repmat(mean(pctpdknow),1,length(pctpdknow)+2),repmat(ste(pctpdknow'),1,length(pctpdknow)+2),{'--','linewidth',2,'color',mycolors(4,:),'markerfacecolor',mycolors(4,:)},1);
        ylabel('Response Percentage (phase3)','fontsize',18);
        xlabel('Subject Number','fontsize',18);
        legend({'old','new','Don''t Know','pair DK'},'fontsize',18);
        set(gca,'fontsize',18);
        box off;
        
    end
    
    
    %output
    respbias.pctold = pctold;
    respbias.pctnew = pctnew;
    respbias.pctdknow = pctdknow;
    respbias.pctpdknow = pctpdknow;
    
    results.hit = hit;
    results.miss = miss;
    results.cr = cr;
    results.fa = fa;
    results.dprime = dprime;
    results.pairhit = phit;
    
    bycond.hit = chit;
    bycond.fa = cfa;
    bycond.conds = conds;
    bycond.dprime = cdprime;
    bycond.miss = cmiss;
    bycond.cr = ccr;
    bycond.pairhit = cphit;
    
    results.bycond = bycond;
    
    results.respbias = respbias;
    
    
    out(iblk).data = logdata;
    out(iblk).results = results;
    
end


