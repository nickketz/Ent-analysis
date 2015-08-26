function [exper,newdata] = ent_makeAccConds(cfg,adFile)
% make correct and incorrect conditions based on phase 3 accuracy
% todo: pair-correct or just old/new correct?
%
% 

%load analysis details
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);
conds = ana.eventValues{1};
subnum = regexp(exper.subjects,'_([0-9]+)$','tokens');
subnum = cellfun(@(x) (x{1}),subnum);
subnum = cellfun(@(x) (str2double(x)),subnum);
to = ana.trl_order;
for icond = 1:length(conds)
    to.(conds{icond}) = cat(2,'cndn',to.(conds{icond})(1:end-1));
end

behdata = out.data.data;

for isub = 1:length(exper.subjects)
    %icond = 1;
    for icond = 1:length(conds)
        %load subject raw ft data
        load(fullfile(dirs)

        subpow = data_pow.ses1.(conds{icond}).sub(isub).data;
        
        
        %find subject column
        tmp = out.data.data.Subj(1,:);
        icol = tmp == subnum(isub);
        
        %iterate through trials and get acc info
        iphs = behdata.Phase(:,icol) == 3;
        ti = subpow.trialinfo;
        cor = nan(size(ti,1),1);
        for itrial = 1:size(ti,1)
            ordr = ti(itrial,ismember(to.(conds{icond}),'ordr'));
            iordr = behdata.Order(:,icol) == ordr;
            
            blkn = ti(itrial,ismember(to.(conds{icond}),'blkn'));
            iblkn = behdata.Block(:,icol) == blkn;
            
            trialind = iordr & iblkn & iphs;
            if sum(trialind) ~= 1   error('non-unique trialindex found');   end
            cor(itrial) = behdata.Correct(trialind,icol);
        end
        subpow.trialinfo = cat(2,ti,cor);
        data_pow.ses1.(conds{icond}).sub(isub).data.trialinfo = subpow;        
    end
end






    



