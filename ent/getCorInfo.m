%add accuracy into trialinfo

%load analysis details
%adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/test/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/pilot/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
load(adFile)

% %remove ENT_20, there's no behav data right now
% badInd = strcmp(exper.subjects,'ENT_20');
% exper = ent_rmSubs(exper,badInd);

% %remove subcorr temporarily for sub20
% for icond = 1:length(conds)
%     ana.trl_order.(conds{icond}) = {'subn'  'blkn'  'phsn'  'trln'  'ordr'  'stmt'  'pfrq'  'itid'  'corr'};
% end

conds = {'flckr0','flckr6','flckr10','flckr20'};
subnum = regexp(exper.subjects,'_([0-9]+)$','tokens');
subnum = cellfun(@(x) (x{1}),subnum);
subnum = cellfun(@(x) (str2double(x)),subnum);
to = ana.trl_order;

mysubs = {'ENT_16'}

%load behav data
out = ent_behavior(dirs,exper);
behdata = out.data.data;

chngto = 0;
%isub = 1;
for isub = find(ismember(exper.subjects,mysubs))'%14:length(exper.subjects)
    %icond = 1;
    fprintf('%s:\n\t',exper.subjects{isub});
    for icond = 1:length(conds)
        fprintf('%s, ',conds{icond});
        
        if exper.nTrials.(conds{icond})(isub) == 0 
            fprintf('\nNo trials for %s %s\n',exper.subjects{isub},conds{icond});
            continue;
        end
        
        %load data
        fname = fullfile(dirs.saveDirProc, exper.subjects{isub}, exper.sessions{1}{1}, ['data_pow_' conds{icond} '.mat']);
        load(fname);
        
        subpow = freq;
        
        %find subject column
        tmp = out.data.data.Subj(1,:);
        icol = tmp == subnum(isub);
        
        
        ti = subpow.trialinfo;
        %is the subno first?
        if ti(1,2) == subnum(isub) && strcmp(ana.trl_order.(conds{icond}){1},'subn') && max(ti(:,1)<5)
            %remove junk first col which gets added during import
            ti = ti(:,2:end);
            subpow.trialinfo = ti;
        elseif length(ana.trl_order.(conds{icond})) ~= size(ti,2)
            fprintf('\ntrial order doesn''t match trialinfo size, skipping %s %s (probably already done)\n',exper.subjects{isub}, conds{icond});
            continue
        end

        
        %iterate through trials and get acc info
        iphs = behdata.Phase(:,icol) == 3;
        cor = nan(size(ti,1),1);
        for itrial = 1:size(ti,1)
            ordr = ti(itrial,ismember(to.(conds{icond}),'ordr'));
            iordr = behdata.Order(:,icol) == ordr;
            
            blkn = ti(itrial,ismember(to.(conds{icond}),'blkn'));
            iblkn = behdata.Block(:,icol) == blkn;
            
            trialind = iordr & iblkn & iphs;
            if sum(trialind) ~= 1   error('non-unique trialindex found');   end
            %check for 'OD' response
            if strcmp(behdata.Answer{trialind,icol}, 'OD')
                cor(itrial) = 12;
            else
                cor(itrial) = behdata.Correct(trialind,icol);
            end
        end
        
        %subpow.trialinfo = cat(2,ti(:,1:end),cor);
        colind = ismember(ana.trl_order.(conds{icond}),'subcorr');
        if sum(colind)>0
            subpow.trialinfo(:,colind) = cor;            
        else
            subpow.trialinfo(:,end+1) = cor;
            chngto = 1;
        end
        freq = subpow;
        fprintf('\b\b\nsaving...');
        save(fname,'freq','-v7.3');
        
    end
    fprintf('done\n\n');
end


if chngto
    for icond = 1:length(conds)
        ana.trl_order.(conds{icond}) = cat(2,ana.trl_order.(conds{icond}),'subcorr');
%         ana.trl_order.(conds{icond}) = {'subn'  'blkn'  'phsn'  'trln'  'ordr'  'stmt'  'pfrq'  'itid'  'corr'};
    end
end
