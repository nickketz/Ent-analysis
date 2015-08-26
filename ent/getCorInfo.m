%add accuracy into trialinfo

%load analysis details
adFile = '/Users/nketz/Documents/Documents/boulder/Entrain_EEG/Analysis/data/ENT/EEG/Sessions/test/ft_data/flckr0_flckr6_flckr10_flckr20_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
load(adFile)

conds = {'flckr0','flckr6','flckr10','flckr20'};
subnum = regexp(exper.subjects,'_([0-9]+)$','tokens');
subnum = cellfun(@(x) (x{1}),subnum);
subnum = cellfun(@(x) (str2double(x)),subnum);
to = ana.trl_order;


%load behav data
out = ent_behavior(dirs,exper);
behdata = out.data.data;

%isub = 1;
for isub = 1:length(exper.subjects)
    %icond = 1;
    fprintf('%s:\n\t',exper.subjects{isub});
    for icond = 1:length(conds)
        fprintf('%s, ',conds{icond});
                
        %load data
        fname = fullfile(dirs.saveDirProc, exper.subjects{isub}, exper.sessions{1}{1}, ['data_pow_' conds{icond} '.mat']);
        load(fname);
        
        subpow = freq;
                
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
            %check for 'OD' response
            if strcmp(behdata.Answer{trialind,icol}, 'OD')
                cor(itrial) = 12;
            else
                cor(itrial) = behdata.Correct(trialind,icol);
            end
        end
        
        %subpow.trialinfo = cat(2,ti(:,1:end),cor);
        subpow.trialinfo(:,ismember(ana.trl_order.(conds{icond}),'corr')) = cor;
        freq = subpow;  
        fprintf('\b\b\nsaving...');
        save(fname,'freq','-v7.3');
        
    end
    fprintf('done\n\n');
    
end