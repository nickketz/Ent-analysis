function entAnalBlkEffects(out,myvar)

nblks = 12;
%myvar = 'hit';


%check for block effects
cfg.filefiltstr = 'EntAssocWordFreq_12blks_1reps_16trials_Sub[0-9]{1,2}\.txt';
cfg.block = 1:nblks;
cfg.doplots = 0;
%out = entAnalFreqBlk(cfg);
nsubs = length(out(1).data.logname);
%outvar = reshape(cell2mat(arrayfun(@(x)(x.results.(myvar)),out,'uniformoutput',false)),nsubs,nblks)';

X = ones(nblks*nsubs,1);
Y = nan(nblks*nsubs,1);
for iblk = cfg.block
    sind = ((iblk-1) * nsubs)+1;
    eind = sind+nsubs-1;
    X(sind:eind,1) = iblk;
    Y(sind:eind,1) = out(iblk).results.(myvar)';
end
fprintf('across conditions test of %s block effects\n---------------',myvar);
fitlm(X,Y)
for isub = 1:nsubs
    lm{isub} = fitlm(X(isub:nsubs:end),Y(isub:nsubs:end));
    beta(isub) = lm{isub}.Coefficients.Estimate(2);
end
fprintf('\n\nheierarchical test of %s block effects\n---------------',myvar);
fitlm(beta,1:nsubs)

% %block effects by condition
% conds = out(iblk).results.bycond.conds;
% nconds = length(conds);
% condvar = nan(nconds,nsubs,length(cfg.block));
% for iblk = cfg.block
%     %conds{iblk} = out(iblk).results.bycond.conds;
%     condvar(:,:,iblk) = cat(1,out(iblk).results.bycond.(myvar){:});
% end
% 
% diff = squeeze(condvar(1,:,:) - condvar(2,:,:));
% Y2 = reshape(diff,nsubs*nblocks,1);
% fprintf('\n\nbetween conditions test of %s block effects\n---------------',myvar);
% fitlm(X,Y2)


