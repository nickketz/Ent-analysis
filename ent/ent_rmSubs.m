function newexper = ent_rmSubs(exper,badsubs)

% remove subjects from exper struct and from subsequent analysis
% input:
%   exper: experiment structure
%   badsubs: binary vector with 1's identifing which subjects to remove 
%
% ouput:
%   newexper: modified exper struct

conds = exper.eventValues;
newexper = exper;
%remove bad subs from nTrials and badEv
tmptrials = exper.nTrials;
for icond = 1:length(conds)
    tmptrials.(conds{icond}) = tmptrials.(conds{icond})(~badsubs);
end
newexper.nTrials = tmptrials;

%remove bad subs from bad chan and bad ev
newexper.badEv = exper.badEv(~badsubs);
newexper.badChan = exper.badChan(~badsubs);

%remove bad subs
newexper.subjects = exper.subjects(~badsubs);

