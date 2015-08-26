tmp = out.results.bycond.pairhit;
%pairhit = reshape(cell2mat(tmp),[size(tmp{1},2) size(tmp,2)]);
pairhit = cell2mat(tmp)';
tmp = out.results.bycond.dprime;
dprime = cell2mat(tmp)';

[h,p,ci,stats] = ttest(dprime,repmat(dprime(:,1),[1 size(dprime,2)]))
[h,p,ci,stats] = ttest(pairhit,repmat(pairhit(:,1),[1 size(pairhit,2)]))

bline = 1;

%dprime by condition, works with multiple conditions
figure('color','white');
hold on
tmp = dprime-repmat(dprime(:,bline),1,size(dprime,2));
tmp = tmp(:,1:size(tmp,2)~=bline);
plot(tmp,'.','markersize',30);
ylabel('baseline corrected dprime','fontsize',18);
xlabel('subject number','fontsize',18);
mycolors = get(gca,'ColorOrder');
conds = out.results.bycond.conds;
myconds = conds(1:length(conds)~=bline);
for icond = 1:length(myconds)
    shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
end
conds = out.results.bycond.conds;
legend({num2str(myconds)},'fontsize',18,'location','best');
set(gca,'fontsize',18);
box off
hold off

%pairhit by condition, works with multiple conditions
figure('color','white');
hold on
tmp = pairhit;
plot(tmp,'.','markersize',30);
ylabel('pair hit-rate','fontsize',18);
xlabel('subject number','fontsize',18);
mycolors = get(gca,'ColorOrder');
myconds = out.results.bycond.conds;
for icond = 1:length(myconds)
    shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
end
legend({num2str(myconds)},'fontsize',18,'location','best');
set(gca,'fontsize',18);
box off
hold off

%pairhit by condition, works with multiple conditions
figure('color','white');
hold on
tmp = pairhit-repmat(pairhit(:,bline),1,size(pairhit,2));
tmp = tmp(:,1:size(tmp,2)~=bline);
plot(tmp,'.','markersize',30);
ylabel('baseline corrected pair hit-rate','fontsize',18);
xlabel('subject number','fontsize',18);
mycolors = get(gca,'ColorOrder');
conds = out.results.bycond.conds;
myconds = conds(1:length(conds)~=bline);
for icond = 1:length(myconds)
    shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
end
legend({num2str(myconds)},'fontsize',18,'location','best');
set(gca,'fontsize',18);
box off
hold off