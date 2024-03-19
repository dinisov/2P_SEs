% make transient figure

close all; clear;

before = load('before_transient.mat');
during = load('during_transient.mat');
after = load('after_transient.mat');

errorbar([before.meanTransient during.meanTransient after.meanTransient],[before.semMeanTransient during.semMeanTransient after.semMeanTransient],'linewidth',1.5);

set(gca,'yticklabels',[],'xtick',0:12); 

xlim([0.5 11.5])

legend({'Baseline','Red Light ON','Recovery'},'box','off','fontsize',12);

xlabel('Time post-stimulus','fontsize',12); ylabel('dF/F','fontsize',12);