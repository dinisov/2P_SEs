close all; clear;

nProfiles = 1000;

randomProfiles = zeros(nProfiles,16);

r_squared_slrp = zeros(1,nProfiles);
r_squared_lrpr = zeros(1,nProfiles);
r_squared_overall = zeros(1,nProfiles);
r_squared_weird = zeros(1,nProfiles);

load slrp_lrpr.mat

options = optimset('Algorithm','interior-point','FinDiffType','central');

for i = 1:nProfiles
    
    randomProfiles(i,:) = rand(1,16);
   
    sse_total = sum((randomProfiles(i,:)-mean(randomProfiles(i,:))).^2);
    
    seq_eff_pattern = randomProfiles(i,:).';

    %fit only to slrp or lrpr
   [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
   [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
   [x_weird,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);

   %fit overall
   [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

   r_squared_slrp(i) = 1-(sse_slrp/sse_total);
   r_squared_lrpr(i) = 1-(sse_lrpr/sse_total);
   r_squared_weird(i) = 1-(sse_weird/sse_total);
   r_squared_overall(i) = 1-(sse_overall/sse_total);

   
end

%% plot the best ones

best_slrp = randomProfiles(r_squared_slrp == max(r_squared_slrp),:);
best_lrpr = randomProfiles(r_squared_lrpr == max(r_squared_lrpr),:);
best_weird = randomProfiles(r_squared_weird == max(r_squared_weird),:);
best_overall = randomProfiles(r_squared_overall == max(r_squared_overall),:);


figure; create_seq_eff_plot(normalize(best_slrp).',normalize(slrp));
figure; create_seq_eff_plot(normalize(best_lrpr).',normalize(lrpr))
figure; create_seq_eff_plot(normalize(best_weird).',-normalize(weird))
% figure; create_seq_eff_plot(normalize(best_overall).',-normalize(slrp))