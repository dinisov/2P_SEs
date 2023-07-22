function R = calculateR2(meanVolData,ephys)
%calculateR2 Summary of this function goes here
%   Detailed explanation goes here

options = optimset('Algorithm','interior-point','FinDiffType','central');

load('slrp_lrpr.mat','slrp','lrpr','weird');

R = struct;

gridSize = [size(meanVolData,1) size(meanVolData,2)];

r_squared_slrp = zeros(gridSize);
r_squared_lrpr = zeros(gridSize);
r_squared_overall = zeros(gridSize);
r_squared_weird = zeros(gridSize);
r_squared_ephys = zeros(gridSize);

xSize = gridSize(1);
ySize = gridSize(2);

parfor xInd = 1:xSize

    for yInd = 1:ySize

        seq_eff_pattern = squeeze(meanVolData(xInd,yInd,:));

        %fit only to slrp or lrpr
       [~,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
       [~,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
       [~,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
       [~,sse_ephys] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,ephys,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);

       %fit overall
       [~,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

       sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);

       r_squared_slrp(xInd,yInd) = 1-(sse_slrp/sse_total);
       r_squared_lrpr(xInd,yInd) = 1-(sse_lrpr/sse_total);
       r_squared_weird(xInd,yInd) = 1-(sse_weird/sse_total);
       r_squared_overall(xInd,yInd) = 1-(sse_overall/sse_total);
       r_squared_ephys(xInd,yInd) = 1-(sse_ephys/sse_total);

    end

end

R.r_squared_slrp = r_squared_slrp;
R.r_squared_lrpr = r_squared_lrpr;
R.r_squared_weird = r_squared_weird;
R.r_squared_overall = r_squared_overall;
R.r_squared_ephys = r_squared_ephys;

end

