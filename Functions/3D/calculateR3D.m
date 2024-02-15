function R = calculateR3D(meanVolData,ephys)
%calculateR Summary of this function goes here
%   Detailed explanation goes here

    R = struct;

    load('slrp_lrpr.mat','slrp','lrpr','weird');
    
    gridSize = [size(meanVolData,1) size(meanVolData,2) size(meanVolData,3)];
    
    r_slrp = zeros(gridSize);
    r_lrpr = zeros(gridSize);
    r_weird = zeros(gridSize);
    r_ephys = zeros(gridSize);

    p_slrp = zeros(gridSize);
    p_lrpr = zeros(gridSize);
    p_weird = zeros(gridSize);
    p_ephys = zeros(gridSize);
    
    xSize = gridSize(1);
    ySize = gridSize(2);
    zSize = gridSize(3);

    parfor xInd = 1:xSize

        for yInd = 1:ySize
            
            for zInd = 1:zSize

               seq_eff_pattern = squeeze(meanVolData(xInd,yInd,zInd,:));

               [r,p] = corrcoef(seq_eff_pattern,slrp);
               r_slrp(xInd,yInd,zInd) = r(2); p_slrp(xInd,yInd,zInd) = p(2);

               [r,p] = corrcoef(seq_eff_pattern,-lrpr);
               r_lrpr(xInd,yInd,zInd) = r(2); p_lrpr(xInd,yInd,zInd) = p(2);

               [r,p] = corrcoef(seq_eff_pattern,weird);
               r_weird(xInd,yInd,zInd) = r(2); p_weird(xInd,yInd,zInd) = p(2);

               [r,p] = corrcoef(seq_eff_pattern,ephys);
               r_ephys(xInd,yInd,zInd) = r(2); p_ephys(xInd,yInd,zInd) = p(2);
           
            end

        end

    end
    
    R.r.r_slrp = r_slrp;
    R.r.r_lrpr = r_lrpr;
    R.r.r_weird = r_weird;
    R.r.r_ephys = r_ephys;
    
    R.p.p_slrp = p_slrp;
    R.p.p_lrpr = p_lrpr;
    R.p.p_weird = p_weird;
    R.p.p_ephys = p_ephys;

end

