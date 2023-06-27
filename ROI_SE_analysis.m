close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');

% where the 2P data is located
twoPDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\15Jun23\fly2_exp1_15Jun23';

%where the sequence data is located (stimulus setup files)
sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data';

%%

data = readtable(fullfile(twoPDirectory,'measurements'));
% data = readtable('Results_23Feb_fly2_exp2');
% data = readtable('test');

startVol = 1;

meanFB = data.Mean_FB_(startVol:end);
% meanEB = data.Mean_EBOrFB_(startVol:end);
% meanEBTOP = data.Mean_EB_TOP_(startVol:end);
% meanMBL = data.Mean_MB_L_(startVol:end);
% meanMBR = data.Mean_MB_R_(startVol:end);
% meanOverall = data.Mean1(startVol:end);

ROIs = [meanFB.'];
% ROIs = [meanEB.'; meanMBL.'; meanMBR.'];
% ROIs = [meanEB.'; meanMBL.'; meanMBR.'];
% ROIs = meanOverall.';

nROI = size(ROIs,1);

nVol = 11;

nBack = 5; nSeq = 2^nBack;

ROIs = normalize(ROIs,2);

%% 

% stimuli_filenames = {'SE_fly1_exp4_24Mar2023_1000_trials'};
% 
% % load sequence data
% randomSequence = load(stimuli_filenames{1},'randomSequence');

randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs','fly2_exp1_15Jun23.csv')).';

% randomSequence = randomSequence.randomSequence-1;

sequenceLength = length(randomSequence);

%% global response trace

nTrials = 1000;

allTraces = zeros(nROI, nTrials, nVol);

for roi = 1:nROI

    for trial = 0:nTrials-1

        allTraces(roi, trial+1, :) = ROIs(roi, (nVol*trial + 1):(nVol*trial + nVol));

    end

end

% onion plot
% for roi = 1:nROI
%     figure; plot(squeeze(allTraces(roi,:,:)).'); 
% end

meanTrace = squeeze(mean(allTraces,2));

for roi = 1:nROI
    figure; plot(meanTrace(:,roi)); 
end

%plot over time
for roi = 1:nROI
    figure; plot(mean(squeeze(allTraces(roi,:,2:end)).'));
end

%% sort according to sequence

roiSeq = zeros([nROI nVol nSeq sequenceLength/nBack]);

%sort images according to sequence
for n = 0:(sequenceLength/nBack-1)
    
    % decimal value of binary sequence of length n_back
    seq = bin2dec(num2str(randomSequence(n*nBack + 1:n*nBack+(nBack-1)))) + 1;

    for roi = 1:nROI
    
        for vol = 1:nVol

            % this should automatically discard extra volumes at end of
            % session
            roiSeq(roi, vol, seq, n+1) = ROIs(roi,n*nVol+vol);
        
        end
    
    end

end

% join sequences corresponding to =the same pattern (01001 and 10110 and so on)
roiSeq = roiSeq + flip(roiSeq,3);
roiSeq = roiSeq(:,:,1:16,:);
roiSeq = roiSeq(:,:,seq_eff_order(nBack),:);

% in preparation for calculating the nan mean
roiSeq(roiSeq == 0) = nan;

% mean across fourth dimension
meanROISeq = mean(roiSeq, 4, 'omitnan');

stdROISeq = std(roiSeq,[],4,'omitnan');

nERPs = sum(~isnan(roiSeq), 4);

semROISeq = stdROISeq./(sqrt(nERPs));

% reorder according to the literature
% meanROISeq = meanROISeq(:,:,seq_eff_order(nBack));

%% sequential effects profiles per ROI and time step (volume)

% figure_names = {'FB' 'MB_left' 'MB_right'};
% 
% for roi = 1:1
%    
%     for vol = 1:nVol
%     
%         figure; create_seq_eff_plot(-squeeze(meanROISeq(roi,vol,:)),[]);
%     
%     end
%     
%     saveas(gcf,[figure_names{roi} '.png']);
%     
% end

%% trace as a function of sequence per ROI

for roi = 1:nROI
   figure; plot(squeeze(meanROISeq(roi,:,:))); 
end

%% plot AAAA vs AAAR and RRRR-RRRA with error bars

seqs_aux = squeeze(meanROISeq(1,:,:));

std_aux = squeeze(semROISeq(1,:,:));

%AAAA vs AAAR
figure; errorbar(seqs_aux(:,16),std_aux(:,16)); hold on; errorbar(seqs_aux(:,8),std_aux(:,8));

%RRRR vs RRRA
figure; errorbar(seqs_aux(:,1),std_aux(:,1)); hold on; errorbar(seqs_aux(:,9),std_aux(:,9));