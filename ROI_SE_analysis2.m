close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');

flyID = 'fly1_exp1_15Jun23';

% where the 2P data is located
twoPDirectory = ['\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\15Jun23\' flyID];

%where the sequence data is located (stimulus files)
sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data';

%%

nVol = 10;

slices = 5;

% data = readtable(fullfile(twoPDirectory,'measurements'));

green_channel = readtable(fullfile(twoPDirectory,'green_channel'));
red_channel = readtable(fullfile(twoPDirectory,'red_channel'));

% data = readtable('Results_23Feb_fly2_exp2');
% data = readtable('test');

% startVol = 1;

% meanFB = green_channel.Mean_FB_;
% meanRedChannel = red_channel.Mean_Whole_;

% allVolumes = zeros(55,1000);
% 
% for i = 0:1000-1
%     allVolumes(:,i+1) = meanFB((i*slices*nVol)+1:(i+1)*slices*nVol).';
% end
% 
% grandMean = mean(allVolumes,2);
% 
% grandMeanVector = repmat(grandMean.',[1 1000]);
% 
% plot(grandMeanVector);

% meanEB = data.Mean_EBOrFB_(startVol:end);
% meanEBTOP = data.Mean_EB_TOP_(startVol:end);
% meanMBL = data.Mean_MB_L_(startVol:end);
% meanMBR = data.Mean_MB_R_(startVol:end);
% meanOverall = data.Mean1(startVol:end);
% figure; plot(meanFB);
 
% meanFB = normalize(meanFB)-normalize(red_channel.Mean_Whole_);
% hold on; plot(green_channel.Mean_Whole_,'r');
% figure; plot(green_channel.Mean_Whole_);

% meanFB = normalize(meanFB)-normalize(grandMeanVector).';

% figure; plot(meanFB);
% meanFB = reshape(meanFB,[5 length(meanFB)/5]);
% meanFB = mean(meanFB);

degrees = 3;
filt_width = 55;

% plot(green_channel.Mean_FB_(1:220));

% remove first volume of every series for green and red channels
meanFB = reshape(green_channel.Mean_FB_,[5 11000]);
meanFB(:,1:11:end) = [];
meanFB = reshape(meanFB,[1 50000]);

wholeGreen = reshape(green_channel.Mean_Whole_,[5 11000]);
wholeGreen(:,1:11:end) = [];
wholeGreen = reshape(wholeGreen,[1 50000]);

wholeRed = reshape(red_channel.Mean_FB_,[5 11000]);
wholeRed(:,1:11:end) = [];
wholeRed = reshape(wholeRed,[1 50000]);

% try just applying an sgolay filter
filtFB = sgolayfilt(meanFB,degrees,filt_width);
figure; plot(meanFB(1:330)); hold on; plot(filtFB(1:330));
blah = meanFB-sgolayfilt(meanFB,degrees,filt_width);
figure; plot(blah(1:330));
meanFB = meanFB-sgolayfilt(meanFB,degrees,filt_width);

% try subtracting the red channel then apply sgolay
% meanFB = normalize(meanFB)-normalize(wholeRed);
% meanFB = meanFB-sgolayfilt(meanFB,degrees,filt_width);

% try subtracting the whole green channel then apply sgolay
% meanFB = normalize(meanFB)-normalize(wholeGreen);
% meanFB = meanFB-sgolayfilt(meanFB,degrees,filt_width);

meanFB = reshape(meanFB,[5 length(meanFB)/5]);
meanFB = mean(meanFB);

figure; plot(meanFB);

ROIs = [meanFB];
% ROIs = [meanEB.'; meanMBL.'; meanMBR.'];
% ROIs = [meanEB.'; meanMBL.'; meanMBR.'];
% ROIs = meanOverall.';

nROI = size(ROIs,1);

nBack = 5; nSeq = 2^nBack;

ROIs = normalize(ROIs,2);

%% 

% stimuli_filenames = {'SE_fly1_exp4_24Mar2023_1000_trials'};
% 
% % load sequence data
% randomSequence = load(stimuli_filenames{1},'randomSequence');

randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';

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
for roi = 1:nROI
    figure; plot(squeeze(allTraces(roi,:,:)).'); 
end

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

%remove outliers on a sequence basis
for i = 1:16
    this_std = std(roiSeq(:,:,i,:),[],4,'omitnan');
end

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

figure_names = {'FB' 'MB_left' 'MB_right'};

for roi = 1:1
   
    for vol = 1:nVol
    
        figure; create_seq_eff_plot(squeeze(meanROISeq(roi,vol,:)),[]);
    
    end
    
    saveas(gcf,[figure_names{roi} '.png']);
    
end

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