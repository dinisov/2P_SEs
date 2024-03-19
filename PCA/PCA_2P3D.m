close all; clear;
% 
addpath('..\..\Scripts\Global functions\');
addpath('..\..\Scripts\Indexes and legends\');
addpath('..\..\2P SEs\Functions\');
addpath('..\..\2P SEs\Functions\3D');
addpath('../../Extracted datasets/');

resultsDirectory = '../../2P Results 3D';
dataDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded blocks
blocks = blocks(~logical(blocks.Exclude),:);

% chosenFlies = 25:35; %cholinergic CC LED
% chosenFlies = [19:24 57:63];% pan-neuronal LED
chosenFlies = [95:102 104:106];

imageSize = [32 32];

trim = 3;

imageSize = imageSize-2*trim;

n_pixels = prod(imageSize);

load six_hertz.mat

% choose from ['time','sequence']
pcaType = {'sequence','time'};

% number of components to retain for time and sequence
n_comp_seq = 5;
n_comp_t = 3;

%%

FLIES = struct;

for fly = chosenFlies
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
    
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)],'results.mat'));

        n_z = thisFlyBlocks(b,:).Steps;

        %remove first time point (often it is anomalous)
%         results.meanDataSeq = results.meanDataSeq(1:end,:,:,:);

        % if blank blocks were use as pedestal, otherwise use mean
        % transient
        if thisFlyBlocks(b,:).BlankBlocks
%             results.meanBlankTransient = results.meanBlankTransient(:,:,:);
            sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([4 1 2 3]); sizeAux = [sizeAux(1) 1 sizeAux(2:4)];
    
            % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[4 1 2 3]),sizeAux),[1 16 1 1 1]);
        else
            sizeAux = size(results.meanTransient); sizeAux = sizeAux([4 1 2 3]); sizeAux = [sizeAux(1) 1 sizeAux(2:4)];
%             results.meanTransient = results.meanTransient(:,:,1:end);
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanTransient,[4 1 2 3]),sizeAux),[1 16 1 1 1]);
        end

        if any(strcmp(pcaType,'sequence'))
            % construct a matrix of SE profiles averaged across time
            SEProfiles = permute(squeeze(sum(results.meanDataSeq,1)),[2 3 4 1]);
    
            SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:,:);
    
            %data matrix for SEs
            XSeq = reshape(SEProfiles,[imageSize(1)*imageSize(2)*n_z 16]);
            
    %         FLIES(fly).BLOCK(b).XSeq = XSeq - repmat(mean(XSeq,2),[1 16]);
            FLIES(fly).BLOCK(b).XSeq = XSeq.';

        end
        
        if any(strcmp(pcaType,'time'))
            % construct a matrix of activity over time
            activities = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 4 1]);
            
            % remove sides from images
            activities([1:trim end-(trim-1):end],:,:,:) = [];
            activities(:, [1:trim end-(trim-1):end],:,:) = [];

            %data matrix for activities
            XAct = reshape(activities,[imageSize(1)*imageSize(2)*n_z size(activities,4)]);

            FLIES(fly).BLOCK(b).XAct = XAct.';
        end
        
    end
    
end

%% grouped analysis

% X_All = [];
% 
% %which blocks to use for each fly
% % blocksFlies = {[1],[1]}; %#ok<NBRAK>
% 
% % grouped cluster analysis
% for fly = 1:length(chosenFlies)
% 
%     for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
%         X_All = [X_All; (FLIES(chosenFlies(fly)).BLOCK(b).XSeq)]; %#ok<AGROW>
%     end
% 
% end

%% PCA per fly

% close all

% load jentzsch_data.mat

for fly = chosenFlies
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);

    for b = 1:length(FLIES(fly).BLOCK)
        
        thisBlock = thisFlyBlocks(b,:);

        n_z = thisBlock.Steps; flyback = thisBlock.FlybackFrames;

        currentDate = char(datetime(thisBlock.Date,'Format','dMMMyy'));
        currentFlyDirectory = ['fly' num2str(thisBlock.FlyOnDay) '_exp' num2str(thisBlock.Block) '_' currentDate];

        % might save a volume to make this faster
        load(fullfile(dataDirectory,currentDate,currentFlyDirectory,'green_channel_aligned.mat'),'green_channel_aligned');
        
        % get rid of flyback frames
        green_channel_aligned = green_channel_aligned(:,:,1:(size(green_channel_aligned,3) - flyback),:);
        
        brain = mean(green_channel_aligned,4);
        
        trimmedBrain = brain(2*trim*4+1:end-(2*trim*4),2*trim*4+1:end-(2*trim*4),:);% this assumes a 128x128 image
        
        % normalize the brain for greater contrast (after trimming so no edges)
        trimmedBrain = trimmedBrain - min(trimmedBrain,[],'all');
        trimmedBrain = trimmedBrain ./ max(trimmedBrain,[],'all');
        
        if any(strcmp(pcaType,'sequence'))

            thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)],'PCA');

            if ~exist(thisFlyDirectory,'dir')
               mkdir(thisFlyDirectory); 
            end
        
            [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq);
            
            for i = 1:n_comp_seq
               plot3D(reshape(coeff(:,i),[imageSize n_z]),'on');
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;

               sign_ephys = sortOrientation(score(:,i),normalize(six_hertz));

%                figure; create_seq_eff_plot(normalize(score(:,i)),normalize(sign_ephys*six_hertz));
               figure; create_seq_eff_plot(normalize(score(:,i)),[]); ylabel('PCA score');

               saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;
               
               % overlay plot on brain
               plotBrainPCA3D(imresize3(reshape(coeff(:,i),[imageSize n_z]),size(trimmedBrain),'box'),trimmedBrain,'on');
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '_overlay.png']));
               close;
            end
            
            figure; plot(explained);
            saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
            close;
            
            save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained');

        end

        if any(strcmp(pcaType,'time'))

            thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)],'PCAact');
            if ~exist(thisFlyDirectory,'dir')
               mkdir(thisFlyDirectory); 
            end
            
            [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XAct);

            for i = 1:n_comp_t
               plot3D(reshape(coeff(:,i),[imageSize n_z]),'off'); colorbar; colormap(jet(256));
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;

               plot(normalize(score(:,i)));
               saveas(gcf,fullfile(thisFlyDirectory,['c_act' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;
               
               % overlay plot on brain
               plotBrainPCA3D(imresize3(reshape(coeff(:,i),[imageSize n_z]),size(trimmedBrain),'box'),trimmedBrain,'on');
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '_overlay.png']));
               close;
            end
            
            figure; plot(explained);
            saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
            close;
            
            save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained');

        end

    end

end

% for cp = 1:3
%     se_profile = zeros(16,1);
%     end
%     figure; create_seq_eff_plot(se_profile,[]);
% end

% [coeff,score,latent,tsquared,explained,mu] = pca(XSeq);
% 
% for i = 1:4
%    figure; create_seq_eff_plot(coeff(:,i),[]); 
% end
% 
% figure; plot(explained);

% [L1,T] = rotatefactors(coeff(:,1:2),'Method','procrustes','Target',[SLRP -LRPR],'Maxit',15000,'type','oblique');
% % [L1,T] = rotatefactors(coeff(:,1:2),'Method','equamax');
%  
% for i = 1:2
%    figure; create_seq_eff_plot(L1(:,i),[]); 
% end

%% PCA for all flies

% load six_hertz.mat
% 
% [coeff,score,latent,tsquared,explained,mu] = pca(X_All);
% 
% % for i = 1:5
% %    figure; create_seq_eff_plot(normalize(coeff(:,i)),normalize(six_hertz));
% % %    saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '.png']));
% % %            close all;     
% % end
%  
% figure; plot(explained);
% 
% for c = 1:1
%    
%    for fly = 1:size(X_All,1)/prod(imageSize)
%         figure; imagesc(reshape(score((fly-1)*prod(imageSize)+1:fly*prod(imageSize),c),imageSize)); colorbar; colormap(jet(256));
%    end
% %    saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
% %    close;     
% %    figure; create_seq_eff_plot(score(:,i),[]);
%     
% end

% [L1,T] = rotatefactors(coeff(:,1:2),'Method','pattern','Target',[six_hertz ones(16,1)]);
% % [L1,T] = rotatefactors(coeff(:,1:2),'Method','equamax');
%   
% for i = 1:2
%    figure; create_seq_eff_plot(normalize(L1(:,i)),normalize(six_hertz)); 
% end

%% average some components of interest

close all;

% % the component numbers for each fly/block
% flies = [80 94 98 101 103];
% % normalised by mean transient
% chosen_blocks = {[1],[1],[1],[1],[1]};
% cp_num = {[1],[1],[1],[1],[1]};

% the component numbers for each fly/block
% flies = [80 94 98 101 103 84 83 95 100];
% % normalised by mean transient
% chosen_blocks = {[1],[1],[1],[1],[1],[1],[1],[1],[1]};
% cp_num = {[1],[1],[1],[1],[1],[1],[2],[1],[2]};

% % the component numbers for each fly/block
% flies = [106 107 115 116 118 120];
% % normalised by mean transient
% chosen_blocks = {[1],[1],[1],[1],[1],[1]};
% cp_num = {[1],[1],[1],[1],[1],[1]};

% % the component numbers for each fly/block
% flies = [80,94,98,101,103,106,107,114,115,116,118, 120];
% % normalised by mean transient
% chosen_blocks = {[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};
% cp_num = {[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]};

% minor components 7 stacks
flies = [107 114 116 118 120];
% normalised by mean transient
chosen_blocks = {[1],[1],[1],[1],[1],[1]};
cp_num = {[-3],[-3],[4],[3],[-4],[4]};

% flies = [95:102 104:106];
% component = 5;
% block = 3;
% % % normalised by mean transient
% chosen_blocks = {[block],[block],[block],[block],[block],[block],[block],[block],[block],[block],[block]};
% % chosen_blocks = {[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2]};
% % chosen_blocks = {[3],[3],[3],[3],[3],[3],[3],[3],[3],[3],[3]};
% cp_num = {[component],[component],[component],[component],[component],[component],[component],[component],[component],[component],[component]};

% flies = [80 94];
% % % normalised by mean transient
% chosen_blocks = {[3],[3]};
% cp_num = {[1],[1]};

% flies = [27];
% % % normalised by mean transient
% blocks = {[1 2 3]};
% cp_num = {[3 2 -2]};

% flies = [32];
% % % normalised by mean transient
% blocks = {[1 2]};
% cp_num = {[1 1]};

trim = [3 3 3 3 3 3 3 3 3 3 3 3];

profiles = zeros(16,length(cell2mat(cp_num)));

p = 1;

for fly = 1:length(flies)
    disp(fly);
    imageSize = [32 32]-2*trim(fly);
   
    for b = 1:length(chosen_blocks{fly})
        
        thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(flies(fly))],['Block' num2str(chosen_blocks{fly}(b))],'PCA');
        pca_results = load(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score');
        
        profiles(:,p) = sign(cp_num{fly}(b))*pca_results.score(:,abs(cp_num{fly}(b)));
        
%         plot3D(sign(cp_num{fly}(b))*reshape(pca_results.coeff(:,abs(cp_num{fly}(b))),[imageSize 4]),'on');
%         saveas(gcf,['fly_' num2str(flies(fly)) '_block_' num2str(chosen_blocks{fly}(b)) '_cp_' num2str(abs(cp_num{fly}(b))) '.png']);
%         figure; create_seq_eff_plot(normalize(sign(cp_num{fly}(b))*pca_results.score(:,abs(cp_num{fly}(b)))),[]);
%         saveas(gcf,['fly_' num2str(flies(fly)) '_block_' num2str(chosen_blocks{fly}(b)) '_cp_' num2str(abs(cp_num{fly}(b))) '_profile.png']);

        p = p + 1;
    end
    
end

figure; create_seq_eff_plot(mean(profiles,2),[],'error',std(profiles,[],2)/sqrt(size(profiles,2))); ylabel('Mean PCA Score');

names = {'before','during','after'};

saveas(gcf,['Component averages/mean_cp' num2str(component) '_' names{block} '.png']);
% saveas(gcf,'mean_cp1_during.png');
% saveas(gcf,'mean_cp1_after.png');

%%

% load six_hertz.mat
% 
% load behaviour_results.mat

% figure; create_seq_eff_plot(normalize([mean(profiles,2)]),[]); ylabel('Z-score');
% % legend({'Two-photon','Ephys'},'box','off');
% saveas(gcf,'average_cp1_before.png');
% 
% figure; create_seq_eff_plot(normalize([mean(profiles,2) -RTs]),[]); ylabel('Z-score');
% legend({'Two-photon','Behaviour'},'box','off');
% saveas(gcf,'minor_cps_behav.png');

% figure; create_seq_eff_plot(normalize([mean(profiles,2) six_hertz -RTs]),[]); ylabel('Z-score');
% legend({'Two-photon','Ephys','Behaviour'},'box','off');
% saveas(gcf,'minor_cps_behav_ephys.png');

% figure; create_seq_eff_plot(profiles,[]);
% legend({'Baseline','Red light ON','Recovery'},'box','off'); ylabel('PCA Score');

% figure; create_seq_eff_plot(profiles(:,[1 3]),[]);
% legend({'Baseline','Recovery'},'box','off'); ylabel('PCA Score');

% figure; create_seq_eff_plot(profiles(:,[2]),[]);
% legend({'Red light ON'},'box','off'); ylabel('PCA Score');

function sign_ephys = sortOrientation(scores,ephys)

    sse_plus = sum((normalize(scores)-normalize(ephys)).^2);
    sse_minus = sum((normalize(scores)-normalize(-ephys)).^2);

    if sse_plus < sse_minus
        sign_ephys = 1;
    else
        sign_ephys = -1;
    end

end