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
chosenFlies = [107:120];

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
            FLIES(fly).BLOCK(b).XSeq = XSeq;

        end
        
        if any(strcmp(pcaType,'time'))
            % construct a matrix of activity over time
            activities = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 4 1]);
            
            % remove sides from images for clustering
            activities([1:trim end-(trim-1):end],:,:,:) = [];
            activities(:, [1:trim end-(trim-1):end],:,:) = [];

%         activities = makeTransient(activities);

            %data matrix for activities
            XAct = reshape(activities,[imageSize(1)*imageSize(2)*n_z size(activities,4)]);

            FLIES(fly).BLOCK(b).XAct = XAct;
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

for fly = 1:length(FLIES)
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);

    for b = 1:length(FLIES(fly).BLOCK)
        
        thisBlock = thisFlyBlocks(b,:);

        n_z = thisBlock.Steps;

        currentDate = char(datetime(thisBlock.Date,'Format','dMMMyy'));
        currentFlyDirectory = ['fly' num2str(thisBlock.FlyOnDay) '_exp' num2str(thisBlock.Block) '_' currentDate];
        brainImage = imread(fullfile(dataDirectory,currentDate,currentFlyDirectory,'brain.jpg'));
        
        trimmedBrainImg = brainImage(2*trim*16+1:end-(2*trim*16),2*trim*16+1:end-(2*trim*16));

        if any(strcmp(pcaType,'sequence'))

            thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)],'PCA');

            if ~exist(thisFlyDirectory,'dir')
               mkdir(thisFlyDirectory); 
            end
        
            [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq);
            
            for i = 1:n_comp_seq
               plot3D(reshape(score(:,i),[imageSize n_z]),'off'); colorbar; colormap(jet(256));
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;

               figure; create_seq_eff_plot(normalize(coeff(:,i)),normalize(six_hertz));

               saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;
               
               % overlay plot on brain
%                plotBrainPCA(reshape(score(:,i),imageSize),trimmedBrainImg,'on');
%                saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '_overlay.png']));
%                close;
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
               plot3D(reshape(score(:,i),[imageSize n_z]),'off'); colorbar; colormap(jet(256));
               saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;

               plot(normalize(coeff(:,i)));
               saveas(gcf,fullfile(thisFlyDirectory,['c_act' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
               close;
               
               % overlay plot on brain
%                plotBrainPCA(reshape(score(:,i),imageSize),trimmedBrainImg,'on');
%                saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '_overlay.png']));
%                close;
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

% "L/R" component

% the component numbers for each fly/block
flies = [25 27 32 35];
% normalised by mean transient
blocks = {[1 2 3],[1 2 3],[1 2],[1 2 3]};
cp_num = {[-2 -4 -1],[3 2 -2],[1 1],[1 1 2]};

% not normalised
% blocks = {[1],[1 2 3],[1],[2 3],[1 2],[1 2 3]};
% cp_num = {[-3],[-4 2 -3],[-3],[-2 2],[2 -2],[2 -2 2]};

% flies = [25];
% % % normalised by mean transient
% blocks = {[1 2 3]};
% cp_num = {[-2 -4 -1]};

% flies = [35];
% % % normalised by mean transient
% blocks = {[1 2 3]};
% cp_num = {[1 1 2]};

% flies = [27];
% % % normalised by mean transient
% blocks = {[1 2 3]};
% cp_num = {[3 2 -2]};

% flies = [32];
% % % normalised by mean transient
% blocks = {[1 2]};
% cp_num = {[1 1]};

% trim = [3 3 3 3 3 3];
% 
% profiles = zeros(16,length(cell2mat(cp_num)));
% 
% p = 1;
% 
% for fly = 1:length(flies)
%     disp(fly)
%     imageSize = [32 32]-2*trim(fly);
%    
%     for b = 1:length(blocks{fly})
%         
%         thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(flies(fly))],['Block' num2str(blocks{fly}(b))],'PCA');
%         pca_results = load(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score');
%         
%         profiles(:,p) = sign(cp_num{fly}(b))*pca_results.coeff(:,abs(cp_num{fly}(b)));
%         
%         figure; imagesc(reshape(sign(cp_num{fly}(b))*pca_results.score(:,abs(cp_num{fly}(b))),imageSize)); colorbar; colormap(jet(256));
%         figure; create_seq_eff_plot(normalize(sign(cp_num{fly}(b))*pca_results.coeff(:,abs(cp_num{fly}(b)))),normalize(six_hertz));
%         
%         p = p + 1;
%     end
%     
% end
% 
% figure; create_seq_eff_plot(normalize(mean(profiles,2)),[]);