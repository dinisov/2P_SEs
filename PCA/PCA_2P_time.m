close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');
addpath('D:\group_swinderen\Dinis\Extracted datasets\');

resultsDirectory = '../../2P Results';

outputDirectory = '../../2P Clusters';
% outputDirectory = '../../PCA Results';

blocks = readtable('../../2P Record/2P_record');

% blocks_to_exclude = [18 ]

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

flyList = unique(blocks.Fly);

% chosenFlies = 25:35; %cholinergic CC LED
% chosenFlies = [19:24 57:63];% pan-neuronal LED
chosenFlies = [65];

imageSize = [32 32];

trim = 3;

imageSize = imageSize-2*trim;

n_pixels = prod(imageSize);

time_point = 1;

%%

FLIES = struct;

for fly = chosenFlies
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
    
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)],'results.mat'));

        %first time point is usually anomalous
        results.meanDataSeq = results.meanDataSeq(1:end,:,:,:); 
        results.meanBlankTransient = results.meanBlankTransient(:,:,1:end);
%         results.meanTransient = results.meanTransient(:,:,2:end);

        sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
%         sizeAux = size(results.meanTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];

        % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
        results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[3 1 2]),sizeAux),[1 16 1 1]);
%         results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanTransient,[3 1 2]),sizeAux),[1 16 1 1]);

        % construct a matrix of SE profiles for a given time point
        SEProfiles = permute(squeeze(results.meanDataSeq(1,:,:,:)),[2 3 1]);
%         SEProfiles = permute(squeeze(max(results.meanDataSeq,[],1)-min(results.meanDataSeq,[],1)),[2 3 1]); % amplitude

        % construct a matrix of activity over time
        activities = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 1]);
        
        % remove sides from images for clustering
        activities([1:trim end-(trim-1):end],:,:) = [];
        activities(:, [1:trim end-(trim-1):end],:) = [];

%         activities = makeTransient(activities);

        SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
       
        %data matrix for activities
        XAct = reshape(activities,[imageSize(1)*imageSize(2) size(activities,3)]);

        %data matrix for SEs
        XSeq = reshape(SEProfiles,[imageSize(1)*imageSize(2) 16]);
        
        FLIES(fly).BLOCK(b).XAct = XAct;
        FLIES(fly).BLOCK(b).XSeq = XSeq;
        
    end
    
end

%% grouped analysis

X_All = [];

%which blocks to use for each fly
% blocksFlies = {[1],[1]}; %#ok<NBRAK>

% grouped cluster analysis
for fly = 1:length(chosenFlies)

    for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
        X_All = [X_All; (FLIES(chosenFlies(fly)).BLOCK(b).XSeq)]; %#ok<AGROW>
    end

end

%% PCA per fly

% close all

% normalise data

% for i = 1:size(X_All,1)
%    X_All(i,:) = normalize(X_All(i,:));  %#ok<SAGROW>
% end

% load jentzsch_data.mat

% options = statset('MaxIter',10000,'TolFun',1e-8);

for fly = 1:length(FLIES)

    for b = 1:length(FLIES(fly).BLOCK)

        thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)],'PCA');
%         thisFlyDirectory = outputDirectory;
        if ~exist(thisFlyDirectory,'dir')
           mkdir(thisFlyDirectory); 
        end

        [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq.');
        
        for i = 1:5
           figure; imagesc(reshape(coeff(:,i),imageSize)); colorbar; colormap(jet(256));
           saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
           close;     
%            figure; create_seq_eff_plot(-score(:,i),[]);
        end
        
        figure; plot(explained);
        saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
        close;

        [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq);

        for i = 1:5
           figure; create_seq_eff_plot(coeff(:,i),[]);
           saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
           close;     
        end

        figure; plot(explained);
        saveas(gcf,fullfile(thisFlyDirectory,'explained_seq.png')); 
        close;

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

% [L1,T] = rotatefactors(coeff(:,1:2),'Method','procrustes','Target',[SLRP -LRPR],'Maxit',15000,'type','orthogonal');
% [L1,T] = rotatefactors(coeff(:,1:2),'Method','equamax');
% 
% for i = 1:2
%    figure; create_seq_eff_plot(L1(:,i),[]); 
% end

%% PCA for all flies

% [coeff,score,latent,tsquared,explained,mu] = pca((X_All));
% 
% for i = 1:5
%    figure; create_seq_eff_plot(coeff(:,i),[]);
% %    saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '.png']));
% %            close all;     
% end
% 
% figure; plot(explained);