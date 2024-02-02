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
chosenFlies = [19:24 57:63];% pan-neuronal LED
% chosenFlies = [19:24];

imageSize = [32 32];

trim = 3;

imageSize = imageSize-2*trim;

n_pixels = prod(imageSize);

%%

FLIES = struct;

for fly = chosenFlies
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
    
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)],'results.mat'));

        results.meanDataSeq = results.meanDataSeq(2:end,:,:,:); %first time point is usually anomalous
        
        % construct a matrix of SE profiles averaged across time
        SEProfiles = permute(squeeze(sum(results.meanDataSeq,1)),[2 3 1]);

        % construct a matrix of activity iver time
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
        
        % make pairwise distance matrix
%         distMat = squareform(pdist(X,'correlation'));
%         
%         imagesc(distMat);

        %% k means or k meoids

%         n_clusters = 3;

%         idx = spectralcluster(XSeq,n_clusters);
%         idx = kmeans(XSeq,n_clusters,'Distance','correlation');
%         idx = kmedoids(XSeq,n_clusters,'Distance','correlation');
        
%         for cl = 1:n_clusters
            
%             cluster_mask = reshape(idx == cl,imageSize);
            
%             figure; imagesc(cluster_mask); saveas(gcf,['cluster_mask' num2str(cl) '.png']);
            
%             figure; plot(mean(XAct(cluster_mask(:),:))); saveas(gcf,['cluster_activity' num2str(cl) '.png']);
            
%             se_profile = zeros(16,1);
%             %SE profile for cluster
%             for s = 1:16
%                se_profile(s) = mean(XSeq(cluster_mask(:),s));
%             end
            
%             figure; create_seq_eff_plot(se_profile,[]); saveas(gcf,['SE_profile_cluster' num2str(cl) '.png']);

            % plot profiles over time
%             for t = 1:size(activities,3)
%     
%                se_profile = zeros(16,1);
%                %SE profile for cluster
%                for s = 1:16
%                    thisSeq = squeeze(results.meanDataSeq(t,s,:,:));
%                    thisSeq = thisSeq(:);
%                    se_profile(s) = mean(thisSeq(cluster_mask(:)));
%                end
%                figure; create_seq_eff_plot(se_profile,[]); %saveas(gcf,['SE_profile_cluster' num2str(i) '.png']);
%     
%             end

%         end
        
        % hierarchical clustering (probably not relevant)
%         T = clusterdata(X,5);
        
    end
    
end

%% grouped clustering analysis

X_All = [];

%which blocks to use for each fly
% blocksFlies = {[1],[1]}; %#ok<NBRAK>

% grouped cluster analysis
for fly = 1:length(chosenFlies)

    for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
        X_All = [X_All; (FLIES(chosenFlies(fly)).BLOCK(b).XSeq)]; %#ok<AGROW>
    end

end

% n_clusters = 3;

% idx = kmeans(X_All,n_clusters,'Distance','correlation');
% idx = kmedoids(X_All,n_clusters,'Distance','correlation');
% [idx,V,D] = spectralcluster(X_All,n_clusters,'Distance','correlation','ClusterMethod','kmedoids');

% clusterMasks = cell(1,n_clusters*n_flies);

% se_profiles = zeros(16,n_clusters);

%% plot stuff

% for cl = 1:n_clusters
%     
%    subDirectory = fullfile(outputDirectory,['Cluster' num2str(cl)]);
%    if ~exist(subDirectory,'dir')
%       mkdir(subDirectory); 
%    end
%    
%    block_count = 0;
%    
%     for fly = 1:length(chosenFlies)
%         for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
%             cluster_mask = reshape(idx((block_count*n_pixels + 1):n_pixels*(block_count+1)) == cl,imageSize);
%             figure('Visible','off'); 
% %                 figure; 
%             imagesc(cluster_mask); saveas(gcf,fullfile(subDirectory,['fly' num2str(chosenFlies(fly)) '_block' num2str(b) '_cluster' num2str(cl) '.png']));
%             close all;
%             block_count = block_count+1;
%         end
%     end
% 
%    %SE profile for cluster
%    for s = 1:16
%        se_profiles(s,cl) = mean(X_All(idx == cl,s));
%    end
% 
% %    figure;
%    figure('Visible','off'); 
%    create_seq_eff_plot(se_profiles(:,cl),[]); saveas(gcf,fullfile(subDirectory,['SE_profile_cluster' num2str(cl) '.png']));
%    close all;
% end

%% PCA

% close all

% normalise data

% for i = 1:size(X_All,1)
%    X_All(i,:) = normalize(X_All(i,:));  %#ok<SAGROW>
% end

% load jentzsch_data.mat

% options = statset('MaxIter',10000,'TolFun',1e-8);

% for fly = 1:length(FLIES)
% 
%     for b = 1:length(FLIES(fly).BLOCK)
% 
%         thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)],'PCA');
% %         thisFlyDirectory = outputDirectory;
%         if ~exist(thisFlyDirectory,'dir')
%            mkdir(thisFlyDirectory); 
%         end
% 
%         [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq.');
%         
%         for i = 1:10
%            figure; imagesc(reshape(coeff(:,i),imageSize)); colorbar; colormap(jet(256));
%            figure; create_seq_eff_plot(-score(:,i),[]);
% %            saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
% %            close all;     
%         end
%         
%         figure; plot(explained);
% %         saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); close all;
% 
% %         [coeff,score,latent,tsquared,explained,mu] = pca(FLIES(fly).BLOCK(b).XSeq);
% 
% %         for i = 1:5
% %            figure; create_seq_eff_plot(coeff(:,i),[]);
% % %            saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(fly) '_' num2str(b) '.png']));
% % %            close all;     
% %         end
% 
% %         figure; plot(explained);
% %         saveas(gcf,fullfile(thisFlyDirectory,'explained_seq.png')); %close all;
% 
%     end
% 
% end

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

[coeff,score,latent,tsquared,explained,mu] = pca((X_All));

for i = 1:10
   figure; create_seq_eff_plot(coeff(:,i),[]);
%    saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '.png']));
%            close all;     
end

figure; plot(explained);

%% factor analysis

% n_fac = 2;
% 
% % options = statset('MaxIter',15000);
% 
% [lambda,psi] = factoran(X_All,n_fac,'maxit',1500000);
% 
% for i = 1:n_fac
%    figure; create_seq_eff_plot(lambda(:,i),[]); 
% end

        %% DBSCAN

% minpts = 20;
% 
% kD = pdist2(X_All,X_All,'correlation','Smallest',minpts);
% 
% plot(sort(kD(end,:)));
% title('k-distance graph')
% xlabel('Points sorted with 50th nearest distances')
% ylabel('Nth nearest distances')
% grid
% 
% epsilon = 0.12;
% 
% idx = dbscan(X_All,epsilon,minpts,'Distance','correlation');
% 
% for cl = 1:max(unique(idx))
%     
%    subDirectory = fullfile(outputDirectory,['Cluster' num2str(cl)]);
%    if ~exist(subDirectory,'dir')
%       mkdir(subDirectory); 
%    end
% 
%     block_count = 0;
%     
%     for fly = 1:length(chosenFlies)
%         for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
%             cluster_mask = reshape(idx((block_count*n_pixels + 1):n_pixels*(block_count+1)) == cl,imageSize);
%             figure('Visible','off'); 
% %             figure; 
%             imagesc(cluster_mask); saveas(gcf,fullfile(subDirectory,['fly' num2str(chosenFlies(fly)) '_block' num2str(b) '_cluster' num2str(cl) '.png']));
%             close all;
%             block_count = block_count+1;
%         end
%     end
% 
%    se_profile = zeros(16,1);
%    %SE profile for cluster
%    for s = 1:16
%        se_profile(s) = mean(X_All(idx == cl,s));
%    end
%    %    figure;
%    figure('Visible','off'); 
%    create_seq_eff_plot(se_profile,[]); saveas(gcf,fullfile(subDirectory,['SE_profile_cluster' num2str(cl) '.png']));
%    close all;
% 
% end

function data = makeTransient(data)

%     data = data(:,:,2:end); %get rid of first time point
    data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1
    
end