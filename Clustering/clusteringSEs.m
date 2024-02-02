close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');
addpath('D:\group_swinderen\Dinis\Extracted datasets\');

resultsDirectory = '../../2P Results';

outputDirectory = '../../2P Clusters';

blocks = readtable('../../2P Record/2P_record');

% blocks_to_exclude = [18 ]

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

flyList = unique(blocks.Fly);

% chosenFlies = 6:24;%57C10 x GCamp7s
% chosenFlies = 25:35; %cholinergic CC LED
% chosenFlies = [2:5 25:35];
% chosenFlies = 19:24;%57C10 x GCamp7s LED
% chosenFlies = [19:21 23];
% chosenFlies = [19 21 23]; %three best flies
% chosenFlies = 19:23; % usable 57C10 flies 
% chosenFlies = 6:18; %57C10x GCamp7s Projector
chosenFlies = [21 23];

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
    
    for b = 1:1%nBlocks
    
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)],'results.mat'));
        
        % construct a matrix of SE profiles averaged across time
        SEProfiles = permute(squeeze(sum(results.meanDataSeq,1)),[2 3 1]);
        
%         remove sides from images for clustering
        SEProfiles([1:trim end-(trim-1):end],:,:) = [];
        SEProfiles(:, [1:trim end-(trim-1):end],:) = [];
       
        %data matrix
        X = reshape(SEProfiles,[imageSize(1)*imageSize(2) 16]);
        
        FLIES(fly).BLOCK(b).X = normalize(X,2);
        
        % make pairwise distance matrix
%         distMat = squareform(pdist(X,'correlation'));
%         
%         imagesc(distMat);

        %% k means or k meoids

        n_clusters = 2;

        idx = kmeans(X,n_clusters,'Distance','correlation');
%         idx = kmedoids(X,n_clusters,'Distance','correlation');
        
%         for i = 1:n_clusters
%             
%            cluster_mask = reshape(idx == i,imageSize);
%            
% %            cluster_mask([1:3 end-3:end],:) = 0;
% %            cluster_mask(:, [1:3 end-3:end]) = 0;
%            
% %            figure; imagesc(cluster_mask); saveas(gcf,['cluster_mask' num2str(i) '.png']);
%            
%            se_profile = zeros(16,1);
%            %SE profile for cluster
%            for s = 1:16
%                se_profile(s) = mean(X(cluster_mask(:),s));
%            end
%            figure; create_seq_eff_plot(se_profile,[]); saveas(gcf,['SE_profile_cluster' num2str(i) '.png']);
%         end
        
        % hierarchical clustering (probably not relevant)
%         T = clusterdata(X,5);
        
    end
    
end

%% grouped clustering analysis

% n_flies = 18;

X_All = [];

%which blocks to use for each fly
% blocksFlies = {[1],[1]}; %#ok<NBRAK>

% grouped cluster analysis
for fly = 1:length(chosenFlies)

    for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
        X_All = [X_All; FLIES(chosenFlies(fly)).BLOCK(b).X]; %#ok<AGROW>
    end

end

n_clusters = 2;

% idx = kmeans(X_All,n_clusters,'Distance','correlation');
idx = kmeans(X_All,n_clusters,'Distance','correlation');
% [idx,V,D] = spectralcluster(X_All,n_clusters,'Distance','correlation','ClusterMethod','kmedoids');

% clusterMasks = cell(1,n_clusters*n_flies);

se_profiles = zeros(16,n_clusters);

%% plot stuff

for cl = 1:n_clusters
    
   subDirectory = fullfile(outputDirectory,['Cluster' num2str(cl)]);
   if ~exist(subDirectory,'dir')
      mkdir(subDirectory); 
   end
   
   block_count = 0;
   
    for fly = 1:length(chosenFlies)
        for b = 1:length(FLIES(chosenFlies(fly)).BLOCK)
            cluster_mask = reshape(idx((block_count*n_pixels + 1):n_pixels*(block_count+1)) == cl,imageSize);
            figure('Visible','off'); 
%                 figure; 
            imagesc(cluster_mask); saveas(gcf,fullfile(subDirectory,['fly' num2str(chosenFlies(fly)) '_block' num2str(b) '_cluster' num2str(cl) '.png']));
            close all;
            block_count = block_count+1;
        end
    end

   %SE profile for cluster
   for s = 1:16
       se_profiles(s,cl) = mean(X_All(idx == cl,s));
   end

%    figure;
   figure('Visible','off'); 
   create_seq_eff_plot(se_profiles(:,cl),[]); saveas(gcf,fullfile(subDirectory,['SE_profile_cluster' num2str(cl) '.png']));
   close all;
end

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