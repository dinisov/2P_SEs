close all; clear;

addpath('..\Functions\');

fly = 25;
b = 1;
trim = 3;

resultsDirectory = '../../2P Results';
thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(fly)],['Block' num2str(b)]);
pca_results = load(fullfile(thisFlyDirectory,'PCA','pca_results_normalised'),'coeff','score');

brainImage = imread(fullfile(thisFlyDirectory,'brain.jpg'));

% imshow(brainImage);

pca_component = reshape(pca_results.score(:,2),[32 32]-2*trim);

% pca_component = (pca_component - min(pca_component,[],'all'));
% pca_component = pca_component./max(pca_component,[],'all');

% figure; imagesc(pca_component); colorbar; colormap(jet(256));

trimmed_brain_img = brainImage(2*trim*16+1:end-(2*trim*16),2*trim*16+1:end-(2*trim*16));

plotBrainPCA(pca_component,trimmed_brain_img,.1,'on');