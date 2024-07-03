function globalTransients(R, chosenFlies, outputDirectory)

for fly = 1:length(R)
    
    for b = [R(fly).BLOCK.blockNum] 
        thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Transients');
        disp(thisFlyDirectory);
        if ~exist(thisFlyDirectory,'dir')
           mkdir(thisFlyDirectory); 
        end
            
        results = R(fly).BLOCK(b);
        
        trim = results.Trim;

        % if blank blocks were use as pedestal, otherwise use mean
        % transient
        if isfield(results,'meanBlankTransient')
            sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
    
            % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[3 1 2]),sizeAux),[1 16 1 1]);
        else
            sizeAux = size(results.meanTransient); %sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
            results.meanTransient = results.meanTransient(:,:,1:end);
        end

        transients = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 1]);
            
        % remove sides from images
        transients([1:trim end-(trim-1):end],:,:) = [];
        transients(:, [1:trim end-(trim-1):end],:) = [];

        % calculate mean and plot
        meanTransient = squeeze(mean(mean(transients,1),2));

%         stdMeanTransient = std(meanTransients,[],2);
% 
%         semMeanTransient = stdMeanTransient./sqrt(totalBlocks);

        plot(meanTransient); xlim([0 length(meanTransient)+1]);

        saveas(gcf,fullfile(thisFlyDirectory,'global_transient.png'));
        close;
    end
        
end

%% group and average transients

% close all


% % this assumes all transient matrices are the same size and were trimmed
% % the same way
% allTransients = zeros([size(FLIES(1).BLOCK(1).transients) totalBlocks]);
% 
% block = 1;
% 
% for fly = 1:length(chosenFlies)
% 
%     for b = 1:length(FLIES(fly).BLOCK)
%         
%         allTransients(:,:,:,:,block) = FLIES(fly).BLOCK(b).transients;
% 
%         block = block + 1;
% 
%     end
% 
% end
% 
% % calculate mean and plot
% meanTransients = squeeze(mean(mean(mean(allTransients,1),2),3));
% 
% meanTransient = mean(meanTransients,2);%singular
% 
% stdMeanTransient = std(meanTransients,[],2);
% 
% semMeanTransient = stdMeanTransient./sqrt(totalBlocks);
% 
% errorbar(meanTransient,semMeanTransient); xlim([0 size(meanTransients,1)+1]);
% 
% save('before_transient','meanTransient','semMeanTransient');
% % save('during_transient','meanTransient','semMeanTransient');
% % save('after_transient','meanTransient','semMeanTransient');