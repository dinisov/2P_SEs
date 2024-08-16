function flyPCA2D(R, chosenFlies, outputDirectory)

load('six_hertz.mat','six_hertz');

% choose from ['time','sequence']
pcaType = {'sequence'};

% number of components to retain for time and sequence
n_comp_seq = 5;
n_comp_t = 3;

%%

FLIES = struct;

for fly = 1:length(R)
    
    for b = [R(fly).BLOCK.blockNum]    
        
        results = R(fly).BLOCK(b);
        
        trim = results.Trim;
        
        imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
    
        % if blank blocks were collected use as pedestal, otherwise use mean
        % transient
        if isfield(results,'meanBlankTransient')
            sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
    
            % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[3 1 2]),sizeAux),[1 16 1 1]);
        else
            sizeAux = size(results.meanTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
            results.meanTransient = results.meanTransient(:,:,1:end);
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanTransient,[3 1 2]),sizeAux),[1 16 1 1]);
        end

        if any(strcmp(pcaType,'sequence'))
            % construct a matrix of SE profiles averaged across time
            SEProfiles = permute(squeeze(sum(results.meanDataSeq,1)),[2 3 1]);
    
            % trim sides
            SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
    
            %data matrix for SEs
            XSeq = reshape(SEProfiles,[imageSize(1)*imageSize(2) 16]);
            
            FLIES(fly).BLOCK(b).XSeq = XSeq.';
            FLIES(fly).BLOCK(b).numInstances = length( find( squeeze( sum( results.dataSeq, [1:4]) ) ~= 0 ) ); %Find non-zero layers and imply size of actual data
        end
        
%         if any(strcmp(pcaType,'time'))
%             % construct a matrix of activity over time
%             activities = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 1]);
%             
%             % trim sides
%             activities([1:trim end-(trim-1):end],:,:) = [];
%             activities(:, [1:trim end-(trim-1):end],:) = [];
% 
%             %data matrix for activities
%             XAct = reshape(activities,[imageSize(1)*imageSize(2) size(activities,3)]);
% 
%             FLIES(fly).BLOCK(b).XAct = XAct.';
%         end

        %And behav separated data, if applicable
        if isfield( R(fly).BLOCK(b) , 'dataSeqBehav' )
            for statInd = 1:size( R(fly).BLOCK(b).dataSeqBehav,2 )
                imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
                if isfield(results,'meanBlankTransient')
                    sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
                    % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
                    results.dataSeqBehav(statInd).meanDataSeqReduced = results.dataSeqBehav(statInd).meanDataSeqReduced./repmat(reshape(permute(results.meanBlankTransient,[3 1 2]),sizeAux),[1 16 1 1]);
                else
                    sizeAux = size(results.meanTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
                    results.meanTransient = results.meanTransient(:,:,1:end);
                    results.dataSeqBehav(statInd).meanDataSeqReduced = results.dataSeqBehav(statInd).meanDataSeqReduced./repmat(reshape(permute(results.meanTransient,[3 1 2]),sizeAux),[1 16 1 1]);
                end
                if any(strcmp(pcaType,'sequence'))
                    % construct a matrix of SE profiles averaged across time
                    SEProfiles = permute(squeeze(sum(results.dataSeqBehav(statInd).meanDataSeqReduced,1)),[2 3 1]);
                    % trim sides
                    SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
                    %data matrix for SEs
                    XSeq = reshape(SEProfiles,[imageSize(1)*imageSize(2) 16]);
                    %FLIES(fly).BLOCK(b).XSeq = XSeq.';
                    FLIES(fly).BLOCK(b).dataSeqBehav(statInd).XSeq = XSeq.';
                    FLIES(fly).BLOCK(b).dataSeqBehav(statInd).state = results.dataSeqBehav(statInd).state;
                    FLIES(fly).BLOCK(b).dataSeqBehav(statInd).numInstances = results.dataSeqBehav(statInd).numInstances; %As above, but using provided information
                end
            end            
        end
        
    end
    
end

%% PCA per fly

for fly = 1:length(R)

    for b = [R(fly).BLOCK.blockNum]
        
        trim = R(fly).BLOCK(b).Trim;
        
        %imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;

        brainImage = R(fly).BLOCK(b).brainImage;
        
        trimmedBrainImg = brainImage(2*trim*16+1:end-(2*trim*16),2*trim*16+1:end-(2*trim*16));
        
        data = struct;
        
        %------
        c = 1;
        data(c).name = 'All';
        data(c).XSeq = FLIES(fly).BLOCK(b).XSeq;
        data(c).imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim; %Moved here to be in setup loop
        data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'PCA');
        %disp(data(c).thisFlyDirectory);
        if ~exist(data(c).thisFlyDirectory,'dir')
           mkdir(data(c).thisFlyDirectory); 
        end
        data(c).numInstances = FLIES(fly).BLOCK(b).numInstances; %Such data not available (currently) for All data
        c = c + 1;
        %The reason for this rigamarole is so behav separated data can be added if existing
        if isfield( R(fly).BLOCK(b) , 'dataSeqBehav' )
            for statInd = 1:size( R(fly).BLOCK(b).dataSeqBehav,2 )
                data(c).name = ['State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state)];
                data(c).XSeq = FLIES(fly).BLOCK(b).dataSeqBehav(statInd).XSeq;
                data(c).imageSize = size(R(fly).BLOCK(b).dataSeqBehav(statInd).meanDataSeqReduced,[3 4]) - 2 * trim;
                data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],['State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state)],'PCA');
                %disp(data(c).thisFlyDirectory);
                if ~exist(data(c).thisFlyDirectory,'dir')
                   mkdir(data(c).thisFlyDirectory); 
                end
                if isfield(FLIES(fly).BLOCK(b).dataSeqBehav(statInd), 'numInstances')
                    data(c).numInstances = FLIES(fly).BLOCK(b).dataSeqBehav(statInd).numInstances;
                else
                    data(c).numInstances = [];
                end
                c = c + 1;
            end
        end
        %------
        
        
        for datInd = 1:size(data,2) %Makes the code run for as many loops as there are data (1 - n)
            thisData = data(datInd).XSeq;
            thisFlyDirectory = data(datInd).thisFlyDirectory;
            imageSize = data(datInd).imageSize;

            if any(strcmp(pcaType,'sequence'))
                tic

                %[coeff,score,~,~,explained,~] = pca(FLIES(fly).BLOCK(b).XSeq);
                [coeff,score,~,~,explained,~] = pca( thisData );
                
                numInstances = data(datInd).numInstances;

                for i = 1:n_comp_seq
                   figure; imagesc(reshape(coeff(:,i),imageSize)); colorbar; colormap(jet(256));
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   close;

                   sign_ephys = sortOrientation(score(:,i),normalize(six_hertz));

                   figure; create_seq_eff_plot(normalize(score(:,i)),normalize(sign_ephys*six_hertz));

                   saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   close;

                   % overlay plot on brain
                   plotBrainPCA(reshape(-coeff(:,i),imageSize),trimmedBrainImg,'on');
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '_overlay.png']));
                   close;
                end

                figure; plot(explained);
                saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
                close;

                %save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained');
                save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained','numInstances');
                disp(['Processed PCA for ',data(datInd).name,' in ',num2str(toc),'s'])
                disp(thisFlyDirectory);

            end

            %The upstream processing for this has been disabled currently
            %{
            if any(strcmp(pcaType,'time'))

                %thisFlyDirectory = fullfile(resultsDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'PCAact');
                thisFlyDirectory = strrep( thisFlyDirectory, 'PCA', 'PCAact' );
                if ~exist(thisFlyDirectory,'dir')
                   mkdir(thisFlyDirectory); 
                end

                %[coeff,score,~,~,explained,~] = pca(FLIES(fly).BLOCK(b).XAct);
                [coeff,score,~,~,explained,~] = pca( <not implemented> );

                for i = 1:n_comp_t
                   figure; imagesc(reshape(coeff(:,i),imageSize)); colorbar; colormap(jet(256));
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   close;
                   figure; plot(normalize(score(:,i)));
                   saveas(gcf,fullfile(thisFlyDirectory,['c_act' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   close;

                   % overlay plot on brain
                   plotBrainPCA(reshape(coeff(:,i),imageSize),trimmedBrainImg,'on');
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '_overlay.png']));
                   close;
                end

                figure; plot(explained);
                saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
                close;

                save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained');

            end
            %}
        
        
        end
        

    end

end

function sign_ephys = sortOrientation(scores,ephys)

    sse_plus = sum((normalize(scores)-normalize(ephys)).^2);
    sse_minus = sum((normalize(scores)-normalize(-ephys)).^2);

    if sse_plus < sse_minus
        sign_ephys = 1;
    else
        sign_ephys = -1;
    end

end

end