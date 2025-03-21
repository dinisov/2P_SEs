function flyPCA2D(R, chosenFlies, outputDirectory, options)
arguments
    R struct
    chosenFlies double
    outputDirectory
    options.showModel double = 1
    options.pcaType string = {'sequence'}
end

if options.showModel
load('six_hertz.mat','six_hertz');
end

% choose from ['time','sequence']
%pcaType = {'sequence'};
pcaType = options.pcaType;

% number of components to retain for time and sequence
n_comp_seq = 5;
n_comp_t = 3;

%%

FLIES = struct;

for fly = 1:length(R)
    
    %for b = [R(fly).BLOCK.blockNum]    
    for b = 1:size( R(fly).BLOCK,2 )
        blockNum = R(fly).BLOCK(b).blockNum;
        
        results = R(fly).BLOCK(b);
        
        %n_
        %trim = results.Trim;
        if ~isfield(results,'TrimCoords') || isempty(results.TrimCoords)
            trim = results.Trim;
        else
            trim = results.TrimCoords;
        end
        
        %imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
        asymmTrim = 0;
        if isequal( size(trim) , [1,1] )
            imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
        elseif isequal( size(trim) , [1,4] )
            imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - [trim(1)+trim(3),trim(2)+trim(4)]; %"Top, Bottom, Right, Left
            asymmTrim = 1;
        else
            ['## Invalid trim parameters specified! ##']
            disp(['Fly #',num2str(chosenFlies(fly)),' Block ',num2str(results.blockNum)])
            crash = yes
        end
    
        % if blank blocks were collected use as pedestal, otherwise use mean
        % transient
        if isfield(results,'meanBlankTransient') && ~isempty( results.meanBlankTransient )
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
            if asymmTrim == 0
                SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
            else
                SEProfiles = SEProfiles(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
            end
    
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
                %imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
                if asymmTrim == 0
                    imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
                else
                    imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - [trim(1)+trim(3),trim(2)+trim(4)]; %"Top, Bottom, Right, Left
                end
                if isfield(results,'meanBlankTransient') && ~isempty( results.meanBlankTransient )
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
                    %SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
                    if asymmTrim == 0
                        SEProfiles = SEProfiles(trim+1:end-trim,trim+1:end-trim,:);
                    else
                        SEProfiles = SEProfiles(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
                    end
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

    %for b = [R(fly).BLOCK.blockNum]
    for b = 1:size( R(fly).BLOCK,2 )
        blockNum = R(fly).BLOCK(b).blockNum;
        
        %d_
        %trim = R(fly).BLOCK(b).Trim;
        if ~isfield(R(fly).BLOCK(b),'TrimCoords') || isempty(R(fly).BLOCK(b).TrimCoords)
            trim = R(fly).BLOCK(b).Trim;
        else
            trim = R(fly).BLOCK(b).TrimCoords;
        end
        
        %imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
        imageSize = []; %Clear just in case
        asymmTrim = 0;
        if isequal( size(trim) , [1,1] )
            imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim;
        elseif isequal( size(trim) , [1,4] )
            imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - [trim(1)+trim(3),trim(2)+trim(4)]; %"Top, Bottom, Right, Left
            asymmTrim = 1;
        else
            ['## Invalid trim parameters specified! ##']
            disp(['Fly #',num2str(chosenFlies(fly)),' Block ',num2str(results.blockNum)])
            crash = yes
        end

        brainImage = R(fly).BLOCK(b).brainImage;
        
        %%trimmedBrainImg = brainImage(2*trim*16+1:end-(2*trim*16),2*trim*16+1:end-(2*trim*16));
        if asymmTrim == 0
            trimmedBrainImg = brainImage(2*trim*16+1:end-(2*trim*16),2*trim*16+1:end-(2*trim*16));
        else
            trimmedBrainImg = brainImage(trim(1)*16+1:end-(trim(3)*16),trim(4)*16+1:end-(trim(2)*16)); %Adjusted to be four-coordinate
        end
        
        data = struct;
        
        %------
        c = 1;
        data(c).name = 'All';
        data(c).XSeq = FLIES(fly).BLOCK(b).XSeq;
        %data(c).imageSize = size(R(fly).BLOCK(b).meanDataSeq,[3 4]) - 2 * trim; %Moved here to be in setup loop
        data(c).imageSize = imageSize; %Slight adjustment in case of asymmetrical trim
        %data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'PCA');
        data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(blockNum)],'PCA');
        %disp(data(c).thisFlyDirectory);
        if ~exist(data(c).thisFlyDirectory,'dir')
           mkdir(data(c).thisFlyDirectory); 
        end
        %And delete old directory, if existing (with care)
        if exist( data(c).thisFlyDirectory ) == 7
            %QA
            temp = dir( data(c).thisFlyDirectory );
            if nansum( [temp.isdir] ) > 2 %This is the expected number of 'directories' for a directory that contains no subdirectories
                    %"." and ".." are special directory markers, not true directories
                ['-# Cannot safely delete directory due to detected presence of subdirectories #-']
            else
                %[data(c).thisFlyDirectory,filesep,'*.png']
                %rmdir( data(c).thisFlyDirectory ,'s' ) %Unsafe
                %delete( [data(c).thisFlyDirectory,filesep,'*.png'], [data(c).thisFlyDirectory,filesep,'*.mat'] )
                delete( [data(c).thisFlyDirectory,filesep,'*.png'], [data(c).thisFlyDirectory,filesep,'*.mat'], [data(c).thisFlyDirectory,filesep,'*.fig']  )
                rmdir( data(c).thisFlyDirectory ) %Will only succeed if empty
                    %If not empty, this may have acted on a folder it wasn't meant to
                mkdir(data(c).thisFlyDirectory); 
                disp(['(Existing directory contents deleted, directory remade)'])
            end            
        end
        
        data(c).numInstances = FLIES(fly).BLOCK(b).numInstances; %Such data not available (currently) for All data
        c = c + 1;
        %The reason for this rigamarole is so behav separated data can be added if existing
        if isfield( R(fly).BLOCK(b) , 'dataSeqBehav' )
            for statInd = 1:size( R(fly).BLOCK(b).dataSeqBehav,2 )
                imageSize = [];
                if asymmTrim == 0
                    imageSize = size(R(fly).BLOCK(b).dataSeqBehav(statInd).meanDataSeqReduced,[3 4]) - 2 * trim;
                else
                    imageSize = size(R(fly).BLOCK(b).dataSeqBehav(statInd).meanDataSeqReduced,[3 4]) - [trim(1)+trim(3),trim(2)+trim(4)];
                end
                data(c).name = ['State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state)];
                data(c).XSeq = FLIES(fly).BLOCK(b).dataSeqBehav(statInd).XSeq;
                %data(c).imageSize = size(R(fly).BLOCK(b).dataSeqBehav(statInd).meanDataSeqReduced,[3 4]) - 2 * trim;
                data(c).imageSize = imageSize;
                %data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],['State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state)],'PCA');
                data(c).thisFlyDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(blockNum)],['State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state)],'PCA');
                %disp(data(c).thisFlyDirectory);
                if ~exist(data(c).thisFlyDirectory,'dir')
                   mkdir(data(c).thisFlyDirectory); 
                end
                %And delete old directory, if existing (with care)
                    %Copy of above
                if exist( data(c).thisFlyDirectory ) == 7
                    %QA
                    temp = dir( data(c).thisFlyDirectory );
                    if nansum( [temp.isdir] ) > 2 %This is the expected number of 'directories' for a directory that contains no subdirectories
                            %"." and ".." are special directory markers, not true directories
                        ['-# Cannot safely delete directory due to detected presence of subdirectories #-']
                    else
                        %rmdir( data(c).thisFlyDirectory ,'s' ) %Unsafe
                        delete( [data(c).thisFlyDirectory,filesep,'*.png'], [data(c).thisFlyDirectory,filesep,'*.mat'] )
                        rmdir( data(c).thisFlyDirectory ) %Will only succeed if empty
                            %If not empty, this may have acted on a folder it wasn't meant to
                        mkdir(data(c).thisFlyDirectory); 
                        disp(['(Existing directory contents deleted, directory remade)'])
                    end            
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
                   %saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(blockNum) '.png']));
                   close;

                   if options.showModel
                   sign_ephys = sortOrientation(score(:,i),normalize(six_hertz));

                   %figure; create_seq_eff_plot(normalize(score(:,i)),normalize(sign_ephys*six_hertz));
                   figure; create_seq_eff_plot_Matt(normalize(score(:,i)),normalize(sign_ephys*six_hertz));
                   else
                   %figure; create_seq_eff_plot(normalize(score(:,i)), []); %Hide model
                   figure; create_seq_eff_plot_Matt(normalize(score(:,i)), []); %Hide model
                    %Matt version exists so that changes can be made outside of Global Functions folder
                   end

                   %saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '.png']));
                   saveas(gcf,fullfile(thisFlyDirectory,['c_seq' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(blockNum) '.png']));
                   close;

                   % overlay plot on brain
                   plotBrainPCA(reshape(-coeff(:,i),imageSize),trimmedBrainImg,'on');
                   %saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(b) '_overlay.png']));
                   saveas(gcf,fullfile(thisFlyDirectory,['c' num2str(i) '_fly_' num2str(chosenFlies(fly)) '_' num2str(blockNum) '_overlay.png']));
                   close;
                end

                figure; plot(explained);
                saveas(gcf,fullfile(thisFlyDirectory,'explained.png')); 
                close;

                %save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained');
                save(fullfile(thisFlyDirectory,'pca_results_normalised'),'coeff','score','explained','numInstances');
                disp(['Processed PCA for ',data(datInd).name,' in ',num2str(toc),'s'])
                disp(thisFlyDirectory);

                %Make/save plot of components in space
                if n_comp_seq >= 3
                    colMap = jet(16);
                    %Purely components
                    figure
                    %scatter3( score(:,1), score(:,2) , score(:,3) )
                    scatter3( score(:,1), score(:,2) , score(:,3), [], colMap, 'filled')
                    xlabel(['PC1 (',num2str( explained(1) ),'%)'])
                    ylabel(['PC2 (',num2str( explained(2) ),'%)'])
                    zlabel(['PC3 (',num2str( explained(3) ),'%)'])
                    title(['PC1 x PC2 x PC3'])
                    saveas(gcf,fullfile(thisFlyDirectory,'PC1x2x3.png')); 
                    saveas(gcf,fullfile(thisFlyDirectory,'PC1x2x3_figure'),'fig'); 

                    %First two PCs and seq.
                    figure
                    %scatter3( score(:,1), score(:,2) , score(:,3) )
                    scatter3( score(:,1), score(:,2) , [1:size(score,1)], [], colMap, 'filled')
                    zlim([0,size(score,1)+1])
                    zticks([1:size(score,1)])
                    xlabel(['PC1 (',num2str( explained(1) ),'%)'])
                    ylabel(['PC2 (',num2str( explained(2) ),'%)'])
                    zlabel(['Seq.'])
                    title(['PC1 x PC2 x Seq.'])
                    saveas(gcf,fullfile(thisFlyDirectory,'PC1x2xSeq.png')); 
                    saveas(gcf,fullfile(thisFlyDirectory,'PC1x2xSeq_figure'),'fig'); 

                end

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