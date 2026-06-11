function [chosenFlies,chosenBlocks] = flyProvider(recordPath,commentKeyword,options)
%Function for providing chosenFlies/Blocks for specified criteria, such as keyword presence

%{{
arguments
    %FLIES struct
    %flyRecord table
    recordPath string = "\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx" %Full path to fly record
    %upMode double = 1 %Update mode (1 - Analysis success, 2 - Experiment details [Unfinished])
    commentKeyword = [] %What to search the comment field for; Accepted as either as single char or cell array
    %options.analysisState double = 0; %For upMode 1, what state to update record to (0 - Incomplete/Failed analysis, 1 - Completed analysis)
    options.lastSearchRow double = 2048; %How many rows of Excel to try read (Note this *must* be longer than excel row length)
    %options.stateColumn char = 'AF' %Which column of excel relates to AnalysisSuccess property
    options.commentColumn char = 'AA' %Which column of excel relates to Comments property
    options.expandedArchitecture double = 0 %Whether to arbitrarily replicate chosenFlies/flatten chosenBlocks for certain scripts necessity
    options.omitBattery double = 1 %Whether to omit Battery design from provided blocks
    options.displayTypeColumn char = 'S'  %Which column relates to 'displayType' (i.e. LEDs, Block, etc)
end
%}
%{
upMode = 1; %1 - Pre/Post appendation of analysis success, 2 - New data autofilling
recordPath = "\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx";
analysisState = 0; %Not a parameter
lastSearchRow = 999;
stateColumn = 'AF'; %What column of excel relates to analysis success state
%}

%analysisState = options.analysisState;
lastSearchRow = options.lastSearchRow;
commentColumn = options.commentColumn;
expandedArchitecture = options.expandedArchitecture;
omitBattery = options.omitBattery;
displayTypeColumn = options.displayTypeColumn;

%%

    %Abort early if not actually asked to do anything
    if isempty(commentKeyword) %Note: Will be expanded soon to exp type etc
        disp(['-- No search requested of flyProvider; Exiting --'])
        return
    end

    commentKeywordActual = {};
    if ischar(commentKeyword)
        disp(['Searching for fly/s with comment: ',commentKeyword])
        commentKeywordActual = {commentKeyword}; %Convert to cell, for simplicity
    elseif iscell(commentKeyword)
        disp(['Searching for fly/s with comments: '])
        for cInd = 1:size(commentKeyword,2)
            disp([commentKeyword{cInd}])
        end    
        commentKeywordActual = commentKeyword;
    else
        ['-# Unknown commentKeyword type #-']
        commentKeyword
        crash = yes
    end

    %Load flyRecord 
    flyRecord = readtable(recordPath);

    %Start excel server
    e=actxserver('excel.application');
    
    eW=e.Workbooks;
    %eF = eW.Open("\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx");
    disp([recordPath])
    eF = eW.Open(recordPath);
    eS=eF.ActiveSheet;
    
    %if upMode == 1

    allFliesBlocks = get(eS,'Range',['B2:D',num2str(lastSearchRow)]);
    allFliesBlocks = allFliesBlocks.Value;
    allFliesBlocks = cell2mat(allFliesBlocks);
    eS=eF.ActiveSheet; %May not be necessary to reopen but seems useful

    %Get all fly/blocks Comments
    allComments = get(eS,'Range',[commentColumn,'2:',commentColumn,num2str(lastSearchRow)]);
    allComments = allComments.Value; %Note: Bunch of cells
    eS=eF.ActiveSheet; %May not be necessary to reopen but seems useful

    %Get all fly/blocks design
    allDesign = get(eS,'Range',[displayTypeColumn,'2:',displayTypeColumn,num2str(lastSearchRow)]);
    allDesign = allDesign.Value; %Note: Bunch of cells
    eS=eF.ActiveSheet; %May not be necessary to reopen but seems useful

    %Get all exclusion
    allExclusion = get(eS,'Range',['W2:','W',num2str(lastSearchRow)]);
    allExclusion = allExclusion.Value; %Note: Bunch of cells
    allExclusion = cell2mat(allExclusion);
    eS=eF.ActiveSheet; %May not be necessary to reopen but seems useful


    %Find
    matchMake = nan(size(commentKeywordActual,2),lastSearchRow-1);

    for commI = 1:size(commentKeywordActual,2)
        thisComm = commentKeywordActual{commI};
        for rowInd = 1:lastSearchRow-1
            if contains(char( allComments{rowInd} ),thisComm)
                matchMake(commI,rowInd) = 1;
            else
                matchMake(commI,rowInd) = 0;
            end
        end
        disp([num2str( nansum(matchMake(commI,:)) ),' matches were found for ',thisComm])
    end
    disp(['Match combinations found: ',num2str(unique(nansum(matchMake,1)))]) 
        %0 means non-matches were found (obvs), 1 means either only 1 comment or no flies/blocks matched both, etc

%NOTE: NEED TO ADD EXCLUSION COLUMN CHECKS

    %Report
    for commI = 1:size(commentKeywordActual,2)
        thisComm = commentKeywordActual{commI};
        disp([char(10),thisComm])
        temp = find(matchMake(commI,:) == 1);

        disp(allFliesBlocks(temp,:))
    end

    %Remove battery/exclusion
    if omitBattery
        omissionCount = 0;
        for commI = 1:size(matchMake,1)
            temp = find(matchMake(commI,:) == 1);
            temp2 = allDesign(temp);
            for subI = 1:size(temp,2)
                if contains(temp2{subI},'Battery')
                    matchMake(commI, temp(subI) ) = 0;
                    omissionCount = omissionCount + 1;
                    %disp( [num2str(allFliesBlocks(temp(subI),:)),' detected as battery and omitted'] )
                end
            end
        end  
        if omissionCount ~= 0
            disp([num2str(omissionCount),' blocks were detected as battery and omitted'])
        end
    end
    excludeCount = 0;
    for commI = 1:size(matchMake,1)
        temp = find(matchMake(commI,:) == 1);
        temp2 = find(allExclusion == 1);
        temp3 = intersect(temp,temp2); %NOTE: UNTESTED MATHS
        matchMake( commI, temp3 ) = 0;
        excludeCount = excludeCount + numel(temp3);
    end
    if excludeCount ~= 0
        disp([num2str(excludeCount),' blocks were requested to be excluded'])
    end


    %Assemble chosenFlies/etc
    chosenFlies = cell(1,nanmax(allFliesBlocks(:,1))); %Note slightly different structure, for parallelism with chosenBlocks for the moment
    chosenBlocks = cell(1,nanmax(allFliesBlocks(:,1)));

    for commI = 1:size(commentKeywordActual,2)
        temp = find(matchMake(commI,:) == 1);
        for subI = 1:numel(temp)
            chosenFlies{ allFliesBlocks( temp(subI) , 1 ) } = allFliesBlocks( temp(subI) , 1 ); %Lots of empty cells; Vaguely Dinisian
            chosenBlocks{ allFliesBlocks( temp(subI) , 1 ) } = [chosenBlocks{ allFliesBlocks( temp(subI) , 1 ) }, ...
                allFliesBlocks( temp(subI) , 3 ) ];
        end
    end
    chosenFlies = cell2mat(chosenFlies); %Removes empty/etc
    chosenBlocks = chosenBlocks(chosenFlies); %This operation only valid because chosenFlies self-descriptive
    for fly = 1:size(chosenBlocks,2)
        chosenBlocks{fly} = unique(chosenBlocks{fly}); %Note: Will sort and also merge, obviously
    end
    %QA
    if size(chosenFlies,2) ~= size(chosenBlocks,2)
        ['## Alert: Critical error in chosenFlies/Blocks assembly ##']
        crash = yes
    end
    disp(['Final selected chosenFlies:'])
    disp([chosenFlies])
    disp(['Final selected chosenBlocks:'])
    disp([chosenBlocks])

    %Adjust if necessary
    if expandedArchitecture
        disp(['Expanding architecture per request'])
        chosenFliesNew = [];
        chosenBlocksNew = [];
        for i = 1:size(chosenBlocks,2)
            chosenFliesNew = [chosenFliesNew,...
                repmat( chosenFlies(i) , 1, size(chosenBlocks{i},2) ) ];
            chosenBlocksNew = [chosenBlocksNew, chosenBlocks{i}];
        end
        chosenFlies = chosenFliesNew;
        chosenBlocks = chosenBlocksNew;

        disp(['Expanded chosenFlies/Blocks:'])
        disp([chosenFlies])
        disp([chosenBlocks])

    end

    %Save and close
    %eF.Save; %Not necessary here
    eF.Close;
        
    %end
    
    %jienyo
    %eF.Save;
    %eF.Close;
    e.Quit;
    e.delete;
    %clear e eF eS eW
    disp(['Flies successfully provided'])
end