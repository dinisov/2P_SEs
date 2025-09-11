function recordUpdater(FLIES,flyRecord,recordPath,upMode,options)
%Function for updating 2P record automatically

%{{
arguments
    FLIES struct
    flyRecord table
    recordPath string = "\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx" %Full path to fly record
    upMode double = 1 %Update mode (1 - Analysis success, 2 - Experiment details [Unfinished])
    options.analysisState double = 0; %For upMode 1, what state to update record to (0 - Incomplete/Failed analysis, 1 - Completed analysis)
    options.lastSearchRow double = 1024; %How many rows of Excel to try read (Note this *must* be longer than excel row length)
    options.stateColumn char = 'AF' %Which column of excel relates to AnalysisSuccess property
end
%}
%{
upMode = 1; %1 - Pre/Post appendation of analysis success, 2 - New data autofilling
recordPath = "\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx";
analysisState = 0; %Not a parameter
lastSearchRow = 999;
stateColumn = 'AF'; %What column of excel relates to analysis success state
%}

analysisState = options.analysisState;
lastSearchRow = options.lastSearchRow;
stateColumn = options.stateColumn;

%%

    disp(['Updating records (Mode ',num2str(upMode),')'])

    %Start excel server
    e=actxserver('excel.application');
    
    eW=e.Workbooks;
    %eF = eW.Open("\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx");
    eF = eW.Open(recordPath);
    eS=eF.ActiveSheet;
    
    if upMode == 1
    
        %allFlies = get(eS,'Range','B2:B999');
        %allFlies = get(eS,'Range',['B2:B',num2str(lastSearchRow)]);
        %allFlies = allFlies.Value;
        %allFlies = cell2mat(allFlies); %Will fail on weird entry
        %allBlocks = get(eS,'Range',['C2:C',num2str(lastSearchRow)]);
        %allBlocks = allBlocks.Value;
        %allBlocks = cell2mat(allBlocks); %Will fail on weird entry
        allFliesBlocks = get(eS,'Range',['B2:D',num2str(lastSearchRow)]);
        allFliesBlocks = allFliesBlocks.Value;
        allFliesBlocks = cell2mat(allFliesBlocks);
        eS=eF.ActiveSheet; %May not be necessary to reopen but seems useful
    
        for thisFly = 1:size(FLIES,2)
            for thisBlock = 1:size(FLIES(thisFly).BLOCKS,2)
                %Collect fly/block details from FLIES structure
                    %Better to use this than chosenFlies cos of exclusion/disregardation/etc
                flyNum = FLIES(thisFly).BLOCKS(thisBlock).flyNum;
                blockNum = FLIES(thisFly).BLOCKS(thisBlock).blockNum;
                disp(['Appending analysis success state data for fly #',num2str(flyNum),' block ',num2str(blockNum)])
    
                %Find fly row/s in excel         
                flyNumRows = find( allFliesBlocks(:,1) == flyNum ); %Which rows of excel contain this fly
    
                %QA fly not found
                if isempty(flyNumRows)
                    ['## Alert: Fly not found in specified range of excel ##']
                    e.Quit; e.delete;
                    crash = yes
                end
    
                %Find block in row/s
                blockNumRows = find( allFliesBlocks(flyNumRows,3) == blockNum );
    
                %QA row overfind
                if size(blockNumRows,1) > 1
                    ['## Alert: Block number overfind ##']
                    e.Quit; e.delete;
                    crash = yes
                end
    
                flyBlockRowActual = flyNumRows(blockNumRows)+1;
                disp(['Associated excel row: ',num2str(flyBlockRowActual),' (',num2str(flyBlockRowActual+1),' inc. header)'])
    
                try
                    %Manual testatory
                    %{
                    stateRange = get(eS,'Range',['AF',num2str(flyBlockRowActual)]);
                    disp(['Previous value: ',num2str(stateRange.Value)])
                    stateRange = get(eS,'Range',['AF909']);
                    stateRange.Value = 96;
                    disp(['New value: ',num2str(stateRange.Value)])
                    %}
        
                    %Write state into appropriate position
                    %stateRange = get(eS,'Range',['AF',num2str(flyBlockRowActual)]);
                    stateRange = get(eS,'Range',strcat([stateColumn,num2str(flyBlockRowActual)]));
                    disp(['Previous value: ',num2str(stateRange.Value)])
                    stateRange.Value = analysisState;
                    disp(['New value: ',num2str(stateRange.Value)])
                catch
                    ['-# Failure in saving Excel data #-']
                    e.Quit;
                    e.delete;
                    return
                end
    
                
            end
        end
    
        %Save and close
        eF.Save;
        eF.Close;
        
    elseif upMode == 2
        %flyRecord = readtable("\\uq.edu.au\uq-inst-gateway1\PHDMVDP002-Q1471\TEMP STORE\2P_record.xlsx");
        numRows = size(flyRecord,1);
        flyrecord wrong (not whole size)
        
        disp(['Fly record loading indicates ',num2str(numRows),' rows'])
        disp(['Last fly currently in list: #',num2str(flyRecord(numRows,:).Fly),' from ',datestr(flyRecord(numRows,:).Date)])
        
        falling
        
        
        
        templateRange = get(eS,'Range','A2:AF2');
        blirg = [86,92,33]
        eActivesheetRange.Value = blirg;
        
        eActivesheetRange = get(eS,'Range','J964:L964');
    end
    
    %jienyo
    %eF.Save;
    %eF.Close;
    e.Quit;
    e.delete;
    %clear e eF eS eW
    disp(['Records successfully updated'])
end