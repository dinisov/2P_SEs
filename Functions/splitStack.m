function BLOCK = splitStack(BLOCK)
%splitStack Split image stack into relevant and blank trial frames
%   Detailed explanation goes here

    % index of blank stimuli/blocks
    indBlank = BLOCK.randomSequence == 5;
    
    % remove blank stimuli from random sequence
    BLOCK.randomSequence = BLOCK.randomSequence(~indBlank);
    %Apply same indices to behavSequence, if applicable
    if isfield(BLOCK,'behavSequence')
        BLOCK.behavSequence = BLOCK.behavSequence(~indBlank); %Will crash if not exact same length
    end
    
    % construct index of frames corresponding to blank blocks of stimuli
    indBlankStack = reshape(indBlank,[BLOCK.nStimuli length(indBlank)/BLOCK.nStimuli]).';
    indBlankStack = repmat(indBlankStack(:,1),[1 BLOCK.nVol]).';
    
    % split between stimuli image and blank image stacks
    BLOCK.blankImageStack = BLOCK.greenChannel(:,:,indBlankStack(:));
    BLOCK.greenChannel = BLOCK.greenChannel(:,:,~indBlankStack(:));

end