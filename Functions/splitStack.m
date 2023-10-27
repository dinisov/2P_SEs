function BLOCK = splitStack(BLOCK)
%splitStack Split image stack into relevant and blank trial frames
%   Detailed explanation goes here

    % index of blank stimuli/blocks
    indBlank = BLOCK.randomSequence == 5;
    
    % remove blank stimuli from random sequence
    BLOCK.randomSequence = BLOCK.randomSequence(~indBlank);
    
    % construct index of frames after blank blocks
    indBlankStack = reshape(indBlank,[BLOCK.nStimuli length(indBlank)/BLOCK.nStimuli]).';
    indBlankStack = repmat(indBlankStack(:,1),[1 BLOCK.nVol]).';
    
    % split between stimuli image and blank image stacks
    BLOCK.blankImageStack = BLOCK.greenChannel(:,:,indBlankStack(:));
    BLOCK.greenChannel = BLOCK.greenChannel(:,:,~indBlankStack(:));

end