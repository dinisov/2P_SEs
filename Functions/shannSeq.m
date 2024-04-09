function [shannTab,shannAll] = shannSeq(exLabels)
%shannSeq - Convert sequential dependencies labels to Shannon entropy values
%   Takes 'exLabels' as an input (cell array of strings)
%   If argument empty will calculate Shannon from first principles (based on stored labels)

    disp(['Calculating Shannon entropy of sequences'])
    
    if isempty(exLabels) == 1
        disp(['Performing pre-label processes'])
        addpath('D:\group_vanswinderen\Dinis\Scripts\Indexes and legends\');
        load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
        ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
        exLabels = [];
        for s = 1:size(binomial_x_labels_latex,2)
            exLabels{s} = binomial_x_labels_latex{s}(ind_horiz);
        end
    end
    
    shannTab = []; %Will store both original sequence and Shannon entropy values
    shannAll = []; %Just the Shannon values, for certain operations
    for s = 1:size(exLabels,2)
        temp = exLabels{s};
        inds{1} = find( temp == 'R' );
        inds{2} = find( temp == 'A' );
        shannTab{1,s} = nan( 1 , size(temp,2) );
        shannTab{1,s}( inds{1} ) = 1;
        shannTab{1,s}( inds{2} ) = 2;
        ent = abs( wentropy(shannTab{1,s},'shannon') ); %Abs to make plots a bit easier
        shannTab{2,s} = ent;
        shannAll(s) = ent;
        disp([num2str(s),' - ',num2str(shannTab{1,s}),' - ',num2str(shannTab{2,s})])
    end
    shannAll = round(shannAll,5); %Round, to make plots a little easier
end