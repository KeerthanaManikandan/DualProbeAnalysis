function [relativeImportance,rsq] = dominance(X,y)
% Perform relative dominance analysis based on this paper: https://nd.psychstat.org/_media/lab/papers/azen_budescu_2003.pdf
% The codes for the paper are here: https://dominance-analysis.github.io/dominance-analysis/#:~:text=Dominance%20Analysis%20%2D%20The%20Math!,predictor%20to%20the%20regression%20model.&text=The%20measure%20for%20proportion%20of,have%20used%20Pseudo%20R2.&text=It%20can%20bee%20seen%20that,individual%20predictors%20within%20the%20model.
% X: standardized predictors
% y: standardized responses
% Keerthana Manikandan
% October 1, 2025

% Remove rows with missing values
rows = ~isnan(y);
X = X(rows, :);
y = y(rows);

% Determine the number of subset models
modelSize = size(X,2);

r2YX = NaN; r2Diff= {NaN}; 
% Determine all submodels
for iModel = 1:modelSize
    combs = nchoosek(1:modelSize,iModel);

    for iComb = 1:size(combs,1)
        clear mdl

        % Subset models without adding the predictor of interest
        mdl   = fitlm(X(:,combs(iComb,:)),y);
        r2YX(iComb,iModel) = mdl.Rsquared.Ordinary;

        % Subset model after adding the predictor of interest
        for iPred = 1:modelSize
            if ~ismember(iPred,combs(iComb,:))
                mdlTest = fitlm(X(:,[iPred combs(iComb,:)]),y);
                r2Diff{iModel}(iPred,iComb) = mdlTest.Rsquared.Ordinary - r2YX(iComb,iModel);
            else
                r2Diff{iModel}(iPred,iComb) = NaN;
            end
        end
    end
end

% Calculate dominance values
rsq = r2YX(1,modelSize);
meanR2Diff = cell2mat(cellfun(@(x)mean(x,2,'omitnan'),r2Diff,'un',0));

% Relative importance
relativeImportance = mean([meanR2Diff(:,1:modelSize-1) r2YX(1:modelSize,1)],2,'omitnan');

end