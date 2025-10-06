function [importanceVals, R2] = rwa(X,y,predictors)
% Relative Weights Analysis
% Refer:
% https://www.listendata.com/2015/05/relative-importance-weight-analysis.html
% Check this for math: https://www.tandfonline.com/doi/full/10.1080/00273171.2014.905766#d1e243
% for more information on the steps
% X          : (n × p) predictor matrix
% y          : (n × 1) response vector
% predictors : cell array of predictor names (length p)
% KM + ChatGPT
% September 30, 2025

% Step 1. Remove rows with missing values
rows = ~isnan(y);
X = X(rows, :);
y = y(rows);

% Step 2. Correlation matrices
corX  = corr(X);                
corXY = corr(X,y);               

% Step 3. Eigen decomposition of corX
[V, D] = eig(corX);
delta  = sqrt(D);
L      = V * delta * V';        % Johnson’s transformation from X--Z
lambdaSq = L.^2;

% Step 4. Compute beta
beta = L \ corXY; % corY = L*beta, beta = inv(L)*corXY

% Step 5. Compute R^2
R2 = sum(beta.^2);

% Step 6. Relative weights
rawWgt = lambdaSq * (beta.^2);
importance = (rawWgt / R2) * 100;

% Step 7. Store results in a structure
importanceVals.rawWeight  = rawWgt;
importanceVals.scores     = importance; 
importanceVals.betaVals   = beta; 
importanceVals.predictors = predictors; 

end

