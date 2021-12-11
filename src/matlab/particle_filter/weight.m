% This function calcultes the weights for each particle based on the
% observation likelihood
%           S_bar(t)            4XM
%           outlier             1Xn
%           Psi(t)              1XnXM
% Outputs: 
%           S_bar(t)            4XM
function S_bar = weight(S_bar, Psi, outlier)

    % YOUR IMPLEMENTATION
    % Find maximum likelihoods for not outliers 
    Psi_excluding_outliers = Psi(1,find(~outlier),:);

    % Get un-normalised weightingsfdd 
    weights_un = prod(Psi_excluding_outliers,2); % reduces Psi to a 1x1xM from a 1x3xM

    % normalise weights (Weights need to add to one)
    weights = (weights_un/sum(weights_un)); % (1x1xM) sum creates a scalar for normalisation 
    
    % Assign weights to 4th row of S_bar
    S_bar(4,:) = weights;

end
