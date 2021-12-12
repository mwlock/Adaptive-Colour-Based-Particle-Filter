% This function performs systematic re-sampling
% Inputs:   
%           S_bar(t):       4XM
% Outputs:
%           S(t):           4XM
function S = systematic_resample(S_bar)
	
    global M % number of particles 
    
    % YOUR IMPLEMENTATION
    S = zeros(4,M);

    CDF = cumsum(S_bar(4,:));

    r_0 = rand/M;

    for m = 1:M
        index = find(CDF >= (r_0 + (m-1)/M), 1, 'first');
        S(:,m) = [S_bar(1:3,index); 1/M];
    end
end