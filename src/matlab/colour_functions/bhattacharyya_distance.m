% Adapted from below
%
% BHATTACHARYYA(histogram1, histogram2)
% compute the BHATTACHARYYA distance between 2 histograms
% where each histogram is a 1xN vector
% 
% Based on the estimation presented in 
% "Real-Time Tracking of Non-Rigid Objects using Mean Shift"
%  D. Comaniciu, V. Ramesh & P. Meer (IEEE CVPR 2000, 142-151)
%
% N.B. both histograms must be normalised
% (i.e. bin values lie in range 0->1 and SUM(bins(i)) = 1
%       for i = {histogram1, histogram2} )
%
% Author / Copyright: T. Breckon, October 2005.
% School of Informatics, University of Edinburgh.
% License: http://www.gnu.org/licenses/gpl.txt

function bdist = bhattacharyya_distance(histogram1, histogram2)

    % Calculate the bhattacharyya co-efficient
     bcoeff = sum(sqrt(histogram1.*histogram2))/3;

    % Calculate the bhattacharyya distances
    bdist = sqrt(1 - bcoeff);

end