function [sortedOnes, sortedTwos] = SortPairsByOnes(ones, twos)
% Sorts both vectors in order of increasing ones.
% Assumes both vectors are the same length. See also the built-in function
% sortrows.
%
%    usage: [sortedOnes, sortedTwos] = SortPairsByOnes(ones, twos)
%
%    example:    [sortedOnes, sortedTwos] = SortPairsByOnes([-5, 7, 6, -10], ...
%                                                            [0; 22; -67; 8])
%                sortedOnes = [-10, -5, 6, 7]
%                sortedTwos = [8; 0; -67; 22]

    [sortedOnes, sortIndices] = sort(ones);
    sortedTwos = twos(sortIndices);

end
