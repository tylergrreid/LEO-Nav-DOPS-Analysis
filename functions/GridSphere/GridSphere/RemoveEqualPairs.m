function [prunedOnes, prunedTwos] = RemoveEqualPairs(ones, twos)
% Returns only the distinct ordered pairs, removing identical pairs.
% Both elements of a pair must be within rounding errors typical of the double
% class of the corresponding elements of another pair for the pairs to be
% considered equal. Sorts output in order of increasing ones. Breaks ties among
% points with ones that differ only by rounding errors with a secondary sort in
% order of increasing twos.
%
%    usage: [prunedOnes, prunedTwos] = RemoveEqualPairs(ones, twos)

    [sortedOnes, sortedTwos] = SortPairsByOnes(ones, twos);
    % Sort the pairs in order of ascending ones.

    uniqueOneIndices = DistinctFromPrevious(sortedOnes);
    % Find the lowest index at which each unique value appears in sortedOnes.

    numUniqueOnes = numel(uniqueOneIndices);

    for i = 1:(numUniqueOnes - 1)

        cluster = RangeToNext(uniqueOneIndices, i);
        [sortedTwos(cluster), sortIndices] = sort(sortedTwos(cluster));
        sortIndicesAsInts = uint32(sortIndices);
        sortedOnes(cluster) = sortedOnes(uniqueOneIndices(i) + ...
                                         sortIndicesAsInts - 1);
        % Within each cluster of equal ones, sort in order of ascending twos.
        % Note that this operation doesn't change the validity of
        % uniqueOneIndices.

    end

    finalCluster = uniqueOneIndices(numUniqueOnes):numel(sortedTwos);
    [sortedTwos(finalCluster), sortIndices] = sort(sortedTwos(finalCluster));
    sortIndicesAsInts = uint32(sortIndices);
    sortedOnes(finalCluster) = sortedOnes(uniqueOneIndices(numUniqueOnes) ...
                                          + sortIndicesAsInts - 1);
    % The last cluster must be handled specially to avoid indexing out of
    % bounds.

    uniqueTwoIndices = DistinctFromPrevious(sortedTwos);

    uniquePairIndicesWithDuplicates = ...
        sort([uniqueOneIndices; uniqueTwoIndices]);
    uniquePairIndices = uniquePairIndicesWithDuplicates(...
            DistinctFromPrevious(uniquePairIndicesWithDuplicates));
    % Once the pairs are sorted within clusters, the unique pairs
    % are those with either a one that differs from the previous
    % one or a two that differs from the previous two.

    prunedOnes = sortedOnes(uniquePairIndices);
    prunedTwos = sortedTwos(uniquePairIndices);
    % Create arrays consisting only of the unique pairs.

end
