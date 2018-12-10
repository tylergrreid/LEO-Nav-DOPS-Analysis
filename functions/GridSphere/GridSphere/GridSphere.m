function [latGridInDegrees, longGridInDegrees] = GridSphere(numPoints)
% Generates a geodesic grid, a set of nearly evenly spaced points on a sphere.
% Points are sorted in order of increasing latitude. Ties between latitudes that
% differ only by floating point errors are broken by a secondary sort in order
% of increasing longitude. Tesselates a sphere with triangles. Starts with an
% icosahedron circumscribed by a unit sphere. Divides the triangles into 4 new
% equal triangles. Projects the new vertices onto the unit sphere. Repeats until
% reaching the desired resolution. Refer to
% http://www.scidacreview.org/0904/images/hardware04.jpg for a graphical
% depiction of the algorithm.
%
% -----------------------------------------------------------------------
% N. Teanby   13-01-04    Original IDL (Interactive Data Language) code
% (available at http://www.atm.ox.ac.uk/user/teanby/software.html#icos)
%
% Kurt von Laven    18-09-13    Ported simplified version to GNU Octave.
% -----------------------------------------------------------------------
%
%    usage: [latGridInDegrees, longGridInDegrees] = GridSphere(numPoints)

    PHI = GoldenRatio();
    % Calculate the golden ratio.

    NUM_TRIANGLES_IN_AN_ICOSAHEDRON = 20;
    % The number of faces in an icosahedron.

    NUM_DIMENSIONS = 3;
    % The dimensionality of the sphere.

    NUM_VERTICES_IN_A_TRIANGLE = 3;
    % The number of vertices in a triangle.

    NUM_NEW_TRIANGLES_PER_TRIANGLE = 4;
    % The number of triangles produced 

    VERTICES = [0, PHI, 1; 0, -PHI, 1; 0, PHI, -1; 0, -PHI, -1; 1, 0, PHI; ...
        -1, 0, PHI; 1, 0, -PHI; -1, 0, -PHI; PHI, 1, 0; -PHI, 1, 0; PHI, -1, ...
        0; -PHI, -1, 0] / norm([1, PHI]);
    % Find the 12 vertices of an icosahedron inscribed in a unit sphere.

    As = [2; 5; 9; 7; 11; 4; 6; 2; 1; 5; 3; 9; 8; 7; 12; 12; 6; 1; 3; 10];
    Bs = [4; 11; 5; 9; 7; 2; 12; 5; 6; 9; 1; 7; 3; 4; 8; 6; 1; 3; 10; 8];
    Cs = [11; 2; 11; 11; 4; 12; 2; 6; 5; 1; 9; 3; 7; 8; 4; 10; 10; 10; 8; 12];
    triangles = reshape([VERTICES(As, :), VERTICES(Bs, :), VERTICES(Cs, :)], ...
                            NUM_TRIANGLES_IN_AN_ICOSAHEDRON, NUM_DIMENSIONS, ...
                            NUM_VERTICES_IN_A_TRIANGLE);
    % Find the vertices of the starting triangles. Each row represents a
    % triangle. Columns represent x, y, and z coordinates. Depth represents the
    % A, B, and C triangle vertices. Ensure that each vertex is an A, B, and C
    % vertex in at least one triangle.

    numDivisions = ceil(Logarithm((numPoints - 2) / 25, ...
                                            NUM_NEW_TRIANGLES_PER_TRIANGLE));
    % Determines the number of divisions that will result in a grid with as
    % close to the number of points as the user requested as possible.

    for i = 1:numDivisions
        triangles = DivideTriangles(triangles, NUM_DIMENSIONS, ...
                    NUM_VERTICES_IN_A_TRIANGLE, NUM_NEW_TRIANGLES_PER_TRIANGLE);
        % Split each triangle into 4 smaller triangles.
    end

    [latGridInDegrees, longGridInDegrees] = ...
                                    Cartesian2Spherical(triangles(:, 1, 1), ...
                                    triangles(:, 2, 1), triangles(:, 3, 1));
    % Calculate latitudes and longitudes based on Cartesian points. All the
    % vertices in the grid are shared by a number of triangles. Only the unique
    % values concern us. DivideTriangles and the original icosahedron
    % construction guarantee that every vertex is an A vertex in at least one
    % triangle, so only A vertices are considered.

    [latGridInDegrees, longGridInDegrees] = ...
                        RemoveEqualPairs(latGridInDegrees, longGridInDegrees);
    % Even though B and C vertices were discarded, many identical points remain,
    % so they must be removed. Sort first by increasing latitude, breaking ties
    % between latitudes differing only because of floating point imprecision by
    % sorting by increasing longitude. RemoveEqualPairs requires that for every
    % cluster of coordinates with equal latitudes, the greatest longitude does
    % not equal the least longitude in the next cluster. This property is
    % satisfied because, for all clusters of equal latitude excluding -90 and
    % +90 degrees, there is at least one negative and one positive longitude.
    % Furthermore, should there be any points at the poles, their longitudes are
    % set to 0 degrees.

end
