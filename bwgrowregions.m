function [ I, distanceTransform ] = bwgrowregions( I, distanceMetric )
%BWGROWREGIONS Multi-class region growing for binary images and volumes.
%   Useful for centerline- or skeleton-based segmentation, by using the 
%   skeleton to define seeds. Pixels are given the label of their closest 
%   seed, according to geodesic distance.
% 
%   SYNTAX
%   - labels = bwgrowregions( I )
%   - labels = bwgrowregions( I, distanceMetric )
%   - [ labels, distanceTransform ] = bwgrowregions( __ )
% 
%   INPUTS
%   - I                  Image defining the traversable pixels and seed
%                        locations, given as a 2D or 3D numeric array. Must 
%                        contain at least 1 seed. Note:
%                         -  Values of 0 represent unlabelled pixels which 
%                            can be traversed. These will be labelled if a 
%                            valid path exists, connecting it to any of the 
%                            seeds.
%                         -  Values of NaN represent pixels which cannot be 
%                            traversed and will not be labelled. Paths 
%                            between unlabelled pixels and seeds therefore 
%                            avoid these regions.
%                         -  All other values represent seeds, with the 
%                            value defining that seed's label. The region
%                            grown from a seed is given that seed's label.
%                            Where an unlabelled pixel has a connected path
%                            to more than one seed, the seed with the 
%                            lowest geodesic distance is selected.
%   - distanceMetric     (Optional) Which metric to use for calculating the 
%                        distance transform. This therefore affects the 
%                        precedence of labels for pixels with a path 
%                        connected to more than one seed. Either:
%                         - "chessboard" : Measures the path between pixels 
%                              based on an 8- or 26-connected neighborhood. 
%                              In 2D, pixels whose edges or corners touch 
%                              are 1 unit apart.
%                         - "cityblock" : Measures the path between pixels 
%                              based on a 4- or 6-connected neighborhood. 
%                              In 2D, pixels whose edges touch are 1 unit 
%                              apart, and pixels diagonally touching are 2 
%                              units apart.
%                         - "quasi-euclidean" : (default) Measures the path 
%                              between neighbouring pixels as the 
%                              straight-line distance between them.
%   OUTPUTS
%   - labels             Label matrix of the segmented region of I, 
%                        according to the traversable region, seed 
%                        locations, and seed labels. Returned as a numeric 
%                        array of the same size as I. Each pixel is given 
%                        of the closest seed. Traversable regions without a 
%                        connected path to any seed have a value of 0.
%                        Untraversable regions have a value of NaN. To plot 
%                        the output, first set all NaN values to 0, i.e., 
%                        labels(isnan(labels)) = 0.
%   - distanceTransform  Geodesic distances to the closest seed, returned 
%                        as a numeric array of the same size as I. Seed 
%                        locations have a value of 0, while pixels which 
%                        are not traversable or cannot be reached have a 
%                        value of Inf.
% 
%   EXAMPLES: Please see the file 'examples.mlx' or 'examples.pdf'.
% 
%   Created in 2022b. Compatible with 2019b and later. Compatible with all 
%   platforms. Please cite George Abrahams 
%   https://github.com/WD40andTape/bwgrowregions.
% 
%   See also BWDISTGEODESIC, IMDILATE, BWSKEL, and BWMORPH.

%   Published under MIT License (see LICENSE.txt).
%   Copyright (c) 2023 George Abrahams.
%   - https://github.com/WD40andTape/
%   - https://www.linkedin.com/in/georgeabrahams/
%   - https://scholar.google.com/citations?user=T_xxZLwAAAAJ

    arguments
        I { mustBeNumeric, mustBe2Dor3D, mustContainSeed }
        distanceMetric { mustBeMember( distanceMetric, ...
            { 'chessboard', 'cityblock', 'quasi-euclidean' } ) } ...
            = 'quasi-euclidean'
    end

    % Pad the image so that the neighbours of all pixel are within the 
    % matrix bounds.
    if ismatrix( I )
        I = padarray( I, [1 1], NaN );
    else % 3D (volume).
        I = padarray( I, [1 1 1], NaN );
    end
    [ h, w, d ] = size( I );
    
    % Base defines the offsets to each neighbour in linear indices. In 2D, 
    % the pixel connectivity is 8, and in 3D, 26. baseDistance defines the 
    % distance to each neighbour according to the chosen distance metric. 
    % Precomputed for efficiency. Using linear indices rather than 
    % subscripts avoids a call to sub2ind.
    if d == 1 % 2D (matrix).
        base = [ -h-1, -h, -h+1, -1, +1, +h-1, +h, +h+1 ]';
        if strncmp( distanceMetric, 'chessboard', 2 )
            baseDistance = [ 1, 1, 1, 1, 1, 1, 1, 1 ]';
        elseif strncmp( distanceMetric, 'cityblock', 2 )
            baseDistance = [ 2, 1, 2, 1, 1, 2, 1, 2 ]';
        else % 'quasi-euclidean'
            baseDistance = [ sqrt(2), 1, sqrt(2), 1, 1, ...
                sqrt(2), 1, sqrt(2) ]';
        end
    else % 3D (volume).
        base = [ -h*w-h-1, -h*w-h, -h*w-h+1, -h*w-1, -h*w, -h*w+1, ...
            -h*w+h-1, -h*w+h, -h*w+h+1, -h-1, -h, -h+1, -1, +1, +h-1, ...
            +h, +h+1, +h*w-h-1, +h*w-h, +h*w-h+1, +h*w-1, +h*w, +h*w+1, ...
            +h*w+h-1, +h*w+h, +h*w+h+1 ]';
        if strncmp( distanceMetric, 'chessboard', 2 )
            baseDistance = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]';
        elseif strncmp( distanceMetric, 'cityblock', 2 )
            baseDistance = [ 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, ...
                1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 ]';
        else % 'quasi-euclidean'
            baseDistance = [ sqrt(3), sqrt(2), sqrt(3), sqrt(2), 1, ...
                sqrt(2), sqrt(3), sqrt(2), sqrt(3), sqrt(2), 1, ...
                sqrt(2), 1, 1, sqrt(2), 1, sqrt(2), sqrt(3), sqrt(2), ...
                sqrt(3), sqrt(2), 1, sqrt(2), sqrt(3), sqrt(2), sqrt(3) ]';
        end
    end
    connectivity = size( base, 1 );

    % Identify the seed locations and labels to initialize region growing.
    sourceIndices = find( ~isnan( I ) & I ~= 0 );
    sourceValues = I(sourceIndices);

    % Check to see if only 1 seed label is used. In this case, the built-in
    % function bwdistgeodesic is faster. bwdistgeodesic can also be used
    % for multiple labels, by running it once for each label, and for each 
    % pixel selecting the label with the minimum distance. However, in my
    % experiments, this was always slower than bwgrowregions, 
    % substantially when there are many labels.
    seedLabels = unique( sourceValues );
    numLabels = numel( seedLabels );
    if numLabels == 1
        bw = ~isnan( I );
        mask = I == seedLabels;
        % Find the components connected to the seed locations.
        distanceTransform = bwdistgeodesic( bw, mask, distanceMetric );
        reachablePixels = ~isnan( distanceTransform ) & ...
            ~isinf( distanceTransform );
        I(reachablePixels) = seedLabels;
        return % EARLY RETURN.
    end
    
    % Initialize the distance transform. Seeds have a distance of 0.
    distanceTransform = inf( [ h w d ] );
    distanceTransform(sourceIndices) = 0;
    sourceDist = zeros( numel( sourceIndices ), 1 );

    % Dilate repeatedly, until no more pixels can be reached.
    while true

        % Add base indices to source indices to find all neighbours.
        numSources = numel( sourceIndices );
        sourceIndices = repelem( sourceIndices, connectivity, 1 );
        neighbourIndices = sourceIndices + repmat( base, numSources, 1 );
        
        % Keep track of the label from which each neighbour came.
        % Reuse values from previous iteration, rather than retrieving from
        % matrices as this is faster.
        sourceValues = repelem( sourceValues, connectivity, 1 );
        sourceDist = repelem( sourceDist, connectivity, 1 );

        % For each neighbour, calculate the distance from the seed.
        neighbourDist = sourceDist + repmat( baseDistance, numSources, 1 );

        % Remove neighbours labelled as do-not-label (NaN) and those 
        % with a distance greater than that of the pre-existing label.
        % Doing this now makes later steps more performant.
        isKeep = ~isnan( I(neighbourIndices) ) & ...
            neighbourDist < distanceTransform(neighbourIndices);
        neighbourIndices = neighbourIndices(isKeep,:);
        neighbourDist = neighbourDist(isKeep,:);
        sourceValues = sourceValues(isKeep,:);

        % For non-unique neighbours, select the neighbour with the minimum
        % distance.
        [ ~, groupName, groupNum ] = unique( neighbourIndices );
        groupValue = inf( size( groupName ) );
        groupIdx = zeros( size( groupName ) );
        for iNeighbour = 1 : numel( neighbourIndices )
            if neighbourDist(iNeighbour) < ...
                    groupValue( groupNum(iNeighbour) )
                groupValue( groupNum(iNeighbour) ) = ...
                    neighbourDist(iNeighbour);
                groupIdx( groupNum(iNeighbour) ) = iNeighbour;
            end
        end
        % These neighbours will be written to the distance transform and 
        % labelled on the image.
        isWrite = false( size( neighbourIndices ) );
        isWrite( groupIdx ) = true;
        neighbourIndices = neighbourIndices(isWrite,:);
        neighbourDist = neighbourDist(isWrite,:);
        sourceValues = sourceValues(isWrite,:);

        % Grow regions.
        I(neighbourIndices) = sourceValues;
        distanceTransform(neighbourIndices) = neighbourDist;

        % Use newly labeled pixels as the input for the next iteration.
        sourceIndices = neighbourIndices;
        sourceDist = neighbourDist;

        % Exit the loop if the regions have stopped growing.
        if isempty( sourceIndices )
            if any( I == 0 )
                warning( [ 'Some pixels were unreachable from the ' ...
                    'seed labels and were left unlabelled.' ] )
            end
            break
        end

    end

    % Remove the padding which we added to the image.
    if d == 1 % 2D (matrix).
        I = I( 2:end-1, 2:end-1 );
    else % 3D (volume).
        I = I( 2:end-1, 2:end-1, 2:end-1 );
    end
    
end

%% Validation functions

function mustBe2Dor3D( a )
    if ndims( a ) > 3 % ndims is always greater than 2.
        id = "bwgrowregions:Validators:ArrayNot2Dor3D";
        msg = "Must be either 2D or 3D.";
        throwAsCaller( MException( id, msg ) )
    end
end

function mustContainSeed( a )
    if all( isnan( a ) | a == 0, "all" )
        id = "bwgrowregions:Validators:NoSeeds";
        msg = "Must contain at least 1 seed, i.e., a nonNaN " + ...
            "nonzero element.";
        throwAsCaller( MException( id, msg ) )
    end
end