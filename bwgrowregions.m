function I = bwgrowregions( I, method )
%BWGROWREGIONS Multi-class, distance-based segmentation of binary images.
%   - I is a 2D image or a 3D volume, where:
%   -- 0 represents unlabelled pixels
%   -- NaN represents pixels not to label
%   -- All other values represent seed labels

    arguments
        I { mustBeNumeric, mustBe2Dor3D, mustContainSeed }
        method { mustBeMember( method, { 'chessboard', 'cityblock', ...
            'quasi-euclidean' } ) } = 'quasi-euclidean'
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
        if strncmp( method, 'chessboard', 2 )
            baseDistance = [ 1, 1, 1, 1, 1, 1, 1, 1 ]';
        elseif strncmp( method, 'cityblock', 2 )
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
        if strncmp( method, 'chessboard', 2 )
            baseDistance = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]';
        elseif strncmp( method, 'cityblock', 2 )
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

    sourceIndices = find( ~isnan( I ) & I ~= 0 );
    sourceValues = I(sourceIndices);

    % Check to see if only 1 seed label is used. In this case, the built-in
    % function bwdistgeodesic is faster. bwdistgeodesic can also be used
    % for multiple labels, by running it once for each label, and for each 
    % pixel selecting the label with the minimum distance. However, in my
    % experiments, this was never faster than bwgrowregions, and 
    % substantially slower when many labels are used.
    seedLabels = unique( sourceValues );
    numLabels = numel( seedLabels );
    if numLabels == 1
        bw = ~isnan( I );
        mask = I == seedLabels;
        geodesicDist = bwdistgeodesic( bw, mask, method );
        I(~isnan( geodesicDist )) = seedLabels;
        return % EARLY RETURN.
    end
    
    geodesicDist = inf( [ h w d ] );
    geodesicDist(sourceIndices) = 0;

    % Dilate repeatedly, until no more pixels can be reached.
    while true

        % Add base indices to source indices to find all neighbours.
        numSources = numel( sourceIndices );
        sourceIndices = repelem( sourceIndices, connectivity, 1 );
        neighbourIndices = sourceIndices + repmat( base, numSources, 1 );
        distToNeighbours = repmat( baseDistance, numSources, 1 );
        geodesicDistAtSource = geodesicDist( sourceIndices );
        geodesicDistAtNeighbours = geodesicDistAtSource + distToNeighbours;
        sourceValues = repelem( sourceValues, connectivity, 1 );

        % Stop assessing trivial cases: do-not-label (NaN) pixels and those 
        % with a greater distance than that of the pre-existing label.
        isKeep = ~isnan( I(neighbourIndices) ) & ...
            geodesicDistAtNeighbours < geodesicDist( neighbourIndices );
        neighbourIndices = neighbourIndices(isKeep,:);
        geodesicDistAtNeighbours = geodesicDistAtNeighbours(isKeep,:);
        sourceValues = sourceValues(isKeep,:);
        % For non-unique neighbours, select the neighbour with the minimum
        % distance.
        [ ~, grName, grNum ] = unique( neighbourIndices );
        grValue = inf( size( grName ) );
        grIdx = zeros( size( grName ) );
        for i = 1 : length( geodesicDistAtNeighbours )
            if geodesicDistAtNeighbours(i) < grValue( grNum(i) )
                grValue( grNum(i) ) = geodesicDistAtNeighbours(i);
                grIdx( grNum(i) ) = i;
            end
        end
        isOverwrite = false( size( neighbourIndices ) );
        isOverwrite( grIdx ) = true;
        neighbourIndices = neighbourIndices(isOverwrite,:);
        geodesicDistAtNeighbours = geodesicDistAtNeighbours(isOverwrite,:);
        sourceValues = sourceValues(isOverwrite,:);

        % Dilate with source values.
        I(neighbourIndices) = sourceValues;
        geodesicDist(neighbourIndices) = geodesicDistAtNeighbours;

        % Use newly labeled pixels as input for next iteration.
        sourceIndices = neighbourIndices;

        % Exit loop if the region has stopped growing.
        if isempty( sourceIndices )
            if any( I == 0 )
                warning( [ 'Some pixels were unreachable from the seed ' ...
                    'labels and were left unlabelled.' ] )
            end
            break
        end

    end

    % Remove padding.
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
        id = "bwgrowregions:Validators:NoSeedLabels";
        msg = "Must contain at least 1 seed label, i.e., a nonNaN " + ...
            "nonzero element.";
        throwAsCaller( MException( id, msg ) )
    end
end