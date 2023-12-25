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
        distanceTransform = bwdistgeodesic( bw, mask, method );
        I(~isnan( distanceTransform ) & ~isinf( distanceTransform )) = seedLabels;
        return % EARLY RETURN.
    end
    
    % Initialize the distance transform. Seeds have a distance of 0.
    distanceTransform = inf( [ h w d ] );
    distanceTransform(sourceIndices) = 0;

    % Dilate repeatedly, until no more pixels can be reached.
    while true

        % Add base indices to source indices to find all neighbours.
        numSources = numel( sourceIndices );
        sourceIndices = repelem( sourceIndices, connectivity, 1 );
        neighbourIndices = sourceIndices + repmat( base, numSources, 1 );
        
        % Keep track of the label from which each neighbour came.
        sourceValues = repelem( sourceValues, connectivity, 1 );

        % For each neighbour, calculate the distance from the seed.
        distAtSource = distanceTransform( sourceIndices );
        distAtNeighbours = distAtSource + ...
            repmat( baseDistance, numSources, 1 );

        % Remove neighbours labelled as do-not-label (NaN) and those 
        % with a distance greater than that of the pre-existing label.
        % Doing this now makes later steps more performant.
        isKeep = ~isnan( I(neighbourIndices) ) & ...
            distAtNeighbours < distanceTransform( neighbourIndices );
        neighbourIndices = neighbourIndices(isKeep,:);
        distAtNeighbours = distAtNeighbours(isKeep,:);
        sourceValues = sourceValues(isKeep,:);

        % For non-unique neighbours, select the neighbour with the minimum
        % distance.
        [ ~, groupName, groupNum ] = unique( neighbourIndices );
        groupValue = inf( size( groupName ) );
        groupIdx = zeros( size( groupName ) );
        for iNeighbour = 1 : numel( neighbourIndices )
            if distAtNeighbours(iNeighbour) < ...
                    groupValue( groupNum(iNeighbour) )
                groupValue( groupNum(iNeighbour) ) = ...
                    distAtNeighbours(iNeighbour);
                groupIdx( groupNum(iNeighbour) ) = iNeighbour;
            end
        end
        % These neighbours will be written to the distance transform and 
        % labelled on the image.
        isWrite = false( size( neighbourIndices ) );
        isWrite( groupIdx ) = true;
        neighbourIndices = neighbourIndices(isWrite,:);
        distAtNeighbours = distAtNeighbours(isWrite,:);
        sourceValues = sourceValues(isWrite,:);

        % Grow regions.
        I(neighbourIndices) = sourceValues;
        distanceTransform(neighbourIndices) = distAtNeighbours;

        % Use newly labeled pixels as the input for the next iteration.
        sourceIndices = neighbourIndices;

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