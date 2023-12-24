function I = bwgrowregions( I, method )
%BWGROWREGIONS
%   - I is a 2D image or a 3D volume, where:
%   -- 0 represents unlabelled pixels
%   -- NaN represents pixels not to label
%   -- All other values represent seed labels

    arguments
        I { mustBeNumeric, mustBe2Dor3D, mustContainSeed }
        method { mustBeMember( method, { 'chessboard', 'cityblock', ...
            'quasi-euclidean' } ) } = 'quasi-euclidean'
    end
    
    % Base defines the subscript offsets to each neighbour. In 2D, the 
    % connectivity is 8, and in 3D, 26. Precomputed for efficiency, but 
    % calculated as such:
    %
    % conn = 4, 6, 8, 18, or 26
    % connMatrix = images.internal.getBinaryConnectivityMatrix(conn);
    % connMatrix(2,2,2) = false; % Set the central index to 0
    % [ baseI, baseJ, baseK ] = ndgrid( -1 : 1 );
    % if ismatrix( I )
    %     baseK(:) = 0;
    % end
    % base = [ baseI(connMatrix) baseJ(connMatrix) baseK(connMatrix) ];
    %
    if ismatrix( I )
        base = [-1,-1,0;0,-1,0;1,-1,0;-1,0,0;1,0,0;-1,1,0;0,1,0;1,1,0];
    else
        base = [-1,-1,-1;0,-1,-1;1,-1,-1;-1,0,-1;0,0,-1;1,0,-1; ...
                -1,1,-1;0,1,-1;1,1,-1;-1,-1,0;0,-1,0;1,-1,0;-1,0,0; ...
                1,0,0;-1,1,0;0,1,0;1,1,0;-1,-1,1;0,-1,1;1,-1,1; ...
                -1,0,1;0,0,1;1,0,1;-1,1,1;0,1,1;1,1,1];
    end
    connectivity = size( base, 1 );
    if strcmp( method, 'chessboard' )
        baseDistance = ones( connectivity, 1 );
    elseif strcmp( method, 'cityblock' )
        baseDistance = sum( abs( base ), 2 );
    else % 'quasi-euclidean'
        baseDistance = vecnorm( base, 2, 2 );
    end

    sourceIndices = find( ~isnan( I ) & I ~= 0 );
    sourceValues = I(sourceIndices);

    % Check to see if only 1 seed label is used. In this case, the built-in
    % function bwdistgeodesic is faster. bwdistgeodesic can also be used
    % for multiple labels, by running it once for each label, and for each 
    % pixel selecting the label with the minimum distance. However, in my
    % experiments, this was never faster than bwgrowregions, and 
    % substantially slower when many labels are used.
    seedLabels = unique( sourceValues );
    numSeeds = numel( seedLabels );
    if numSeeds == 1
        bw = ~isnan( I );
        mask = I == seedLabels;
        geodesicDist = bwdistgeodesic( bw, mask, method );
        I(~isnan( geodesicDist )) = seedLabels;
        return % EARLY RETURN.
    end
        
    sz = size( I, [1 2 3] ); % Specify dims to handle 2D case.
    geodesicDist = inf( sz );
    geodesicDist(sourceIndices) = 0;
    [ I1, I2, I3 ] = ind2sub( sz, sourceIndices );
    source = [ I1, I2, I3 ];

    % Dilate repeatedly, until no more pixels can be reached.
    while true

        % Add base indices to source indices to find all neighbours.
        neighbours = repelem( source, connectivity, 1 ) + ...
            repmat( base, size(source,1), 1 );
        distToNeighbours = repmat( baseDistance, size(source,1), 1 );
        geodesicDistAtSource = geodesicDist(  repelem( sourceIndices, connectivity, 1 ) );
        geodesicDistAtNeighbours = geodesicDistAtSource + distToNeighbours;
        sourceValues = repelem( sourceValues, connectivity, 1 );

        % Remove out of image neighbours.
        valid = all( neighbours > 0 & neighbours <= sz, 2);
        neighbours = neighbours(valid,:);
        geodesicDistAtNeighbours = geodesicDistAtNeighbours(valid,:);
        sourceValues = sourceValues(valid,:);
        % Remove do-not-label (NaN) pixels and those with a greater 
        % distance than the alternatives.
        neighboursIndices = sub2ind( sz, ...
            neighbours(:,1), neighbours(:,2), neighbours(:,3) );
        [ ~, grName, grNum ] = unique( neighboursIndices );
        grValue = inf( size( grName ) );
        grIdx = zeros( size( grName ) );
        for i = 1 : length( geodesicDistAtNeighbours )
            if geodesicDistAtNeighbours(i) < grValue( grNum(i) )
                grValue( grNum(i) ) = geodesicDistAtNeighbours(i);
                grIdx( grNum(i) ) = i;
            end
        end
        isOverwrite = false( size( neighboursIndices ) );
        isOverwrite( grIdx ) = true;
        preGeodesicDistAtNeighbours = geodesicDist( neighboursIndices );
        isOverwrite = isOverwrite & ...
            geodesicDistAtNeighbours < preGeodesicDistAtNeighbours;
        valid = ~isnan( I(neighboursIndices) ); % Neighbours are traversable (nonNan)
        valid = valid & isOverwrite; % Neighbours are min distance.
        neighbours = neighbours(valid,:);
        neighboursIndices = neighboursIndices(valid,:);
        geodesicDistAtNeighbours = geodesicDistAtNeighbours(valid,:);
        sourceValues = sourceValues(valid,:);

        % Dilate with source values.
        I(neighboursIndices) = sourceValues;
        geodesicDist(neighboursIndices) = geodesicDistAtNeighbours;

        % Use newly labeled pixels as input for next iteration.
        source = neighbours;
        sourceIndices = neighboursIndices;

        % Exit loop if the region has stopped growing.
        if isempty( source )
            if any( I == 0 )
                warning( [ 'Some pixels were unreachable from the seed ' ...
                    'labels and were left unlabelled.' ] )
            end
            break
        end

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