% - I is a 2D image or a 3D volume, where:
% -- 0 represents unlabelled pixels
% -- NaN represents pixels not to label
% -- All other values represent seed labels
% - conn is the pixel connectivity (higher values are slightly slower):
% -- Must be 4 or 8 for 2D image
% -- Must be 6, 18, or 26 for 3D image
%
% ** 2D EXAMPLE **
% hands = imread("hands1-mask.png");
% I = nan(size(hands)); I(hands) = 0; I(122,12) = 1; I(56,30) = 2; 
% I(14,67) = 3; I(14,120) = 4; I(70,201) = 5; I(234,187) = 6;
% I = bwgrowregions(I,8);
% I(isnan(I)) = 0; figure, imagesc(I), axis equal
%
% ** 3D EXAMPLE **
% load("mri.mat",'D')
% mri = squeeze(D);
% nonzero = find(mri); randomNonzero = nonzero(randi(numel(nonzero),6,1));
% I = nan(size(mri)); I(mri>0) = 0; I(randomNonzero) = 1:6;
% I = bwgrowregions(I,6);
% figure, volshow(I)

function I = bwgrowregions(I,conn)

    assert(ismember(ndims(I),[2 3]),'I must be a 2D image or a 3D volume.')
    assert(any(~isnan(I)&I~=0,'all'),'I contains no seed labels.')
    assert((ismatrix(I)&ismember(conn,[4 8]))|...
        (ndims(I)==3&ismember(conn,[6 18 26])),...
        'conn not valid for given I dimensionality.')
    
    % Create base, the subscript offsets to each neighbour
    connMatrix = images.internal.getBinaryConnectivityMatrix(conn);
    connMatrix(2,2,2) = false; % Set the central index to 0
    conn = sum(connMatrix,'all');
    [baseI,baseJ,baseK] = ndgrid(-1:1);
    if ismatrix(I) % 2D
        baseK(:) = 0;
    end
    base = [baseI(connMatrix) baseJ(connMatrix) baseK(connMatrix)];

    % Dilate repeatedly, until no more pixels can be reached
    sz = size(I,[1 2 3]); % Specify dims to handle 2D case
    source = find( ~isnan(I) & I~=0 );
    [I1,I2,I3] = ind2sub(sz,source);
    source = [ I1, I2, I3 ];
    while true

        % Add base indices to source indices to find all neighbours
        neighbours = repelem( source, conn, 1 ) + ...
            repmat( base, size(source,1), 1 );
        % Keep track of source values
        sourceIndices = sub2ind( sz, source(:,1), source(:,2), source(:,3) );
        sourceValues = repelem( I(sourceIndices), conn, 1 );

        % Remove out of image neighbours
        valid = all(neighbours>0 & neighbours<=sz, 2);
        neighbours = neighbours(valid,:);
        sourceValues = sourceValues(valid,:);
        % Remove do-not-label (NaN), already labeled, and non-unique neighbours
        % Note, unique statement implictly selects lowest index source value 
        % when there is an overlap.
        neighboursIndices = sub2ind( sz, ...
            neighbours(:,1), neighbours(:,2), neighbours(:,3) );
        valid = ~isnan(I(neighboursIndices));
        valid = valid & I(neighboursIndices)==0;
        [~,firstOccurence] = unique(neighboursIndices);
        valid = valid & ismember(1:numel(neighboursIndices),firstOccurence)';
        neighbours = neighbours(valid,:);
        neighboursIndices = neighboursIndices(valid,:);
        sourceValues = sourceValues(valid,:);

        % Dilate with source values
        I(neighboursIndices) = sourceValues;

        % Use newly labeled pixels as input for next iteration
        source = neighbours;

        % Exit loop if the region has stopped growing
        if isempty(source)
            if any(I==0)
                warning(['Some pixels were unreachable from the seed ' ...
                    'labels and were left unlabelled.'])
            end
            break
        end

    end
    
end