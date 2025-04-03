
function DistAng =CC_DistAng2(SPcroscorr)

aCorr =  SPcroscorr;

    halfSize = ceil(size(aCorr)/2);
    
    half_height = halfSize(1);
    half_width = halfSize(2);
    aCorrRad = min(halfSize);
    aCorrSize = size(aCorr);

maxValue = max(max(aCorr));

    if maxValue ~= 1
        aCorr = aCorr / maxValue;
    end
    
     [rr, cc] = meshgrid(1:size(aCorr, 2), 1:size(aCorr, 1));

    % Define iteration radius step size for the gridness score
    radSteps = 1:aCorrRad;

    numSteps = length(radSteps);

    GNS = zeros(numSteps, 2);

    mainCircle = sqrt((cc - half_height).^2 + (rr - half_width).^2);

    % Define expanding ring of autocorrellogram and do x30 correlations
    for i = 1:numSteps
        ind = mainCircle < radSteps(i);
        tempCorr = reshape(aCorr(ind), 1, [])';

        GNS(i, 2) = radSteps(i);
    end

    gscoreLoc = aCorrRad-5;
   
 bestCorr = (mainCircle < radSteps(gscoreLoc) * 1.25) .* aCorr;
    
    bestCorr(isnan(bestCorr)) = 0;
    regionalMaxMap = imregionalmax(bestCorr, 4);
    se = strel('square', 3);
    im2 = imdilate(regionalMaxMap, se); % dilate map to eliminate fragmentation
    cc = bwconncomp(im2, 8);
    stats = regionprops(cc, 'Centroid');

allCoords = [stats(:).Centroid];
clear centresOfMass;
    centresOfMass(:, 1) = allCoords(1:2:end);
    centresOfMass(:, 2) = allCoords(2:2:end);
       
    centerX = size(aCorr,2)/2;
    centerY = size(aCorr,1)/2; 
    
    DistAng.DistanceX = centresOfMass(:, 1)-centerX; 
    DistAng.DistanceY = centresOfMass(:, 2)-centerY; 
    
    DistAng.DistanceXY = sqrt((centresOfMass(:, 1)-centerX).^2 + (centresOfMass(:, 2)-centerY).^2);    
    [~, sortInd] = sort(DistAng.DistanceXY);
    
    % %     stats = stats(sortInd);
    DistAng.DistanceXY = DistAng.DistanceXY(sortInd);
    DistAng.DistanceX = DistAng.DistanceX(sortInd);
    DistAng.DistanceY = DistAng.DistanceY(sortInd);
    
    
    DistAng.orientationXY = rad2deg(atan2(centresOfMass(:, 2)-centerY, centresOfMass(:, 1)-centerX));
    DistAng.orientationXY = DistAng.orientationXY(sortInd);
    DistAng.orientationXY = DistAng.orientationXY(1:1);
    DistAng.DistanceXY = DistAng.DistanceXY(1:1);
    DistAng.DistanceX = DistAng.DistanceX(1:1);
    DistAng.DistanceY = DistAng.DistanceY(1:1);
    
    
end



  











