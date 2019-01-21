function confusionMat = ConfusionMatrix(mats, predMats)
    [nvars, ~, numMats] = size(mats);
    
    confusionMat = zeros(2^(nvars^2 - nvars));
    for j = 1 : numMats
        mat = mats(:, :, j);
        matIndex = indexMat(mat);
        predMat = predMats(:, :, j);
        predMatIndex = indexMat(predMat);
        
        confusionMat(matIndex, predMatIndex) = confusionMat(matIndex, predMatIndex) + 1;
    end
end
