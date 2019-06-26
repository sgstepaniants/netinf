function parDataSave(fname, noisyData, mat)
    save(fname, 'noisyData', 'mat');
end

function parResultsSave(fname, idx, predMats, est, tprLog, tpr, fprLog, fpr, accLog, acc)
    predMats{idx} = est;
    tprLog(idx) = tpr;
    fprLog(idx) = fpr;
    accLog(idx) = acc;
    save(fname, 'predMats', 'tprLog', 'fprLog', 'accLog');
end
