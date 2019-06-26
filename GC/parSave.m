classdef parSave
    methods(Static)
        function parDataSave(fname, noisyData, mat)
            save(fname, 'noisyData', 'mat');
        end
        
        function parResultsSave(fname, j, k, m, predMats, est, tprLog, tpr, fprLog, fpr, accLog, acc)
            predMats{j, k, m} = est;
            tprLog(j, k, m) = tpr;
            fprLog(j, k, m) = fpr;
            accLog(j, k, m) = acc;
            save(fname, 'predMats', 'tprLog', 'fprLog', 'accLog');
        end
    end
end
