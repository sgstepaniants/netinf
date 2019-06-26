classdef parSave
    methods(Static)
        function parDataSave(fname, noisyData, mat)
            save(fname, 'noisyData', 'mat');
        end
        
        function parPertDataSave(fname, noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat, K)
            save(fname, 'noisyData', 'pertIdx', 'obsIdx', 'pertLength', 'pertTimes', 'mat', 'K');
        end
        
        function parResultsSave(fname, j, k, m, results, est, tpr, fpr, acc)
            predMats = results.predMats;
            predMats{j, k, m} = est;
            
            tprLog = results.tprLog;
            tprLog(j, k, m) = tpr;
            
            fprLog = results.fprLog;
            fprLog(j, k, m) = fpr;
            
            accLog = results.accLog;
            accLog(j, k, m) = acc;
            
            save(fname, 'predMats', 'tprLog', 'fprLog', 'accLog');
        end
    end
end
