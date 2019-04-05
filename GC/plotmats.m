function plotmats(true, pred, showText)
    if nargin == 1
        showText = false;
    end
    
    subplot(1,2,1)
    plotmat(true, showText)
    subplot(1,2,2)
    plotmat(pred, showText)
end
