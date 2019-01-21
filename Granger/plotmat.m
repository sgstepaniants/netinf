function plotmat(mat, showText)
    if nargin == 1
        showText = false;
    end
    
    imagesc(mat)
    %colorbar;
    caxis([0 1]);
    
    if showText
        textStrings = num2str(mat(:), '%0.3f');       % Create strings from the matrix values
        textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
        [x, y] = meshgrid(1 : size(mat, 1), 1 : size(mat, 2));  % Create x and y coordinates for the strings
        hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                    'HorizontalAlignment', 'center', 'FontSize', 20, ...
                    'FontWeight', 'bold');
    end
end
