function analyze3DDistribution(posArray)
    % PLOTXYDISTRIBUTION Takes a 3 x n matrix (x,y,z) and plots the 
    % distribution of (x,y) as a 2D histogram.
    %
    %   positions is a 3 x n matrix, where each column is [x; y; z].
    %   This function will ignore z and just focus on x and y,
    %   plotting the distribution using histcounts2 and imagesc.

    %% 1) Basic checks
    if size(posArray, 1) ~= 3
        posArray = posArray';
    end

    %% 2) Extract x and y
    x = posArray(1,:);
    y = posArray(2,:);

    %% 3) Compute 2D histogram
    nbins = 50;  % Adjust the number of bins if needed
    [counts, xEdges, yEdges] = histcounts2(x, y, nbins);

    %% 4) Plot the 2D histogram as an image
    figure;
    % imagesc flips the axes, so we'll pass the edges directly
    % Note: histcounts2 returns edges for bin boundaries
    imagesc(xEdges, yEdges, counts');
    axis xy;        % Correct the axis orientation
    colorbar;       % show a color scale
    xlabel('X');
    ylabel('Y');
    title('2D Distribution (X vs Y)');

    %% 5) Compute the (x,y) center of mass
    x_cm = mean(x);
    y_cm = mean(y);

    hold on
    plot(x_cm, y_cm, 'o')
    hold off

    %% 6) Compute radial distances from the (x_cm, y_cm)
    dx = x - x_cm;
    dy = y - y_cm;
    r = sqrt(dx.^2 + dy.^2);
    

    %% 7) Plot a histogram of the radial distribution
    figure;
    [data, edges] = histcounts(r, 40);   % 30 bins as an example
    edgesCentre = edges(2:end)-(edges(2)-edges(1));
    areaNorm = edges(2:end).^2-edges(1:end-1).^2;
    dataNorm = data./areaNorm/sum(data./areaNorm);
    plot(edgesCentre, dataNorm);
    xlabel('Radial Distance from Center of Mass (xy-plane)');
    ylabel('Counts');
    title('Radial Distance Distribution');
end
