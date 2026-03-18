function [x, y, area, param] = get_size_image(param, fileName, pathName, um_conv)

    % um_cov represents how many um per pixel.

    if nargin == 3
        param.size_unit = 'pixel';
        um_conv = 1;
    else
        param.size_unit = '?m';
    end
    
    I = imread(fullfile(pathName, fileName),1);
    [x, y] = size(I);
    x = x*um_conv;
    y = y*um_conv;
    area = x*y;
    param.xSize = x;
    param.ySize = y;
    param.pictureArea = area;