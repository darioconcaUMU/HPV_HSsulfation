function showParticleVideo3DMulti( ...
    tArray, posArray, attachmentArray, anchorPos, R )
% SHOWPARTICLEVIDEO3DMULTI
%
% Animate a 3D Brownian tether motion of a sphere (radius R), with multiple tethers.
% Also draws the projection of the sphere center onto the XY, XZ, and YZ planes.
%
% Inputs:
%   tArray        : [nFrames x 1], the time array
%   posArray      : [nFrames x 3], the sphere center coordinates at each frame
%   attachmentArray : [nFrames x 3 x nTethers], sphere anchor points for each tether
%   anchorPos     : [3 x nTethers], plane anchor positions (often z=0, but not required)
%   R             : scalar, the sphere radius
%
% The script:
%   1) Plots the sphere, lines from plane anchors -> sphere anchors,
%   2) Adds a point for each sphere anchor,
%   3) Shows the sphere center's projection onto XY, XZ, and YZ planes,
%   4) Saves an AVI "multiTethers3D_excl.avi".
%
% Example usage:
%   showParticleVideo3DMulti(tArray, posArray, attachmentArray, anchorPos, 25e-9);

    % Convert anchorPos from (3 x nTethers) to expected shape if needed
    % user said: anchorPos = anchorPos';  done outside or inside, depending on usage
    % anchorPos = anchorPos';

    nFrames = length(tArray);
    if size(posArray,1)~=nFrames
        error('posArray rows must match nFrames');
    end
    if size(attachmentArray,1)~=nFrames
        error('attachmentArray first dimension must match nFrames');
    end

    % number of tethers:
    nTethers = size(attachmentArray, 3);

    % We skip frames so we only show about ~100 frames
    nPlotFrames = 100;
    skip = max(1, floor(nFrames / nPlotFrames));

    % --- Figure out bounding box for the axis ---
    % Flatten sphere anchor coords:
    allXAnch = reshape(attachmentArray(:,1,:), [], 1);
    allYAnch = reshape(attachmentArray(:,2,:), [], 1);
    allZAnch = reshape(attachmentArray(:,3,:), [], 1);

    allX = [posArray(:,1); allXAnch];
    allY = [posArray(:,2); allYAnch];
    allZ = [posArray(:,3); allZAnch];

    minX = min(allX) - R;  maxX = max(allX) + R;
    minY = min(allY) - R;  maxY = max(allY) + R;
    minZ = min(allZ) - R;  maxZ = max(allZ) + R;

    if minZ > 0, minZ = 0; end  % ensure we see plane z=0 if relevant

    % --- Prepare figure ---
    fig = figure('Name','3D Multi-Tethers Animation','Color','white');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    xlabel(ax,'x [m]');  ylabel(ax,'y [m]');  zlabel(ax,'z [m]');
    title(ax,'Sphere + Multiple Tethers (3D)');

    xlim(ax,[minX, maxX]);
    ylim(ax,[minY, maxY]);
    zlim(ax,[0, maxZ]);
    view(ax, 3);  % 3D perspective view

    % Plot each plane anchor as a black square
    for k=1:nTethers
        plot3(ax, anchorPos(k,1), anchorPos(k,2), anchorPos(k,3), ...
            'ks','MarkerFaceColor','k');
    end

    % We'll create line objects for each tether, and a marker for each tether's anchor
    tetherLines = gobjects(nTethers,1);
    sphereAnchors = gobjects(nTethers,1);

    % Initialize them for frame=1
    for k=1:nTethers
        xA = attachmentArray(1,1,k);
        yA = attachmentArray(1,2,k);
        zA = attachmentArray(1,3,k);

        px = [anchorPos(k,1), xA];
        py = [anchorPos(k,2), yA];
        pz = [anchorPos(k,3), zA];

        tetherLines(k) = plot3(ax, px, py, pz, 'b-','LineWidth',1.5);
        sphereAnchors(k)= plot3(ax, xA, yA, zA, 'ko','MarkerFaceColor','k');
    end

    % sphere as a surface (parametric)
    sphereRes = 20;
    [theta,phi] = meshgrid(linspace(0,pi,sphereRes), linspace(0,2*pi,sphereRes));
    Xs = zeros(size(theta));
    Ys = zeros(size(theta));
    Zs = zeros(size(theta));

    % set up initial sphere surface at frame=1
    rC = posArray(1,:);  % center
    for idx=1:numel(theta)
        Xs(idx) = rC(1) + R*sin(theta(idx))*cos(phi(idx));
        Ys(idx) = rC(2) + R*sin(theta(idx))*sin(phi(idx));
        Zs(idx) = rC(3) + R*cos(theta(idx));
    end
    sphSurf = surf(ax, Xs, Ys, Zs, 'FaceColor','r','FaceAlpha',0.3,...
                   'EdgeColor','none');

    % --- Add markers for the center's 2D projections on XY, XZ, YZ planes
    % We'll define them as squares or circles with different colors:
    xC = rC(1); yC=rC(2); zC=rC(3);

    % XY-plane: z=0 => (xC, yC, 0)
    projXY = plot3(ax, xC, yC, 0, 'rs','MarkerFaceColor','r','MarkerSize',6, ...
        'DisplayName','ProjXY');
    % XZ-plane: y=maxY => (xC, maxY, zC)
    projXZ = plot3(ax, xC, maxY, zC, 'gs','MarkerFaceColor','g','MarkerSize',6, ...
        'DisplayName','ProjXZ');
    % YZ-plane: x=maxX => (maxX, yC, zC)
    projYZ = plot3(ax, maxX, yC, zC, 'bs','MarkerFaceColor','b','MarkerSize',6, ...
        'DisplayName','ProjYZ');

    % --- Setup VideoWriter ---
    v = VideoWriter('multiTethers3D_real_test2.avi','Motion JPEG AVI');
    v.FrameRate = 10;
    open(v);

    % --- Animation Loop ---
    for i = 1:skip:nFrames

        % Sphere center
        xC = posArray(i,1);
        yC = posArray(i,2);
        zC = posArray(i,3);

        % 1) Update sphere
        for idx=1:numel(theta)
            Xs(idx) = xC + R*sin(theta(idx))*cos(phi(idx));
            Ys(idx) = yC + R*sin(theta(idx))*sin(phi(idx));
            Zs(idx) = zC + R*cos(theta(idx));
        end
        set(sphSurf, 'XData',Xs, 'YData',Ys, 'ZData',Zs);

        % 2) Update each tether line and sphere anchor
        for k=1:nTethers
            xA = attachmentArray(i,1,k);
            yA = attachmentArray(i,2,k);
            zA = attachmentArray(i,3,k);

            px = [anchorPos(k,1), xA];
            py = [anchorPos(k,2), yA];
            pz = [anchorPos(k,3), zA];

            set(tetherLines(k), 'XData',px, 'YData',py, 'ZData',pz);
            set(sphereAnchors(k), 'XData',xA, 'YData',yA, 'ZData',zA);
        end

        % 3) Update the projection markers
        set(projXY, 'XData', xC, 'YData', yC, 'ZData', 0);
        set(projXZ, 'XData', xC, 'YData', maxY,  'ZData', zC);
        set(projYZ, 'XData', maxX,   'YData', yC, 'ZData', zC);

        % Update title with current time
        curTime = tArray(i);
        title(ax, sprintf('Multi-Tethers (t=%.2g s)', curTime));

        drawnow;

        % Capture frame
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    % close video
    close(v);
    pause(1);  % hold final frame
end
