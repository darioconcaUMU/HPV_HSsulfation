
function plot_selected_trajectories(track,kon,binders,n_tracks,color_axis)

if nargin < 5,  color_axis = [0 0 0];  end

sel_track = squeeze(squeeze(track(kon, binders,:)));
side_plot = sqrt(n_tracks);
plot_subplot(floor(side_plot),ceil(side_plot),sel_track, color_axis);

end


function plot_subplot(n1,n2,track,color_axis)

if nargin < 4,  color_axis = [0 0 0];  end

fig = figure('Units','pixels', ...          % full-screen figure
             'OuterPosition',[0 0 500 500]);

tl = tiledlayout(fig,n1,n2, ...
                 'TileSpacing','compact', ...   % or 'none'
                 'Padding','compact');          % or 'none'

for i = 1:n1*n2
    ax = nexttile;                              % replaces subplot
    traj = track{i}.positions;
    x = traj(:,1) - mean(traj(:,1));
    y = traj(:,2) - mean(traj(:,2));
    t = 1:numel(x);

    plot(ax,x,y,'-k','LineWidth',0.2); hold on
    scatter(ax,x,y,2,t,'filled');

    axis(ax,'equal')
    clim(ax,[0 500])
    set(ax,'XColor',color_axis,'YColor',color_axis, ...
           'XTick',[],'YTick',[])
    xlim(ax,[-1 1]); ylim(ax,[-1 1]);
end
end

% function plot_subplot(n1,n2,track,color_axis)
% 
% figure
% if nargin<4
%     color_axis = [0 0 0];
% end
% 
% for i=1:n1*n2
%     subplot(n1,n2,i)
%     traj = track{i}.positions;
%     x = traj(:,1)-mean(traj(:,1));          % x-coordinates  (nm, µm, px …)
%     y = traj(:,2)-mean(traj(:,2));          % y-coordinates
%     t = 1:numel(x);
% 
%     % Scatter plot coloured by time
%     plot(x,y,'-k','LineWidth',0.2);
%     hold on;
%     scatter(x, y, 1, t);   
%     clim([0 500]); 
% 
%     axis equal
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
%     xlim([-1 1]);
%     ylim([-1 1]);
% 
%     ax = gca; 
%     ax.XColor = color_axis;
%     ax.YColor = color_axis;
% 
% end
% end