function [duration,cum_duration,arrival,cum_arrival, gof,...
    fit_dissociation,fit_arrival] = track_analysis(track_props,metadata, label)

frames = metadata.totalFrames;
dt = metadata.deltaT;
tot_time = frames*dt;

if nargin<3
    label = 'Full video';
end

% Analysis attachment

try
    temp = groupcounts(struct2table(track_props),'beginning');
    arrival = temp.beginning*dt;
    arrival_count = temp.GroupCount;
    cum_arrival = cumsum(arrival_count);
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = [cum_arrival(1) (cum_arrival(end)-cum_arrival(1))/(arrival(end)-arrival(1))];
    opts.Lower = [0 0];
    [fit_arrival, gof_arrival] = fit(arrival(arrival>tot_time/10), cum_arrival(arrival>tot_time/10), 'poly1');
    adjR2_arrival = gof_arrival.adjrsquare;
    temp_arrival = confint(fit_arrival);
    conf_int_arrival = temp_arrival(:,1);
    figure
    subplot(4,1,1:3);
    plot(fit_arrival,arrival,cum_arrival,'.');
    title(['Attachment - ' label]);
    legend(gca,'off');
    ylabel('Comulative number attached particles');
    txt_attachment = {['Fit function: ' num2str(fit_arrival.p1,3) '*x+' num2str(fit_arrival.p2,3)],...
        ['Arrival rate = ' num2str(fit_arrival.p1,'%.2e') char(177) num2str(fit_arrival.p1-conf_int_arrival(1),'%.2e') ' particles/s'],...
        ['Adjusted R squared: ' num2str(adjR2_arrival,3)]};
    text(tot_time*.95,cum_arrival(1),txt_attachment,'VerticalAlignment','bottom',...
        'HorizontalAlignment','right')
    xlim([0 tot_time]);
    subplot(4,1,4);
    bar(arrival,arrival_count,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
    %ylim([0 0.5*(max(arrival_count(2:end)))])
    ylim([0 max(fit_arrival.p1*dt*5,2.5)]);
    xlim([0 tot_time]);
    xlabel('Time (s)');
    drawnow
    
    gof.association = gof_arrival;
    gof.association.arr_rate_95conf = conf_int_arrival';
    
catch
    warning('Attachment analysis unsuccessful. No data from it will be saved');
    gof.attachement = [];
    arrival = [];
    cum_arrival = [];
    fit_arrival = [];
end


% Analysis detachment

detach_last = floor(max([track_props(:).ending])/2)*dt;
if isfield(track_props,'bleached')
    detach_tracks = track_props([track_props(:).beginning]>0 &...
        [track_props(:).beginning]*dt<detach_last & ~[track_props(:).bleached]);
else
    detach_tracks = track_props([track_props(:).beginning]>0 &...
        [track_props(:).beginning]*dt<detach_last);
end

if length(detach_tracks)<2
    warning('No sufficient detachment events detected to perform the analysis. No plot will be displayed.');
    cum_duration = [];
    duration = [];
    fit_dissociation = [];
else
    
    temp = groupcounts(struct2table(detach_tracks),'duration');
    duration = temp.duration*dt;
    dur_count = temp.GroupCount;
    last = sum(dur_count(duration>detach_last));
    dur_count(duration>detach_last) = [];
    duration(duration>detach_last) = [];
    
    try
        if size(dur_count,1) < 3
            warning('No sufficient detachment events detected to perform the analysis. No plot will be displayed.');
            cum_duration = [];
            gof.dissociation = [];
            fit_dissociation = [];
        else
            dur_count(end) = dur_count(end)+last;
            cum_duration = flip(cumsum(flip(dur_count)));
            ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.StartPoint = [cum_duration(1) 0.01 cum_duration(end)];
            opts.Lower = [0 0 0];
            [fit_dissociation, gof_dissociation] = fit(duration, cum_duration, ft, opts);
            adjR2_dissociation = gof_dissociation.adjrsquare;
            temp_dissociation = confint(fit_dissociation);
            conf_int_dissociation = temp_dissociation(:,2);
            irr_fraction = fit_dissociation.c/(fit_dissociation.a+fit_dissociation.c)*100;
            % Calculation of error propagation from x=c/(a+c): dx =
            % x*sqrt((dx/da)^2+(dx/dc)^2) -> dx = 1.41*ac/(a+c)^2
            irr_error = irr_fraction*fit_dissociation.a*sqrt(2)/(fit_dissociation.a+fit_dissociation.c)^2;
            figure
            subplot(4,1,1:3);
            plot(fit_dissociation,duration,cum_duration,'.');
            title(['Dissociation - ' label]);
            legend(gca,'off');
            ylabel('Bound particles');
            xlim([0 tot_time/2]);
            txt_dissociation = {['Fit function: ' num2str(fit_dissociation.a,3)...
                '*exp(-' num2str(fit_dissociation.b,'%.2e') '*x)+' num2str(fit_dissociation.c,3)],...
                ['Koff = ' num2str(fit_dissociation.b,'%.2e') char(177) num2str(fit_dissociation.b-conf_int_dissociation(1),'%.2e') ' s^-1'],...
                ['Irreversible fraction (IF) = ' num2str(irr_fraction,2) char(177) num2str(irr_error,2) '%'],...
                ['Adjusted R squared: ' num2str(adjR2_dissociation,3)]};
            text(tot_time/2*.95,(cum_duration(1)-cum_duration(end))*.95+cum_duration(end),...
                txt_dissociation,'HorizontalAlignment','right','VerticalAlignment','top')
            subplot(4,1,4);
            bar(duration,dur_count);
            ylim([0.5 1.5*(max(dur_count(1:end-1)))])
            xlim([0 tot_time/2]);
            set(gca, 'YScale', 'log')
            xlabel('Time since attachment (s)');
            drawnow
            
            gof.dissociation = gof_dissociation;
            gof.dissociation.koff_95conf = conf_int_dissociation';
            gof.dissociation.irr_frac_95conf = [irr_fraction-irr_error irr_fraction+irr_error];
            
        end
    catch
        warning('Dissociation analysis unsuccessful. No data from it will be saved');
        gof.dissociation = [];
        fit_dissociation = [];
    end
    
end

