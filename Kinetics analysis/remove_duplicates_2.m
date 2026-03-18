function [tracks, track_props] = remove_duplicates_2(tracks, track_props, distance)

% This function merges spots that are present symultaneously and are very
% close to eachoter. This is donebecause the software can double spots with
% high intensity.

disp('Starting process to remove duplicate spots');

partition = 20;

% Extracts data from track_props struction for easier handling
x_tracks = [track_props(:).pos_mean_x];
y_tracks = [track_props(:).pos_mean_y];
start_tracks = [track_props(:).beginning];
end_tracks = [track_props(:).ending];

% Get the maximum position in x and y direction
xlim = ceil(max(x_tracks+1));
ylim = ceil(max(y_tracks+1));
track_num = 1:length(track_props);
segment_x = ceil(xlim/partition);
segment_y = ceil(ylim/partition);

% To avoid problems with running out of memory, similarly to
% blink_detection_4 the script splits the tracks according to their x
% position.

ind_del = [];

for ii=1:partition
    for jj=1:partition
        
        if ii == 1 && jj==2
            start_time = tic;
        end
        if ii == 1 && jj==5
            time_loop = toc(start_time);
            predicted_time = time_loop*partition^2/3;
            if predicted_time>10
                disp(['The extimated completion time in seconds is: ' ...
                num2str(predicted_time)]);
            end
            if predicted_time>600
                disp('Relax and go get a coffee, I''m working hard on it');
            end
        end
        lim_xlow = max((ii-1)*segment_x-2*distance, 1);
        lim_ylow = max((jj-1)*segment_y-2*distance, 1);
        lim_xhigh = min(ii*segment_x+2*distance, xlim);
        lim_yhigh = min(jj*segment_y+2*distance, ylim);
        sel_seg = x_tracks>lim_xlow & x_tracks<=lim_xhigh &...
            y_tracks>lim_ylow & y_tracks<=lim_yhigh;
        index_seg = track_num(sel_seg);
        x_part = x_tracks(sel_seg);
        y_part = y_tracks(sel_seg);
        start_seg = start_tracks(sel_seg);
        end_seg = end_tracks(sel_seg);
        
        % Check the distance between any two tracks detected withing
        [X,Y] = meshgrid(x_part, y_part);
        dist_matrix = hypot(X-X',Y-Y');
        % Removes the possibility that a track is matched with itself
        temp = ones(length(dist_matrix), 1);
        dist_matrix = dist_matrix + diag(temp)*distance;
        % Creates a logic matrix showing all the tracks that on average are closer
        % than the set threshold. The timing is considered later
        short_dist = dist_matrix<distance;
        
        % Checks which tracks are symultaneously present in the video. The
        % beginning of one is before the end of the other
        start_matrix = meshgrid(start_seg);
        end_matrix = meshgrid(end_seg);
        start_diff = start_matrix' - start_matrix;
        late_start = start_diff>=0;
        end_diff = start_matrix' - end_matrix;
        before_end = end_diff<=0;
        
        % Creates a matrix of all the track pairs that fullfil the requirements and
        % need to be merged
        remove_matrix = short_dist & late_start & before_end;
        
        % Finds all the indeces or the tracks to remove or connect.
        j = 0;
        ind_seg = [];
        for i=1:length(remove_matrix)
            current = remove_matrix(i,:);
            if sum(current)>0
                j = j+1;
                ind_seg(j,:) = [index_seg(i), index_seg(find(current,1))];
            end
        end
        
        % appends to the tracks found in the other areas
        ind_del = [ind_del; ind_seg];
        
    end
end

if ~isempty(ind_del)
    % Flips the links so the track with the lower index is always first
    for i=1:size(ind_del,1)
        if ind_del(i,1)>ind_del(i,2)
            ind_del(i,:) = fliplr(ind_del(i,:));
        end
    end
    
    % Sorts the indeces so the links are ordered starting with the one with the
    % first track with the lowest index first
    ind_del = sortrows(ind_del);
    duplicate = [];
    
    % Removes possible duplicats due to the tracks appearing at the same time
    for i=2:size(ind_del,1)
        if sum(ind_del(i,:) == ind_del(i-1,:)) == 2
            duplicate = [duplicate, i];
        end
    end
    ind_del(duplicate,:) = [];
    
    % Creates a cell array with all the indeces of the tracks to join
    k = 1;
    while true
        if k > size(ind_del, 1)
            break
        end
        cel_del{k} = ind_del(k,:);
        % Check if the same indeces are present in the rest of the ind_del
        temp_el = cel_del{k};
        z=1;
        while z<=length(temp_el)
            % Creates a temporary version of ind_del contaning only the lines
            % after the one considered
            ind_temp = ind_del(k+1:end,:);
            temp = find(ind_temp(:) == temp_el(z));
            if ~isempty(temp)
                    to_delete = [];
                    for ii = 1:length(temp)
                        i = temp(ii);
                        if i>size(ind_temp,1)
                            i = i - size(ind_temp,1);
                            cel_del{k} = [cel_del{k} ind_temp(i,1)'];
                        else
                            cel_del{k} = [cel_del{k} ind_temp(i,2)'];
                        end
                        % Removes the link that has been appended to the
                        to_delete = [to_delete i+k];
                        %ind_del(i+k,:) = [];
                    end
                    ind_del(to_delete, :) = [];
            end
            z = z+1;
        end
        k = k+1;
    end
    
    to_remove = [];
    
    for i=1:length(cel_del)
        current = cel_del{i};
        % Finds the earliest beginning and the latest end
        start_arr = start_tracks(current);
        end_arr = end_tracks(current);
        [~,min_start] = min(start_arr);
        [~,max_end] = max(end_arr);
        % In this case there is not a single spot that is present during all
        % the particle presence and the track has to be stictched together
        if min_start ~= max_end
            % orders tracks by starting time
            [~, temp_sort] = sort(start_arr);
            start_sort = current(temp_sort);
            temp_end = end_tracks(start_sort(1));
            for j=start_sort(2:end)
                if end_tracks(j)>temp_end
                    tracks{start_sort(1)}(tracks{start_sort(1)}(:,6)>0,:) = [];
                    cut_second = tracks{j}(tracks{j}(:,2)>temp_end,:);
                    tracks{start_sort(1)} = [tracks{start_sort(1)}; cut_second];
                    track_props(start_sort(1)) = calculate_properties(tracks{start_sort(1)},start_sort(1));
                    track_props(start_sort(1)).spot_index = unique(tracks{start_sort(1)}(:,7));
                    temp_end = end_tracks(j);
                end
                to_remove = [to_remove j];
            end
        else
            current(min_start) = [];
            to_remove = [to_remove current];
        end
    end
    
    track_props(to_remove) = [];
    tracks(to_remove) = [];
end


% for i=1:length(remove_matrix)
%     current = remove_matrix(i,:);
%     if sum(current)>0
%         first = find(current,1);
%         if track_props(i).ending > track_props(first).ending
%             tracks{first}(tracks{first}(:,6)>0,:) = [];
%             cut_second = tracks{i}(tracks{i}(:,2)>track_props(first).ending,:);
%             tracks{first} = [tracks{first}; cut_second];
%             track_props(first) = calculate_properties(tracks{first},first);
%             track_props(first).spot_index = unique(tracks{first}(:,7));
%         end
%         to_remove = [to_remove i];
%     end
% end





