load fall_data_unlabeled
load falllabels_101816

num_falls = length(falllabels.types);

% initialize arrays to keep falls data
acce = cell(num_falls,1);      %sensor data for each window
gyro = cell(num_falls,1);
baro = cell(num_falls,1);

% form 10 sec clips
curr_data_ind = 1;
for i=1:length(falllabels.types)
    % current labels data
    fall_start = falllabels.start_true(i);
    fall_end = falllabels.start_end_marked(i,2);

    clip_found = 0;
    for j=curr_data_ind:length(labels.winsize)

        if ~clip_found
            % current sensors data
            %curr_clip_start = labels.evalstart(j);        
            %curr_clip_end = curr_clip_start + labels.winsize(j)*1000;
            curr_clip_start = labels.timestampSTART_END(j,1);        
            curr_clip_end = labels.timestampSTART_END(j,2);
            
            % fall starts in this clip
            if fall_start > curr_clip_start && fall_start < curr_clip_end
                clip_found = 1;
                curr_data_ind = j+1;
                
                acce_10sec = [];
                gyro_10sec = [];
                baro_10sec = [];
                
                % get 10 seconds around the fall_start event
                new_clip_start = fall_start - 4000;
                new_clip_end = fall_start + 6000;
                
                % current clip data
                acce_curr = labels.acce{j};
                gyro_curr = labels.gyro{j};
                baro_curr = labels.baro{j};
                    
                % convert timestamps to absolute
                acce_curr(:,1) = labels.timestampSTART_END(j,1) + 1000*acce_curr(:,1);
                gyro_curr(:,1) = labels.timestampSTART_END(j,1) + 1000*gyro_curr(:,1);
                baro_curr(:,1) = labels.timestampSTART_END(j,1) + 1000*baro_curr(:,1);
                
                
                if new_clip_start < curr_clip_start % take data from previous clip
                    % previous clip data
                    acce_prev = labels.acce{j-1};
                    gyro_prev = labels.gyro{j-1};
                    baro_prev = labels.baro{j-1};
                    
                    % convert timestamps to absolute
                    acce_prev(:,1) = labels.timestampSTART_END(j-1,1) + 1000*acce_prev(:,1);
                    gyro_prev(:,1) = labels.timestampSTART_END(j-1,1) + 1000*gyro_prev(:,1);
                    baro_prev(:,1) = labels.timestampSTART_END(j-1,1) + 1000*baro_prev(:,1);
                    
                    start_ind_acce = find(acce_prev(:,1) > new_clip_start, 1);
                    start_ind_gyro = find(gyro_prev(:,1) > new_clip_start, 1);
                    start_ind_baro = find(baro_prev(:,1) > new_clip_start, 1);
                    
                    acce_10sec = [acce_10sec; acce_prev(start_ind_acce:end,:)];
                    gyro_10sec = [gyro_10sec; gyro_prev(start_ind_gyro:end,:)];
                    baro_10sec = [baro_10sec; baro_prev(start_ind_baro:end,:)];
                    
                    % start for the current clip
                    start_ind_acce = 1;
                    start_ind_gyro = 1;
                    start_ind_baro = 1;
                else
                    start_ind_acce = find(acce_curr(:,1) > new_clip_start, 1);
                    start_ind_gyro = find(gyro_curr(:,1) > new_clip_start, 1);
                    start_ind_baro = find(baro_curr(:,1) > new_clip_start, 1);
                end
                
                % concatenate current clip data
                acce_10sec = [acce_10sec; acce_curr(start_ind_acce:end,:)];
                gyro_10sec = [gyro_10sec; gyro_curr(start_ind_gyro:end,:)];
                baro_10sec = [baro_10sec; baro_curr(start_ind_baro:end,:)];
                
                % next clip data
                acce_next = labels.acce{j+1};
                gyro_next = labels.gyro{j+1};
                baro_next = labels.baro{j+1};
                    
                % convert timestamps to absolute
                acce_next(:,1) = labels.timestampSTART_END(j+1,1) + 1000*acce_next(:,1);
                gyro_next(:,1) = labels.timestampSTART_END(j+1,1) + 1000*gyro_next(:,1);
                baro_next(:,1) = labels.timestampSTART_END(j+1,1) + 1000*baro_next(:,1);                
                
                if new_clip_end < labels.timestampSTART_END(j+1,2)
                    take_one_more_clip = 0;
                    
                    end_ind_acce = find(acce_next(:,1) > new_clip_end, 1);
                    end_ind_gyro = find(gyro_next(:,1) > new_clip_end, 1);
                    end_ind_baro = find(baro_next(:,1) > new_clip_end, 1);
                    
                    acce_10sec = [acce_10sec; acce_next(1:end_ind_acce,:)];
                    gyro_10sec = [gyro_10sec; gyro_next(1:end_ind_gyro,:)];
                    baro_10sec = [baro_10sec; baro_next(1:end_ind_baro,:)];
                else
                    take_one_more_clip = 1;
                    
                    acce_10sec = [acce_10sec; acce_next];
                    gyro_10sec = [gyro_10sec; gyro_next];
                    baro_10sec = [baro_10sec; baro_next];
                end
                
                if take_one_more_clip
                    % one more clip data
                    acce_onemore = labels.acce{j+2};
                    gyro_onemore = labels.gyro{j+2};
                    baro_onemore = labels.baro{j+2};
                    
                    % convert timestamps to absolute
                    acce_onemore(:,1) = labels.timestampSTART_END(j+2,1) + 1000*acce_onemore(:,1);
                    gyro_onemore(:,1) = labels.timestampSTART_END(j+2,1) + 1000*gyro_onemore(:,1);
                    baro_onemore(:,1) = labels.timestampSTART_END(j+2,1) + 1000*baro_onemore(:,1);
                    
                    end_ind_acce = find(acce_onemore(:,1) > new_clip_end, 1);
                    end_ind_gyro = find(gyro_onemore(:,1) > new_clip_end, 1);
                    end_ind_baro = find(baro_onemore(:,1) > new_clip_end, 1);
                    
                    acce_10sec = [acce_10sec; acce_onemore(1:end_ind_acce,:)];
                    gyro_10sec = [gyro_10sec; gyro_onemore(1:end_ind_gyro,:)];
                    baro_10sec = [baro_10sec; baro_onemore(1:end_ind_baro,:)];
                end
                
                acce{i} = acce_10sec;
                gyro{i} = gyro_10sec;
                baro{i} = baro_10sec;

                if fall_end > curr_clip_end
                    curr_data_ind = j+2;
                end
                
                clip = acce_10sec;
                t = clip(:,1);
                clip_duration = clip(end,1)-clip(1,1);
                plot_title = strcat(falllabels.types{i},{' '},falllabels.location{i});
                figure, plot(t,clip(:,2), t,clip(:,3), t,clip(:,4)), legend('X','Y','Z')
                title(plot_title), xlabel(clip_duration)
                y1=get(gca,'ylim'); hold on, plot([fall_start fall_start],y1)
                y1=get(gca,'ylim'); hold on, plot([fall_end fall_end],y1)

            end
            
        else
            % if clip found, move to the next label
            break;
        end
        
    end
end

%put all data into one structure
data.acce = acce;      %sensor data for each window
data.gyro = gyro;
data.baro = baro;

data.type_str = falllabels.types;
data.subject = falllabels.subject;
data.location = falllabels.location;

save falls_data_10sec data