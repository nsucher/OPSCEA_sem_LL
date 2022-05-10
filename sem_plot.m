% SEM_PLOT 
%           Function that inputs semiology matrix csv (sem_matrix_filename)
%           time duration csv (time_mat_filename) and y limit (ylim)

function sem_plot(sem_matrix_filename,time_matrix_filename,ylim

%   Created by Natalia Sucher and Jon Kleen May 10 2022

            %-----------------
            % 1. Load data and find number of columns
            sem_matrix = readtable(sem_matrix_filename); %JK
            
            [~,t_cols] = size(sem_matrix); %check if transposed properly 
            
            % ----------------------------------------------------------------
            % 2. For loop to separate string features 
            full_name_vec = [];            

            for i = 1:t_cols                                
                feature_el=sem_matrix.Properties.VariableNames{i}; %JK
                                       
                anatomy = feature_el(1:2);
                position = feature_el(3);
                motor = feature_el(4);
            
                % Anatomy
                switch anatomy 
                    case 'lu'; full_anat = 'Left Arm';
                    case 'ru'; full_anat = 'Right Arm';
                    case 'lh'; full_anat = 'Left Hand';
                    case 'rh'; full_anat = 'Right Hand';
                    case 'll'; full_anat = 'Left Leg';
                    case 'rl'; full_anat = 'Right Leg';
                    case 'lf'; full_anat = 'Left Foot';
                    case 'rf'; full_anat = 'Right Foot';
                    case 'ed'; full_anat = 'Eye Deviation';
                    case 'eb'; full_anat = 'Eye Blink';
                    case 'lm'; full_anat = 'Left Mouth';
                    case 'rm'; full_anat = 'Right Mouth';
                    case 'ht'; full_anat = 'Head Turn';
                    case 'tt'; full_anat = 'Torso Turn';
                    case 'vx'; full_anat = 'Voice';
                    case 'gm'; full_anat = 'Gyratory Movement';
                    case 'rx'; full_anat = 'Rocking';
                    case 'bb'; full_anat = 'Bimanual Bipedal Automatism';
                    case 'wx'; full_anat = 'Walking';
                    case 'fx'; full_anat = 'Falling';
                    case 'px'; full_anat = 'Pedaling';
                    case 'ba'; full_anat = 'Behavioral Arrest';
                    case 'cg'; full_anat = 'Chapeau de Gendarme';
                    case 'fe'; full_anat = 'Facial Expression';
                    case 'oa'; full_anat = 'Oral Automatism';
                end
            
                % Position
                switch position 
                    case 'p'; full_pos = 'Proximal';
                    case 'd'; full_pos = 'Distal';
                    case 'l'; full_pos = 'Left';
                    case 'r'; full_pos = 'Right';
                    case 'c'; full_pos = 'Center';
                    case 't'; full_pos = 'Twitch';
                    case 'y'; full_pos = 'Pull'; %y stands for yank so pull doesn't get confused with proximal
                    case 's'; full_pos = 'Superior';
                    case 'i'; full_pos = 'Inferior';
                    case 'f'; full_pos = 'Forward';
                    case 'b'; full_pos = 'Backward';
                    case 'n'; full_pos = 'Nonverbal';
                    case 'v'; full_pos = 'Verbal';
                    case 'x'; full_pos = [];
                end
            
                % Motor
                switch motor
                    case 's'; full_mot = 'Simple';
                    case 'c'; full_mot = 'Complex';
                    case 'n'; full_mot = 'Nonfluent';
                    case 'f'; full_mot = 'Fluent';
                    case 'x'; full_mot = [];
                end
            
                full_name_vec = [string(full_name_vec); string([full_anat ' ' full_pos ' ' full_mot])];
            
            end           
            % -----------------------------------------------------------------
            % 3. Separate numerical elements 
            nums_t_mat=sem_matrix{:,:}; % Separate numerical elements
            % -----------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%
            % 4. Align time series 
            % exact times to the millisecond of video
            time_mat = readtable(time_matrix_filename);
            table2cell(time_mat);
            vid_start = time_mat{1,1};
            vid_end = time_mat{1,2};
            mat_start = time_mat{1,3};
            mat_end = time_mat{1,4};
            %vid_duration = vid_end - vid_start; %result in minutes
            vid_vec = [vid_start vid_end];
            mat_vec = [mat_start mat_end];
% 
            [~,~,~,~,mat_s_min,mat_s_sec] = datevec(mat_start);
            [~,~,~,~,mat_e_min,mat_e_sec] = datevec(mat_end);
            [~,~,~,~,vid_s_min,vid_s_sec] = datevec(vid_start);
            [~,~,~,~,vid_e_min,vid_e_sec] = datevec(vid_end); 

            %find duration of semiology matrix
            if mat_s_min > mat_e_min             % Captures minute time differences such as 8:59 to 9:01, ignores hours
                mat_dur_m = mat_s_min - mat_e_min;
            else
                mat_dur_m = mat_e_min - mat_s_min;
            end

            mat_dur_s = mat_e_sec - mat_s_sec; % can be negative -- make sure!

            
            %find duration of video 
            if vid_s_min > vid_e_min             % Captures minute time differences such as 8:59 to 9:01, ignores hours
                vid_dur_m = vid_s_min - vid_e_min;
            else
                vid_dur_m = vid_e_min - vid_s_min;
            end

            vid_dur_s = vid_e_sec - vid_s_sec; % can be negative -- make sure!

            %convert durations to seconds
            total_mat_s = mat_dur_m * 60 + mat_dur_s;
            total_vid_s = vid_dur_m * 60 + vid_dur_s;
            
            mat_start_sec = (mat_s_min * 60 + mat_s_sec); 
            vid_start_sec = (vid_s_min * 60 + vid_s_sec);

            mat_end_sec = (mat_e_min * 60 + mat_e_sec); 
            vid_end_sec = (vid_e_min * 60 + vid_e_sec);

            SEMfirst = mat_start_sec - vid_start_sec; % To put into S later
            SEMlast = SEMfirst + total_mat_s; % To put into S later
            SEMperiod = [SEMfirst SEMlast];
            setGlobalSEMperiod(SEMperiod);

            SEMperiod = getGlobalSEMperiod; % 1x2 vector of semiology start in relation to video length-- what time the symptoms start in the video, as video start as as 0
            SEMstart = SEMperiod(1);
            SEMend = SEMperiod(2);

%%%%%%%%%%%%%%%%%%%
            % 5. Bin present symptoms (clean mat) and first appearance of
            % symptoms (ll_weight_times)

            y_label_names = full_name_vec;
            clean_mat = [];  %initialize empty matrix to collect cols with 1,2, or 3 present
            [rows,cols] = size(nums_t_mat);

            ll_weight_times = []; %initialize empty matrix to store times to analyze weight difference
            auto_lwt = [];


            x = 2; %desired number of seconds 
            x = x * 5; %multiply by 5 to convert to milliseconds
            
            for n = 1:cols
                if ismember(1,nums_t_mat(:,n)) | ismember(2,nums_t_mat(:,n)) | ismember(3,nums_t_mat(:,n)) % ignore 0 and 4 values
                    clean_mat(:,n) = nums_t_mat(:,n); %collect all semiology with 1, 2, or 3 values in col vector

                    %trying to find symptom time start + sem time start 
                    
                    disp(y_label_names(n)) % for testing the code

                    %add first symptom start time
                    first_auto = find(clean_mat(:,n)==1,1,'first');
                    first_tonic =  find(clean_mat(:,n)==2,1,'first');
                    first_clonic = find(clean_mat(:,n)==3,1,'first');
                    
                    %initialize empty vectors to hold before, after, and start of first symptom
                    %3 page matrix: first page is auto, second page is
                    %tonic, third page is clonic
                    
                    % automatism
%                     ll_weight_times(1,:,1) % first symptom time
%                     ll_weight_times(2,:,1) ; % symptom time - x seconds
%                     ll_weight_times(3,:,1) ; % symptom time + x seconds

                    % tonic
%                     ll_weight_times(1,:,2) % first symptom time
%                     ll_weight_times(2,:,2) % symptom time - x seconds
%                     ll_weight_times(3,:,2) % symptom time + x seconds

                    % clonic
%                     ll_weight_times(1,:,3) % first symptom time
%                     ll_weight_times(2,:,3) % symptom time - x seconds
%                     ll_weight_times(3,:,3) % symptom time + x seconds

                    %find times for future ll weight calculation
                    if any(first_auto)
                        ll_weight_times(1,n,1) =  first_auto / 5  + SEMstart; % divide by 5 milliseconds and add to sem time start
                        %base = first_auto / 5  + SEMstart;
                        %disp(["auto: " num2str(first_auto) num2str(base)])
                        disp(["auto: " num2str(first_auto) num2str(ll_weight_times(1,n,1))])
                        %before = base - x;
                        ll_weight_times(2,n,1) = ll_weight_times(1,n,1) - x; 
                        %after = base + x;
                        ll_weight_times(3,n,1) = ll_weight_times(1,n,1) + x;                        
                    end
                    if any(first_tonic)
                        ll_weight_times(1,n,2) =  first_tonic / 5  + SEMstart; % divide by 5 milliseconds and add to sem time start
                        disp(["tonic: " num2str(first_tonic) num2str(ll_weight_times(1,n,2))])
                        ll_weight_times(2,n,2) = ll_weight_times(1,n,2) - x; 
                        ll_weight_times(3,n,2) = ll_weight_times(1,n,2) + x;
                    end
                    if any(first_clonic)
                        ll_weight_times(1,n,3) =  first_clonic / 5  + SEMstart; % divide by 5 milliseconds and add to sem time start
                        disp(["clonic: " num2str(first_clonic) num2str(ll_weight_times(1,n,3))])
                        ll_weight_times(2,n,3) = ll_weight_times(1,n,3) - x; 
                        ll_weight_times(3,n,3) = ll_weight_times(1,n,3) + x; 
                    end
                else
                    y_label_names(n)=''; %Delete y labels of row vectors without 1,2, or 3 value 
                end
            end
            setGlobal_ll_w_t(ll_weight_times)
            bin_any = any(clean_mat);


            % -----------------------------------------------------------------
            % 6. Plot with ImageSC
            m=clean_mat(:,bin_any~=0)';
            ts=getts(size(clean_mat(:,bin_any~=0)',2),5);
            imagesc(ts,1:size(m,1),m)
            colorbar('location','southoutside','Ticks',[.5:1:4.5],'TickLabels',{'No Motion','Automatism','Tonic','Clonic','Out of Video'}); %Place colorbar beneath graph; Hard code tick location by summing up # of values used (0-4 in this case)            
            cmap = [0.8,0.8,0.8;1,1,0;0,.4,1;.9,0,0;1,1,1]; % 0=Grey, 1=yellow, 2=blue, 3=red, 4=white
            colormap(gca,cmap);
            xlabel('Time (seconds)')
            %set(gca,'ytick',1:length(find(bin_any~=0))
            %for loop to set y labels inside bars to left
            y_ticks = 1:length(find(bin_any~=0));
            y_labels = y_label_names(bin_any~=0);
            set(gca,'ytick',y_ticks,'CLimMode', 'manual', 'CLim', [0 5]);
            for n = 1:length(y_labels)
                %text(0,y_ticks(n),y_labels(n));
                textborder(0,y_ticks(n),y_labels(n),'k','w')
            end
            title('Semiology')
            xline([1.5:1:max(get(gca,'ytick'))],'k-',2.5)
            %iter = 10; %for skinny plot
            iter = 5; %for wide plot
            yline([min(xlim):iter:max(xlim)],'k-',.01)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Time Stamp
            % 7. Time Marking
            hold on; 
            plot([ts(i) ts(i)],ylim,'r-')            