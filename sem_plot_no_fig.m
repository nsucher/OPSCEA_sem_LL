% SEM_PLOT_NO_FIG
%           Function that inputs semiology matrix csv (sem_matrix_filename)
%           time duration csv (time_mat_filename) and y limit (ylim)

function [SEMperiod,ll_weight_times,ll_w_t_labels] = sem_plot_no_fig(sem_matrix_filename,time_matrix_filename,vid_period,perdur)

%   Created by Natalia Sucher and Jon Kleen May 10 2022, Updated by NS May
%   26 2022
% note: semiology time series is in sampling frequency of 5 Hz
            %-----------------
            % 1. Load data and find number of columns
            sem_matrix = readtable(sem_matrix_filename); %JK
            
            [~,t_cols] = size(sem_matrix); %check if transposed properly 
            
            % ----------------------------------------------------------------
            % 2. For loop to separate string features 
            y_label_names = {};      
            feature_el_vec = []; 
            %feature_el_vec(1,1:t_cols) = nan;
            %[1,t_cols] = size(feature_el_vec) 

            for i = 1:t_cols                                
                feature_el=sem_matrix.Properties.VariableNames{i}; %JK
                feature_el_vec = [feature_el_vec;sem_matrix.Properties.VariableNames{i}];
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
            
                y_label_names = [(y_label_names); ([full_anat ' ' full_pos ' ' full_mot])];
            
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

%                    time_vid_start = seconds(time_mat{1,1})
%                    time_vid_end = seconds(time_mat{1,2})
%                    time_mat_start = seconds(time_mat{1,3})
%                    time_mat_end = seconds(time_mat{1,4})
        
            dur_vid = seconds(time_mat{1,2}) - seconds(time_mat{1,1}); %vid end - vid start = video duration
            dur_mat = seconds(time_mat{1,4}) - seconds(time_mat{1,3}); %mat end - mat start = symptom matrix duration
            mat_start = seconds(time_mat{1,3}) - seconds(time_mat{1,1}); % mat start - vid start = symptom matrix start in entire video
            mat_end = seconds(time_mat{1,4}) - seconds(time_mat{1,1}); % mat end - vid start = symptom matrix end in entire video
                        

            %EC 133-03 (w3)
            %ONLY IF END OF SEMIOLOGY MATRIX CORRESPONDS WITH END OF VIDEO
            vid_start = mat_start + dur_mat - vid_period(2);
            SEMfirst = mat_start - vid_start;  %symptom start in shortened video
            SEMlast = SEMfirst + dur_mat; %symptom end in shortened video
            SEMperiod = [SEMfirst SEMlast];
%             setGlobalSEMperiod(SEMperiod);
           


%%%%%%%%%%%%%%%%%%%
            % 5. Bin present symptoms (clean mat) and first appearance of
            % symptoms (ll_weight_times)

            [rows,cols] = size(nums_t_mat);
            clean_mat = [];  %initialize empty matrix to collect cols with 1,2, or 3 present
                        
            ll_w_t_labels = {}; %initialize labels for ll weight times
            ll_weight_times = []; %initialize empty matrix to store times to analyze weight difference

            count_num = 0;
            for n = 1:cols
                if ismember(1,nums_t_mat(:,n)) | ismember(2,nums_t_mat(:,n)) | ismember(3,nums_t_mat(:,n)) % ignore 0 and 4 values
                    clean_mat(:,n) = nums_t_mat(:,n); %collect all semiology with 1, 2, or 3 values in col vector
                    count_num = count_num + 1;

                    %add first symptom start time
                    first_auto = find(clean_mat(:,n)==1,1,'first');
                    first_tonic =  find(clean_mat(:,n)==2,1,'first');
                    first_clonic = find(clean_mat(:,n)==3,1,'first');

                    %initialize empty vectors to hold before, after, and start of first symptom
                    %3 page matrix: first page is auto, second page is
                    %tonic, third page is clonic
                    
                    % automatism
%                     ll_weight_times(1,:,1) % first symptom time
%                     ll_weight_times(2,:,1) ; % symptom time - perdur seconds
%                     ll_weight_times(3,:,1) ; % symptom time + perdur seconds

                    % tonic
%                     ll_weight_times(1,:,2) % first symptom time
%                     ll_weight_times(2,:,2) % symptom time - perdur seconds
%                     ll_weight_times(3,:,2) % symptom time + perdur seconds

                    % clonic
%                     ll_weight_times(1,:,3) % first symptom time
%                     ll_weight_times(2,:,3) % symptom time - perdur seconds
%                     ll_weight_times(3,:,3) % symptom time + perdur seconds

                    %find times for future ll weight calculation
                    if any(first_auto)
                        
                        ll_weight_times(1,n,1) =  first_auto / 5  + SEMfirst; % divide by 5 milliseconds and add to sem time start
                        ll_weight_times(2,n,1) = ll_weight_times(1,n,1) - perdur; 
                        ll_weight_times(3,n,1) = ll_weight_times(1,n,1) + perdur;  
                    end
                    if any(first_tonic)
                        ll_weight_times(1,n,2) =  first_tonic / 5  + SEMfirst; % divide by 5 milliseconds and add to sem time start
                        ll_weight_times(2,n,2) = ll_weight_times(1,n,2) - perdur; 
                        ll_weight_times(3,n,2) = ll_weight_times(1,n,2) + perdur;
                    end
                    if any(first_clonic)
                        ll_weight_times(1,n,3) =  first_clonic / 5  + SEMfirst; % divide by 5 milliseconds and add to sem time start
                        ll_weight_times(2,n,3) = ll_weight_times(1,n,3) - perdur; 
                        ll_weight_times(3,n,3) = ll_weight_times(1,n,3) + perdur; 
                    end
                    ll_w_t_labels{1,n} = y_label_names{n};

                else
                    y_label_names{n}=''; %Delete y labels of row vectors without 1,2, or 3 value 
                end
            end
            bin_any = any(clean_mat);


            % -----------------------------------------------------------------
%             % 6. Plot with ImageSC
%             m=clean_mat(:,bin_any~=0)'; 
%             semts=linspace(SEMfirst,SEMlast,rows);
%             imagesc(semts,1:size(m,1),m)
% 
%             colorbar('location','southoutside','Ticks',[.5:1:4.5],'TickLabels',{'No Motion','Automatism','Tonic','Clonic','Out of Video'}); %Place colorbar beneath graph; Hard code tick location by summing up # of values used (0-4 in this case)            
%             cmap = [0.8,0.8,0.8;1,1,0;0,.4,1;.9,0,0;1,1,1]; % 0=Grey, 1=yellow, 2=blue, 3=red, 4=white
%             colormap(gca,cmap);
% 
%             xlabel('Time (seconds)')
%             y_ticks = 1:length(find(bin_any~=0));
%             y_labels = y_label_names(bin_any~=0);
%             yticklabels('manual'); %remove 
%             set(gca,'ytick',y_ticks,'CLimMode', 'manual', 'CLim', [0 5]);
% 
%             for n = 1:length(y_labels) %for loop to set y labels inside bars to left
%                 textborder(SEMfirst,y_ticks(n),y_labels(n),'k','w')
%             end
% 
%             title('Semiology')
%             xline([1.5:1:max(get(gca,'ytick'))],'k-',2.5)
%             iter = 100; %for really scrunched plot
%             %iter = 10; %for skinny plot
%             %iter = 5; %for wide plot
%      
%             yline([min(xlim):iter:max(xlim)],'k-',.01)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Time Stamp
            %1