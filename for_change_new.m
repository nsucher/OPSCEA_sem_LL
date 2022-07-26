function for_change_new

% 1. set path
% 2. initialize variables
% 3. drop down menu

%% 1. set path 

loop_path=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
loop_datapath=[loop_path 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
loop_filepath=[loop_datapath 'avg_change_folders'];

cd(loop_filepath)

%% 2. initialize variables

files = ls(loop_filepath);

let = lettersPattern;
dig = digitsPattern;
pat = let + dig; 
new_pat = let + dig + '_' + digitsPattern(2);
file_list = extract(files,pat);

cell_list = [];

avg_fig = 1;

anat_count = 0;
sx_count = 0;

weight_list = {};

symptom_list = {'Left Arm Proximal Simple',...
    'Left Arm Proximal Complex',...
    'Left Arm Distal Simple',...
    'Left Arm Distal Complex',...
    'Right Arm Proximal Simple',...
    'Right Arm Complex',...
    'Right Arm Distal Simple',...
    'Right Arm Complex',...
    'Left Hand Proximal Simple',...
    'Left Hand Proximal Complex',...
    'Left Hand Distal Simple',...
    'Left Hand Distal Complex',...
    'Right Hand	Proximal Simple',...
    'Right Hand	Proximal Complex',...
    'Right Hand	Distal Simple',...
    'Right Hand	Distal Complex',...
    'Left Leg Proximal Simple',...
    'Left Leg Proximal Complex',...
    'Left Leg Distal Simple',...
    'Left Leg Distal Complex',...
    'Right Leg Proximal Simple',...
    'Right Leg Proximal Complex',...
    'Right Leg Distal Simple',...
    'Right Leg Distal Complex',...
    'Left Foot Proximal Simple',...
    'Left Foot Proximal Complex',...
    'Left Foot Distal Simple',...
    'Left Foot Distal Complex',...
    'Right Foot	Proximal Simple',...
    'Right Foot	Proximal Complex',...
    'Right Foot	Distal Simple',...
    'Right Foot	Distal Complex',...
    'Eye Deviation Left',...
    'Eye Deviation Right',...
    'Eye Deviation Center',...
    'Left Mouth	Twitch',...
    'Left Mouth	Pull',...
    'Right Mouth Twitch',...
    'Right Mouth Pull',...
    'Head Turn Left',...
    'Head Turn Right',...
    'Head Turn Superior',...
    'Head Turn Inferior',...
    'Head Turn Forward',...
    'Head Turn Backward',...
    'Torso Turn	Left',...
    'Torso Turn	Right',...
    'Torso Turn	Forward',...
    'Torso Turn	Backward',...
    'Voice Nonverbal Simple',...
    'Voice Nonverbal Complex',...
    'Voice Verbal Nonfluent',...
    'Voice Verbal Fluent',...
    'Gyratory Movement',...
    'Rocking',...
    'Bimanual Bipedal Automatisms',...		
    'Walking',...
    'Falling',...
    'Pedaling',...
    'Behavioral Arrest',...	
    'Chapeau de Gendarme',...	
    'Facial Expression',...
    'Oral Automatism'};

neuroanat_list = {'frontalpole',...% FRONTAL LOBE
    'parstriangularis',...
    'parsopercularis',...
    'parsorbitalis',...
    'rostralmiddlefrontal',...
    'caudalmiddlefrontal',...
    'lateralorbitofrontal',...
    'superiorfrontal',...
    'medialorbitofrontal',...
    'precentral',...
    'postcentral',... % PARIETAL LOBE
    'inferiorparietal',...   
    'superiorparietal',...
    'supramarginal',...
    'temporalpole',...% TEMPORAL LOBE
    'middletemporal',...
    'superiortemporal',...
    'inferiortemporal',...
    'parahippocampal',...               
    'Right-Hippocampus',...
    'Left-Hippocampus',...
    'Right-Amygdala',...
    'Left-Amygdala',...
    'entorhinal',...
    'bankssts',...
    'fusiform',...% OCCIPITAL LOBE                
    'lingual',...
    'Right-Inf-Lat-Vent',...% OTHER
    'Right-Cerebral-White-Matter',...
    'Left-Cerebral-White-Matter',...
    'Right-choroid-plexus',...
    'Right-Putamen',...
    'Right-VentralDC'};
[~,n_col] = size(neuroanat_list);

%% 3. Drop Down Menu

set(0,'DefaultFigureVisible','on')

fig = uifigure('Position',[100 300 300 500]);

        %symptom drop down
        sx_list = unique(symptom_list);
        default_sx = sx_list{1};
        symptom = uidropdown(fig,'Position',[84 204 150 20], ...
            'Items',sx_list, ...
            'Value',default_sx,...
            'ValueChangedFcn',@(dd,event) sx_selection(dd,default_sx));
        
        %mode drop down
        mx_list = {'Automatism', 'Tonic', 'Clonic'};
        default_mode = 1;
        mode = uidropdown(fig,'Position',[84 184 150 20],...
            'Items',mx_list,...
            'ItemsData',1:3,...
            'Value',default_mode,...
            'ValueChangedFcn',@(dd,event) mx_selection(dd,default_mode));

        %duration of average drop down 
        %   (# of seconds before and after symptom onset)
        perdur_range = 2:2:20;
        perdur_list = {};
        for i = 1:length(perdur_range)
            perdur_list{i} = num2str(perdur_range(i));
        end


        default_perdur = perdur_list{1};
        perdur_dd = uidropdown(fig,'Position',[84 164 150 20],...
            'Items',perdur_list,...
            'Value',default_perdur,...
            'ValueChangedFcn',@(dd,event) mx_selection(dd,default_perdur));

        disp('To submit values press return: ')
        pause
        
        chosen_mode = mode.Value
        chosen_symptom = symptom.Value
        chosen_perdur = str2num(perdur_dd.Value)

%% 4. Collect symptoms of all seizures  
for i = 1:length(file_list)
    each_filepath = [loop_filepath '/' file_list{i}];
    cd(each_filepath)
    each_file = ls(each_filepath);
    cell_list = [cell_list extract(each_file,new_pat)];
    pt_sz = split(cell_list{i},'_');
    pt = pt_sz{1};
    sz = pt_sz{2};

    ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
    ptpath=[loop_filepath '/' pt '/']; % patient's folder
   
    opsceapath=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
    opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
        if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
    cd(opsceapath);

    % Import parameters for specific seizure 
    [~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
        fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
        prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
        if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end
        cd 

    % Get time segments within the ICEEG file to use
    VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
    S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];

    vid_period = S.VIDperiod; %NS
    
    load([ptpath ptsz '/' ptsz '.mat'])

    
    mat_name = [study '_mat.csv'];   %string for study_mat.csv
    mat_time_name = [study '_time_mat.csv']; %string for study_time_mat.csv


    [~,~,ll_w_t_labels] = sem_plot_no_fig(mat_name,mat_time_name,vid_period,chosen_perdur); %insert filename of semiology matrix as first parameter       
% 
%     u_labels = unique(ll_w_t_labels);


    [~,sx_col] = size(ll_w_t_labels);


    %collect labels from each file into a cell string
    for j = 1:sx_col
        if ~isnan(ll_w_t_labels{1,j})
           sx_count = sx_count + 1;
           all_label_1d{i,sx_count} = ll_w_t_labels{1,j};
        else
            continue
        end
    end
end

%% 5. Match symptoms to drop down selection
chosen_symptom = cellstr(chosen_symptom);
sx_exist = strcmpi(chosen_symptom,all_label_1d); %find if neuroanatomy exists in each seizure

if sum(sx_exist) < 1
    disp('ERROR: No matching symptoms') %most likely due to space in front of full_mot 
else
end

[l_row,~] = size(all_label_1d);

sz_to_plot = zeros(1,l_row);
for l = 1:l_row
    if sum(sx_exist(l,:)) >= 1 
        sz_to_plot(l) = 1;
    else
        sz_to_plot(l) = 0;
    end
end

n_label = {};
sz_count = 0;
n_count = 0;



figure


%% 6. Collect neuroanatomy labels and weight change from each seizure
for sz_i = sz_to_plot % for each seizure with symptom present
    sz_count = sz_count + 1;
    % load seizure-specific neuroanatomy labels
    if sz_i > 0
        each_filepath = [loop_filepath '/' file_list{sz_count}];
        cd(each_filepath)
        each_file = ls(each_filepath);
        cell_list = [cell_list extract(each_file,new_pat)];
        pt_sz = split(cell_list{sz_count},'_');
        pt = pt_sz{1};
        sz = pt_sz{2};
    
        % collect electrode matrix, weights, difference, unique anat (AKA present neuroanatomy), pval
        [clean_ll_w_t_labels,anatstructureselec_weights,EM,ua,pval] = avg_change(pt,sz,avg_fig,chosen_symptom,chosen_mode,chosen_perdur);
        
        weight_count = 0;
        
        %collect outputs from each seizure into array
        [an_row,~] = size(anatstructureselec_weights);
        [ua_row,~] = size(ua);
        [pv_row,~] = size(pval);
        [~,em_col] = size(EM);
        [~,cl_col] = size(clean_ll_w_t_labels);

        for an = 1:an_row
            anatstructureselec_weights_2d{an,sz_count} = anatstructureselec_weights{an};
        end

        for ua_i = 1:ua_row
            ua_2d{ua_i,sz_count} = ua{ua_i};            
        end



        % for each neuroanat_item
        for neuroanat_item = 1:length(neuroanat_list)
            % check to see if neuroanatomy is present in seizure  
            n_exist = strcmpi(neuroanat_list{neuroanat_item},ua);
            if sum(n_exist) > 0 
                n_label{neuroanat_item,sz_count} = neuroanat_list{neuroanat_item};

%                 weight_count = weight_count + 1;

                %put weights in cell array weight list
%                 weight_list{neuroanat_item,sz_count} = anatstructureselec_weights{weight_count};
                weight_list{neuroanat_item,sz_count} = anatstructureselec_weights{n_exist};

            else
                n_label{neuroanat_item,sz_count} = nan;
                weight_list{neuroanat_item,sz_count} = nan;
            end
        end
    end
end

%% 7. Plot
A = cellfun(@ischar,n_label,'UniformOutput',true);
B = sum(A,2)>0; %index of present neuoranatomy to use in plotting data -- follow C format
plot_neuroanat = neuroanat_list(B);
for f = 1:length(file_list)
    plot_weights(:,f) = weight_list(B,f);
end

[plot_row,plot_col] = size(plot_weights);

figure(1)

for sp = 1:length(file_list)
    subplot(1,length(file_list),sp)
    for p_r = 1:plot_row
        plot(plot_weights{p_r,sp},p_r*ones(length(plot_weights{p_r,sp}),1),'ko')
        xlim(max(abs(xlim))*[-1 1]) 
        ylim([0 length(plot_neuroanat)])
        yline(0,'k-'); 
        xline(p_r,'G:',.25);
        set(gca,'ytick',1:plot_row)
    end
    if sp == 1
        subplot(1,length(file_list),1)
        set(gca,'YTickLabel',plot_neuroanat,'fontsize',14)
    else
        set(gca,'YTickLabel',[],'fontsize',14)
    end
    cell_name = join(split(cell_list{sp},'_'),'-');
    xlabel(cell_name)
%     hold on;
end


% 
1
1
