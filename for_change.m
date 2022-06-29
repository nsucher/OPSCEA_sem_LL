% FOR_CHANGE        Run all selected seizure data in for loop to display
%                   average electrode weight change as heatmap with
%                   distribution
% enter perdur, symptom, mode and then put into avg change

% Natalia Sucher
% Dr. Jon Kleen
% Created:      06/13/2022
% Last Edited:  06/19/2022
% UCSF Neurology, Epilepsy Department
% for_change.m

function for_change

% PSEUDO-CODE
% 1. Set path
% 2. Initialize variables
% 3. Run for loop

% 1. set path for folders to loop through

loop_path=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
loop_datapath=[loop_path 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
loop_filepath=[loop_datapath 'avg_change_folders'];

cd(loop_filepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. initialize variables

files = ls(loop_filepath);

let = lettersPattern;
dig = digitsPattern;
pat = let + dig; 
new_pat = let + dig + '_' + digitsPattern(2);
file_list = extract(files,pat);


cell_list = [];

avg_fig = 1;
perdur = 10;

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

%NEED TO ORGANIZE BY BRAIN PART

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

sx_count = 0;
anat_count = 0;
sz_count = 0;


EM_2d = {};
clean_ll_w_t_label_2d = {};
anatstructureselec_weights_2d = {};
pval_2d = {};
ua_2d = {};


neuro_i_vec = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. run for loop

for i = 1:length(file_list)
%     set(0,'DefaultFigureVisible','on')
    subplot(1,length(file_list),i)

    each_filepath = [loop_filepath '/' file_list{i}];
    cd(each_filepath)
    each_file = ls(each_filepath);
    cell_list = [cell_list extract(each_file,new_pat)];
    pt_sz = split(cell_list{i},'_');
    pt = pt_sz{1};
    sz = pt_sz{2};

    ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
    ptpath=[loop_filepath '/' pt '/']; % patient's folder

    if exist([ptpath 'Imaging/Elecs/Electrodefile.mat'])
      load([ptpath 'Imaging/Elecs/Electrodefile.mat']); 
    elseif exist([ptpath 'Imaging/elecs/clinical_elecs_all.mat']) % access variables in old format NS
       load([ptpath 'Imaging/elecs/clinical_elecs_all.mat']); % NS
    end;
    if ~exist('anatomy','var'); 
        anatomy=cell(size(elecmatrix,1),4); 
    end
    if size(anatomy,1)>size(elecmatrix,1); 
        anatomy(size(elecmatrix,1)+1:end,:)=[]; 
    end

    anat=anatomy; 
    clear anatomy; 
    if size(anat,2)>size(anat,1); 
        anat=anat'; 
    end
    if size(anat,2)==1; 
        anat(:,2)=anat(:,1); 
    end; 

    new_anat = anat;


    opsceapath=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
    opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
        if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
    cd(opsceapath);

    %% Import parameters
    % for specific seizure 
    [~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
        fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
        prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
        if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end
    
        cd 

 

%% Get time segments within the ICEEG file to use
    VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
    S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NEUROANATOMY

  
%get ua without running all of avg_change
    [u_anatomy,~,~]=unique(new_anat(:,4));

    [anat_row,anat_col] = size(u_anatomy);
    for a = 1:anat_row
        anat_count = anat_count + 1;
        all_anat_1d{anat_count} = u_anatomy{a,1};
    end

%%%%%%%%%%%%%
% SYMPTOMS




    vid_period = S.VIDperiod; %NS


%     set(0,'DefaultFigureVisible','off')
    
    load([ptpath ptsz '/' ptsz '.mat'])

    
    mat_name = [study '_mat.csv'];   %string for study_mat.csv
    mat_time_name = [study '_time_mat.csv']; %string for study_time_mat.csv


    [~,~,ll_w_t_labels] = sem_plot_no_fig(mat_name,mat_time_name,vid_period,perdur); %insert filename of semiology matrix as first parameter       
% 
%     u_labels = unique(ll_w_t_labels);


    [~,sx_col] = size(ll_w_t_labels);



    %collect labels from each file into a cell string
    for j = 1:sx_col
        if ~isnan(ll_w_t_labels{1,j})
           sx_count = sx_count + 1;
           all_label_1d{sx_count} = ll_w_t_labels{1,j};
        else
            continue

        end

%         all_label_3d{1,j,i} = clean_ll_w_t_labels{1,j};
%         all_label_2d{i,j} = clean_ll_w_t_labels{1,j};

    end

    if i == length(file_list)
        u_anat = unique(all_anat_1d);
        
        % overlap = strcmpi(u_all_label,symptom_list);
        set(0,'DefaultFigureVisible','on')
        
        fig = uifigure('Position',[100 300 300 500]);

        %symptom drop down
        sx_list = unique(all_label_1d);
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

%         neuroanat = uidropdown(fig,'Position',[84 164 200 20],'Items',u_anat);
%         chosen_anat = neuroanat.Value;

        for item = 1:length(cell_list)

            pt_sz = split(cell_list{item},'_');
            pt = pt_sz{1};
            sz = pt_sz{2};
            set(0,'DefaultFigureVisible','on')

            [clean_ll_w_t_labels,anatstructureselec_weights,EM,ua,pval] = avg_change(pt,sz,avg_fig,chosen_symptom,chosen_mode,chosen_perdur);


%             EM_list() = {};
            [~,cl_col] = size(clean_ll_w_t_labels);
            [an_row,~] = size(anatstructureselec_weights);
            [~,em_col] = size(EM);
            [ua_row,~] = size(ua);
            [pv_row,~] = size(pval);

%%%%%%%%%%%
% MAKE FOR LOOP STANDARDIZE INDEX FOR Y LABELS
    % WHAT STRUCTURES TO PUT IN ORDER
    % DEFINE IN CODE
    % RUN FOR LOOP
    % GO THROUGH AND FIND ALL ELECTRODE THAT MATCH ANATOMY OF PREFERRED
    % ORDER
    % GIVE INDICES OF 1 THROUGH X
    % NEXT STRUCTURE 
    % X+1 THROUGH Y 

    %DECIDE STRUCTURE (LOOK AT ANAT LIST) 

%%%%%%%%%


            %cols
%             for label = 1:cl_col
%                 clean_ll_w_t_label_2d{label,item} = clean_ll_w_t_labels{1,label};
%             end
% 
%             for e = 1:em_col
%                 EM_2d{e,item} = EM{1,e};
%             end

            %rows
            for an = 1:an_row
                anatstructureselec_weights_2d{an,item} = anatstructureselec_weights{an};
            end


            for ua_i = 1:ua_row
                ua_2d{ua_i,item} = ua{ua_i};
            end

%             for pv = 1:pv_row
%                 pval_2d{pv,item} = pval(pv);
%             end

                    [u_row,~] = size(ua_2d);


            figure(1)
            subplot(1,length(cell_list),item) 
              if isempty(pval) == 1
                  continue
              else
               for u=1:u_row %this is all using SORTED orders so be cautious!
%                       neuro_i = strcmpi(ua_2d{u,item},neuroanat_list);
                      if u <= length(ua)
                          neuro_i = strcmpi(ua{u},neuroanat_list);
                          neuro_i_num = find(neuro_i==1);
                          neuro_i_vec{u} = neuro_i_num;
    %                       neuro_i_vec{u,item} = neuro_i_num;
                          
                          neuro_var = neuroanat_list{neuro_i};
                      else
                          continue
                      end

                      if pval(u)<.05
                             mrkr='r*'; 
                      elseif pval(u)>.05
                             mrkr='ko';
                      else isnan(pval(u))
                            continue
                      end
                      
                     plot(anatstructureselec_weights_2d{u,item},u*ones(length(anatstructureselec_weights_2d{u,item}),1),mrkr)
                     xlim([max(abs(xlim))*[-1 1]]); 
                     yline(0,'k-'); 
                     for u=1:u_row;
                         xline(u,'G:',.25); 
                     end
                     
                     set(gca,'ytick',1:u_row)
                     if item == 1
                         set(gca,'YTickLabel',ua,'fontsize',18)
                     end
                     hold on;
                     end
                  end

                end
             end
end



1
1