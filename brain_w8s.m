
%% get the weights
function brain_w8s(pt,sz,clean_ll_w_t_l,clean_LL_diff,new_anat,noneed,em)
% the values in this variable are what we need -->  LL_meandiff %find clean version
% for example, to plot tonic head turn, get the weights

clean_diff_3d = clean_LL_diff;
list_symp = clean_ll_w_t_l;

% LL_q8 = readmatrix('q8_LL.csv');
% clean_diff_2d = readmatrix('q8_clean_LL_diff.csv');
% list_symp = readcell('q8_clean_ll_w_t_l.csv');
% 
% [row_2d,col_2d] = size(clean_diff_2d);

%convert clean_diff_2d to 3d matrix
% clean_diff_3d = nan(row_2d,(col_2d/3),3);
% 
% %page 1
% clean_diff_3d(:,:,1) = clean_diff_2d(1:row_2d,1:(col_2d/3)); %w8s is a vector of all electrodes
% %page 2
% clean_diff_3d(:,:,2) = clean_diff_2d(1:row_2d,((col_2d/3)+1):((col_2d/3)*2)); %w8s is a vector of all electrodes
% %page 3
% clean_diff_3d(:,:,3) = clean_diff_2d(1:row_2d,((col_2d/3)*2+1):col_2d); %w8s is a vector of all electrodes



%ask for user input (what mode to plot)
%modee=2;


%ask for user input (what symptom to plot)
  % *will need to list only the ones that are not nans
disp('List of symptoms for this seizure: ')
disp(list_symp)
prompt = 'Column number of symptom to plot: ';
col_num = input(prompt);

disp('List of modes for symptoms: ')
list_mode = ['1. automatism ','2. tonic ','3. clonic'];
disp(list_mode)
prompt2 = 'Column number of mode to plot: ';
col_num2 = input(prompt2);




% in case any of these electrodes have nan values, make nns to index the non-nan electrodes
w8s = squeeze(clean_diff_3d(:,col_num,col_num2));


nns=~isnan(w8s);



%% get the meshes for both hemisphered ready

opsceapath=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
    if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
cd(opsceapath);


%%HARDCODED Path
% pt = 'EC96';
% sz = '01';

ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
ptpath=[opsceadatapath pt '/']; % patient's folder
szpath= [ptpath ptsz '/']; % specific seizure's folder
disp(['Running ' pt ', seizure ' sz '...']);
cd(szpath) % Save video in the same data folder for that seizure

meshpath='Imaging/Meshes/';
Rcortex=load([ptpath meshpath pt '_rh_pial.mat']); 
    loaf.rpial=Rcortex; 
    Rcrtx=Rcortex.cortex; 
    clear Rcortex
Lcortex=load([ptpath meshpath pt '_lh_pial.mat']); 
    loaf.lpial=Lcortex; 
    Lcrtx=Lcortex.cortex; 
    clear Lcortex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
    fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
    prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
    if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end
% Import parameters for patient's specific plot (layout of video frame)
[~,plt]=xlsread([opsceapath 'OPSCEAparams'],pt); 
    fields_PLOT=plt(1,:); plt(1,:)=[]; % header for columns of plotting parameters
    plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)

cd 
% prepare subplot specifications
%     subplotrow=str2double(plt(:,strcmpi(fields_PLOT,'subplotrow')));
%     subplotcolumn=str2double(plt(:,strcmpi(fields_PLOT,'subplotcolumn')));
     subplotstart=plt(:,strcmpi(fields_PLOT,'subplotstart')); subplotstop=plt(:,strcmpi(fields_PLOT,'subplotstop')); 
     for j=1:length(plottype); subplotnum{j,1}=str2double(subplotstart{j}):str2double(subplotstop{j});
     end
     surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
     surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
%     viewangle=lower(plt(:,strcmpi(fields_PLOT,'viewangle')));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if exist([ptpath 'Imaging/Elecs/Electrodefile.mat'])
%         load([ptpath 'Imaging/Elecs/Electrodefile.mat']); 
%     elseif exist([ptpath 'Imaging/elecs/clinical_elecs_all.mat']) % access variables in old format NS
%          load([ptpath 'Imaging/elecs/clinical_elecs_all.mat']); % NS
%     end
% 
%     if ~exist('anatomy','var'); 
%         anatomy=cell(size(elecmatrix,1),4); 
%     end
%     if size(anatomy,1)>size(elecmatrix,1); 
%         anatomy(size(elecmatrix,1)+1:end,:)=[]; 
%     end
%     anat=anatomy; 
%     clear anatomy; 
%     if size(anat,2)>size(anat,1); 
%         anat=anat'; 
%     end
%     if size(anat,2)==1; 
%         anat(:,2)=anat(:,1); 
%     end; 
%     if ~exist('eleclabels','var'); 
%         eleclabels=anat(:,1); 
%     end
% 
% %% load ICEEG data, and the bad channels verified for that specific data
% load([szpath ptsz])
% load([szpath ptsz '_badch']); 
% % rename and clear old format of electrode files -NS
% if exist('bad_chs','var')
%     badch = bad_chs; clear bad_chs; end; 
% 
%    em=elecmatrix; 
%    clear elecmatrix;
%    emnan=isnan(mean(em,2)); 
%    badch(emnan)=1; 
%    em(emnan,:)=0; 


%PULL IN NO NEED AND MATCH EM(NNS) WITH W8S(NNS)

%define anat for noneed
%     noneed=false(size(anat,1),4);
%     anat_rows = size(anat,1); 
%     badch = badch(1:anat_rows); % cut down badch to eliminate extra bad channels after anatomy rows size
%     for num_rows=1:size(anat,1)
%       noneed(num_rows,1) = contains(lower(anat{num_rows,4}),'ctx'); % cell array containing row index of strings with "ctx" in u1
%       noneed(num_rows,2) = contains(lower(anat{num_rows,4}),'wm'); % cell array containing row index of strings with "wm" in u1
%       noneed(num_rows,3) = contains(lower(anat{num_rows,4}),'white-matter'); % cell array containing row index of strings with "Right-Cerebral-White-Matter" in u1             
%       noneed(num_rows,4) = contains(lower(anat{num_rows,4}),'unknown'); % cell array containing row index of strings with "Unknown" in u1
%       noneed(num_rows,5) = contains(lower(anat{num_rows,4}),'vent'); % cell array containing row index of strings with "Unknown" in u1             
%     end
%     noneed=any(noneed,2);
%     
%     noneed = noneed | badch; % now is either uncessary (noneed) or bad channels
%     
%     new_LL=LL_q8;
%     new_anat=anat;
%     new_LL(noneed,:)=[];

   em=em(~noneed,:);
   emnan=isnan(mean(em,2)); 
   badch(emnan)=1; 
   em(emnan,:)=0; 


%%%%%% INITIALIZE VARIABLES FOR CTMR_GAUSS
S.iceeg_scale=prm{strcmp('iceeg_scale',fields_SZ)}; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
    if ischar(S.iceeg_scale); S.iceeg_scale=str2double(S.iceeg_scale); end 

S.fps=str2double(prm{strcmp('fps',fields_SZ)});             %frames per sec of ICEEG (default 15)

S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap

S.gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
    params={'iceeg_scale','fps','cax','gsp'}; 
    paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]); 
    if any(paramsnans); error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']); 
    end

    %%%%%% 
for j=1:size(plt,1) 
          srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
          srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
          if length(srf)~=length(srfalpha); msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet'); return; end
          acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
            for s=1:length(srf)
                srf{s}=lower(srf{s}); %convert to lower case for easier string matching
              if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                switch srf{s}; %see below for case "wholebrain"s
                    case 'rcortex'; srfplot=Rcrtx; 
                    case 'lcortex'; srfplot=Lcrtx; 
%                     case 'rhipp';   srfplot=Rhipp; 
%                     case 'lhipp';   srfplot=Lhipp; 
%                     case 'ramyg';   srfplot=Ramyg; 
%                     case 'lamyg';   srfplot=Lamyg; 
                end
              end
            end
end



%% plot the weights on the brain
hh=ctmr_gauss_plot_edited(srfplot,em(nns,:),w8s(nns),S.cax,0,cmocean('balance'),S.gsp); 

%     alpha(hh,str2double(srfalpha{s})); % Adjust opacity specified for that row
lightsout
litebrain('r',1)
litebrain('l',1)

alpha(1)




