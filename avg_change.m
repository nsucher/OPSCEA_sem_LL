% AVG_CHANGE            Plot average change of electrode weights as a heatmap projected onto a 3D 
%                       reconstructed brain x seconds before and after symptom onset                      
% Input: 
%       pt: patient filename
%       sz: seizure folder
%       fig: display figure (1) or not (0) of brain weights
%       sx: symptom
%       mx: mode
%       perdur: period of duration/number of seconds before and after symptom onset        

% Output: 
%       anatstructureselec_weights:
%       EM:
%       ua: 
%
% Example:
%       avg_change('EC107','01',1,'Head Turn Left',2,10)

% Natalia Sucher
% Dr. Jon Kleen
% Created:      06/08/2022
% Last Edited:  06/13/2022
% UCSF Neurology, Epilepsy Department
% avg_change.m

function [clean_ll_w_t_l,anatstructureselec_weights,EM,ua,pval] = avg_change(pt,sz,fig,sx,mx,perdur)


% PSEUDOCODE
% 1. Set path 
% 2. Initialize global variables 
% 3. Initialize local variables
% 4. Get variables from sem & LL functions
% 5. Plot average changes 
% 6. Graph distribution



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. SET PATH FROM PARAMS


opsceapath=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
    if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
cd(opsceapath);

ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
ptpath=[opsceadatapath pt '/']; % patient's folder
szpath= [ptpath ptsz '/']; % specific seizure's folder
disp(['Running ' pt ', seizure ' sz '...']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Initialize global variables 

%% Initiate global variables
  global S; % holds general parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. Initialize local variables


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
    vid_period = S.VIDperiod; %NS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BLstart=prm(:,strcmpi(fields_SZ,'BLstart')); BLstop=prm(:,strcmpi(fields_SZ,'BLstop')); %chunk of data (seconds into ICEEG data file) to use for baseline (for z-score step)
    S.BLperiod=[str2double(BLstart{1}) str2double(BLstop{1})];

%transform, scaling, and display options
S.llw=str2double(prm{strcmp('llw',fields_SZ)}); %default linelength window (in seconds)
S.iceeg_scale=prm{strcmp('iceeg_scale',fields_SZ)}; %percentile (number >50 and <100), used here similar to gain ICEEG waveform display, usually 95
     if ischar(S.iceeg_scale); S.iceeg_scale=str2double(S.iceeg_scale); end 
S.fps=str2double(prm{strcmp('fps',fields_SZ)});             %frames per sec of ICEEG (default 15)
S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap
S.gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
    params={'iceeg_scale','fps','cax','gsp'}; 
    paramsnans=isnan([(isnan(S.iceeg_scale) | S.iceeg_scale<=50 | S.iceeg_scale>=100)   S.fps   any(isnan(S.cax)) S.gsp]); 
    if any(paramsnans); error(['ATTENTION OPSCEA USER: The "' params{paramsnans} '" term(s) is/are in an incorrect format (perhaps number instead of string), check excel seizure parameter sheet']); 
    end 
  cm=prm{strcmp('cm',fields_SZ)};
  switch cm; case 'cmOPSCEAcool'; cm=cmOPSCEAcool; 
             case 'cmOPSCEAjet'; cm=cmOPSCEAjet; 
  end
S.cm=cm; %colormap to use for heatmap
S.iceegwin=str2double(prm{strcmp('iceegwin',fields_SZ)}); %how much trace-based ICEEG to view at a time in the ICEEG window
S.marg=str2double(prm{strcmp('marg',fields_SZ)}); %offset of real-time LL txform from beginning of viewing window (in sec; converts to samples below)
% S.slicebright=str2double(prm{strcmp('slicebright',fields_SZ)}); if isnan(S.slicebright); S.slicebright=0; end %brighten up slices (usually 0 to 50)


% additional adjustment for display window
S.VIDperiod=[S.VIDperiod(1)-S.marg   S.VIDperiod(2)+S.iceegwin-S.marg]; 
S.fields=fields_SZ; clear fields

S.prm=prm; clear prm
S.prm_allPtSz=prm_allPtSz; clear prm_allPtSz
S=orderfields(S); %alphabetize the structure fields for ease of use/search

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load ICEEG data, and the bad channels verified for that specific data
load([szpath ptsz])
load([szpath ptsz '_badch']); 
% rename and clear old format of electrode files -NS
if exist('ppEEG','var')
    d = ppEEG; clear ppEEG; end;
if exist('fs','var')
    sfx = fs; clear fs; end;
if exist('bad_chs','var')
    badch = bad_chs; clear bad_chs; end; 
% end of edit -NS
if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
[~,ntp]=size(d); f=1;  %EDITED JK 5/4
% disp(['Length of data to play for video is ' num2str(round(ntp/sfx)) 'sec'])

% error checks for selected time periods
if any([S.VIDperiod(1) S.BLperiod(1)]<0)
    error('VIDperiod is out of bounds of file (time < 0). Check both VIDstart and BLstart times and make sure the "marg" value (subtracted from VIDstart and BLstart), cannot be < 0'); 
elseif any([S.VIDperiod(2) S.BLperiod(2)]>ntp)
    error('VIDperiod is beyond the length of the file. Check VIDstop and BLstop times vs. length of actual ICEEG data'); 
end 

%% locate and load electrode file for labels and XYZ coordinates
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
    if ~exist('eleclabels','var'); 
        eleclabels=anat(:,1); 
    end
   em=elecmatrix; 
   clear elecmatrix; 
   emnan=isnan(mean(em,2)); 
   em(emnan,:)=0; 
   EKGorREF=strcmpi('EKG1',anat(:,1))|strcmpi('EKG2',anat(:,1))|strcmpi('EKG',anat(:,2))|strcmpi('EKGL',anat(:,2))|strcmpi('REF',anat(:,1)); anat(EKGorREF,:)=[]; em(EKGorREF,:)=[]; eleclabels(EKGorREF,:)=[]; 
   d(size(anat,1)+1:end,:)=[];
   badch(size(anat,1)+1:end)=[];
   [nch,~]=size(d); %%%%%%%% edited JK 5/4

isR=nansum(em(:,1))>0; isL=isR~=1; %handy binary indicators for laterality


%% formatting checks, and consolidation of bad channels
ns=unique( [find(badch);   find(isnan(mean(em,2)));   find(isnan(mean(d,2)))]  ); % bad channels: those that are pre-marked, or if NaNs in their coordinates or ICEEG data traces
nns=true(nch,1); nns(ns)=0; %nns=find(nns); %consolidate bad channels and those with NaNs
%remove data channels previously clipped at end. Only include that which has electrode coordinates (intracranial)
if size(em,1)>size(d,1); nch=size(em,1); d(nch+1:end,:)=[]; LL(nch+1:end,:)=[]; nns(nch+1:end)=[]; ns(ns>size(em,1))=[]; 
   fprintf(2, 'ALERT: Clipping off extra bad channel entries (make sure you have the right ICEEG and bad channel files loaded)\n');
end


%% ICEEG data processing and transform

% Filter out < 1 Hz (and up to nyquist) out to help decrease movement
% artifact and improve stability during ICEEG trace plotting
[b,a]=butter(2,[1 round(sfx/2-1)]/(sfx/2),'bandpass'); 
d(nns,:) = filtfilt(b,a,d(nns,:)')';

% Line-length transform
L=round(S.llw*sfx)-1; % number of samples to calculate line length
LL=nan(size(d)); 
for i=1:size(d,2)-L; LL(:,i)=sum(abs(diff(d(:,i:i+L),1,2)),2); end

%Normalize LL (to baseline period) as z-scores, "zLL"
BLstartsample=max([sfx*S.BLperiod(1) 1]); BLendsample=sfx*S.BLperiod(2); 
for i=1:nch % z-score using channel-specific baseline
    LL(i,:)=(LL(i,:)-nanmean(LL(i,BLstartsample:BLendsample)))/nanstd(LL(i,BLstartsample:BLendsample)); 
end 

% make a timtestamp vector
ts=0:1/sfx:size(d,2)*(1/sfx)-1/sfx;


% Extract the period of data to be used for the video (remove flanking data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vidperiodidx=round(S.VIDperiod(1)*sfx+1):S.VIDperiod(2)*sfx;

d=d(:,vidperiodidx); ntp=length(vidperiodidx);
LL=LL(:,vidperiodidx); 
ts=ts(vidperiodidx); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ytl=eleclabels(nns,1); 
nch=length(find(nns)); 


%Clear non-labeled channels beyond size of number of electrode rows
%Vectorize unwanted channels 
noneed=false(size(anat,1),4);
anat_rows = size(anat,1); 
badch = badch(1:anat_rows); % cut down badch to eliminate extra bad channels after anatomy rows size
for num_rows=1:size(anat,1)
noneed(num_rows,1) = contains(lower(anat{num_rows,4}),'ctx'); % cell array containing row index of strings with "ctx" in u1
noneed(num_rows,2) = contains(lower(anat{num_rows,4}),'wm'); % cell array containing row index of strings with "wm" in u1
noneed(num_rows,3) = contains(lower(anat{num_rows,4}),'white-matter'); % cell array containing row index of strings with "Right-Cerebral-White-Matter" in u1             
noneed(num_rows,4) = contains(lower(anat{num_rows,4}),'unknown'); % cell array containing row index of strings with "Unknown" in u1
noneed(num_rows,5) = contains(lower(anat{num_rows,4}),'vent'); % cell array containing row index of strings with "Unknown" in u1             
end
noneed=any(noneed,2);

noneed = noneed | badch; % now is either uncessary (noneed) or bad channels


new_LL=LL; %do not clear LL, used to construct w8s
new_LL(noneed,:)=[];

new_anat=anat; %do not clear anat, used to construct noneed 
new_anat(noneed,:)=[];

new_em=em; %do not clear anat, used to construct noneed 
new_em(noneed,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4. Get variables from sem & LL functions

set(0,'DefaultFigureVisible','off')


mat_name = [study '_mat.csv'];   %string for study_mat.csv
mat_time_name = [study '_time_mat.csv']; %string for study_time_mat.csv

[SEMperiod,ll_w_t,ll_w_t_labels] = sem_plot_no_fig(mat_name,mat_time_name,vid_period,perdur); %insert filename of semiology matrix as first parameter       

[si,LL_s,ytl_LL,yt_LL,u2_s,u3_s]=LL_no_plot(new_anat,new_LL,ts,S,SEMperiod);

[clean_ll_w_t_l, clean_LL_diff] = sem_w8s_no_plot(1,ll_w_t,ll_w_t_labels,LL_s,ytl_LL,yt_LL,u2_s,sfx,perdur,vid_period,S,u3_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 5. Plot average changes 

% [symptomm,modee] = brain_w8s(pt,sz,clean_ll_w_t_l,clean_LL_diff,new_em(si,:)); 

% if fig 
% 
% figure('color','w','position',[440,348,836,449]); subplot(1,3,1:2)
% 
% lightsout; if isR; litebrain('r',1); else; litebrain('l',1); end
%     hold on; for i=1:size(new_em,1); plot3(new_em(i,1),new_em(i,2),new_em(i,3),'k.','markersize',15); end
% cmocean('balance')
% 
% colorbar('Ticks',caxis,'fontsize',18,'Position',[.65,.2,.012,.6])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% 
% % 6. Graph distribution
% if fig
% 
% subplot(1,4,4); hold on;
% end
% 
%                   symptomm=19; %make this automatically indexed later
%                   modee=2; 


[ua,~,~]=unique(new_anat(si,4));
pval=nan(length(ua),1);
anatstructureselec_weights={};
sx_i = strcmpi(clean_ll_w_t_l,sx);
chosen_sx = find(sx_i==1);
if sum(sx_i) > 0 

    for u=1:length(ua) %this is all using SORTED orders so be cautious!
      anatstructureindex=strcmpi(new_anat(si,4),ua(u));
      anatstructureselec_weights{u,1}=clean_LL_diff(anatstructureindex,chosen_sx,mx);
      EM{u} = new_em(anatstructureindex,:);
      [~,pval(u)]=ttest(anatstructureselec_weights{u,1});
%       if fig
%           if pval(u)<.05; mrkr='r*'; else; mrkr='ko'; end
%           plot(anatstructureselec_weights{u,1},u*ones(length(anatstructureselec_weights{u,1}),1),mrkr)
%       end
    end
%     if fig
%         xlim([max(abs(xlim))*[-1 1]]); yline(0,'k-'); for u=1:length(ua); xline(u,'G:',.25); end
%         set(gca,'ytick',1:length(ua),'YTickLabel',ua,'fontsize',18)
%     end
else
    anatstructureselec_weights = [];
    EM = [];
    ua = [];
    pval = [];
end
