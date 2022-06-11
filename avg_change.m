% AVG_CHANGE        Plot average change of electrode weights as a heatmap projected onto a 3D reconstructed
%                   brain x seconds before and after symptom onset
%                   Display brain weights (brain_w8s) given patient filename (pt) and seizure folder (sz)                   

% Natalia Sucher
% Dr. Jon Kleen
% 06/11/2021
% UCSF Neurology, Epilepsy Department
% avg_change.m

function avg_change(pt,sz)


% PSEUDOCODE
% 1. Set path
% 2. Initialize variables
% 3. Get variables from other functions
% 4. Plot average changes 


%1. SET PATH FROM PARAMS
    opsceapath=['/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/'];   %path for parameters sheet
    opsceadatapath=[opsceapath 'OPSCEADATA/'];   %path for OPSCEA ICEEG and imaging data
        if ~exist(opsceadatapath,'dir'); error('Directory for your data needs to be corrected'); end
    cd(opsceapath);
    
    ptsz=[pt '_' sz]; % prefix for filenames of specific seizure
    ptpath=[opsceadatapath pt '/']; % patient's folder
    szpath= [ptpath ptsz '/']; % specific seizure's folder
    disp(['Running ' pt ', seizure ' sz '...']);


    [~,prm_allPtSz]=xlsread([opsceapath 'OPSCEAparams'],'params'); 
% 2. Initialize variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % S

    global S; % holds general parameters



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fields_SZ
    fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prm
    prm=prm_allPtSz(strcmp(pt,prm_allPtSz(:,1))&strcmp(sz,prm_allPtSz(:,2)),:);
    if isempty(prm); error(['ATTENTION: No entry exists for ' pt ' seizure ' sz ' in the params master sheet']); end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VIDperiod
    VIDstart=prm(:,strcmpi(fields_SZ,'VIDstart')); VIDstop=prm(:,strcmpi(fields_SZ,'VIDstop')); %chunk of data (seconds into ICEEG data file) to use from the whole ICEEG data clip for the video
    S.VIDperiod=[str2double(VIDstart{1}) str2double(VIDstop{1})];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BLperiod
    BLstart=prm(:,strcmpi(fields_SZ,'BLstart')); BLstop=prm(:,strcmpi(fields_SZ,'BLstop')); %chunk of data (seconds into ICEEG data file) to use for baseline (for z-score step)
    S.BLperiod=[str2double(BLstart{1}) str2double(BLstop{1})];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% load ICEEG data, and the bad channels verified for that specific data
    load([szpath ptsz])
    load([szpath ptsz '_badch']); 
    % rename and clear old format of electrode files -NS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %d
    if exist('ppEEG','var')
        d = ppEEG; clear ppEEG; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sfx
    if exist('fs','var')
        sfx = fs; clear fs; end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %badch

    if exist('bad_chs','var')
        badch = bad_chs; clear bad_chs; end; 
    if size(d,1)>size(d,2); d=d'; end % orient to channels by samples
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %anat

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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %em

    if ~exist('eleclabels','var'); 
        eleclabels=anat(:,1); 
    end
    em=elecmatrix; 
    clear elecmatrix; 
    emnan=isnan(mean(em,2)); 
    em(emnan,:)=0; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %llw
    S.llw=str2double(prm{strcmp('llw',fields_SZ)}); %default linelength window (in seconds)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %cax
    S.cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update d

    d(size(anat,1)+1:end,:)=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update badch

    badch(size(anat,1)+1:end)=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %nch 
    [nch,~]=size(d); %%%%%%%% edited JK 5/4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LL 

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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ts
    
    % make a timtestamp vector
    ts=0:1/sfx:size(d,2)*(1/sfx)-1/sfx;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %update d, LL, and ts
        % Extract the period of data to be used for the video (remove flanking data)
    
    vidperiodidx=round(S.VIDperiod(1)*sfx+1):S.VIDperiod(2)*sfx;
    
    d=d(:,vidperiodidx); ntp=length(vidperiodidx);
    LL=LL(:,vidperiodidx); 
    ts=ts(vidperiodidx); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %perdur
    perdur=16;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
% 3. Get variables from functions


    %   SEMIOLOGY PLOT
    set(0,'DefaultFigureVisible','off')
    [SEMperiod,ll_w_t,ll_w_t_labels] = sem_plot('f7_mat.csv','f7_time_mat.csv',S.VIDperiod,perdur); %insert filename of semiology matrix as first parameter       

    
    %   LL PLOT    
    [si,LL_s,ytl_LL,yt_LL,u2_s,u3_s]=LL_plot(new_anat,new_LL,ts,S,SEMperiod);
    
    
    %   SEM W8S PLOT
    [clean_ll_w_t_l, clean_LL_diff] = sem_w8s(1,ll_w_t,ll_w_t_labels,LL_s,ytl_LL,yt_LL,u2_s,sfx,perdur,S.VIDperiod,S,u3_s);
      
% 4. Plot changes

    %BRAIN W8S
    set(0,'DefaultFigureVisible','on')
    figure

    brain_w8s(pt,sz,clean_ll_w_t_l,clean_LL_diff,new_em(si,:)); 
