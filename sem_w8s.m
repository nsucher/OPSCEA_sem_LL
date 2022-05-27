function sem_w8s(cut_or_not,ll_w_t,ll_w_t_labels,SEMperiod,LL_s,ytl_LL,yt_LL,u2_s,sfx,perdur,vid_period)
format longG

%Created by NS & Dr. Jon Kleen, Updated 5/26 by NS


%for entirety of sem period
LL_start_col = round(SEMperiod(1) * sfx);
LL_end_col = round(SEMperiod(2) * sfx);

%for only duration of each symptom, each element = 1 time point 
ll_true = nan(1,length(ll_w_t_labels));
ll_start = nan(1,length(ll_w_t_labels));
%ll_start = []; %2 sec before
ll_end = nan(1,length(ll_w_t_labels));
% ll_end = [];

% LL_vec = [];
LL_vec=nan(size(LL_s,1),length(ll_w_t_labels),9); % Initialize 3D Matrix

% LL_diff = [];
% LL_diff=nan(size(LL_s,1),length(ll_w_t_labels),3); % Initialize 3D Matrix 
% 
% LL_mean = nan(size(LL_vec,1),length(ll_w_t_labels),6);


clean_ll_w_t_l = {};
%   Page 1: auto diff
%   Page 2: tonic diff
%   Page 3: clonic diff
label_count = 0;

ll_w_t(ll_w_t == 0) = NaN;

[LL_row, LL_col]= size(LL_s);

%subtract start of video from symptom time onset
ll_w_t = ll_w_t - vid_period(1);

%align ll_w_t with ll_w_t_labels
for n = 1:length(ll_w_t_labels) % loop for symptom
    if isstring(ll_w_t_labels{n})
        LL_meandiff(LL_row,n,3) = nan; 

        label_count = label_count + 1' %collect number of present labels
        clean_ll_w_t_l{label_count} = ll_w_t_labels{n}; % collect present symptoms for x labels
        %disp(' ')
        %disp(ll_w_t_labels{n})
        %disp(' ')
        if any(ll_w_t(:,n,1)) % electrode X semiology X AUTOMATISM
            disp(['Auto =', ll_w_t_labels{n}])
            disp(['n =', num2str(n)])
            disp(['ll_w_t(1,n,1) = ', num2str(ll_w_t(1,n,1))])
            
            LL_meandiff(:,n,1)=mean(LL_s(:,round([ll_w_t(1,n,1)       ]*sfx):round([ll_w_t(1,n,1)+perdur]*sfx)),2) ...
                             - mean(LL_s(:,round([ll_w_t(1,n,1)-perdur]*sfx):round([ll_w_t(1,n,1)       ]*sfx)),2);

        elseif any(ll_w_t(:,n,2)) % electrode X semiology X TONIC

            disp(['Tonic =', ll_w_t_labels{n}])
            disp(['n =', num2str(n)])
            disp(['ll_w_t(1,n,2) = ', num2str(ll_w_t(1,n,2))])
            %differences of the average!!! -- 
            %mean of true start and end start - mean of true start and before start
            
            LL_meandiff(:,n,2)=mean(LL_s(:,round([ll_w_t(1,n,2)       ]*sfx):round([ll_w_t(1,n,2)+perdur]*sfx)),2) ...
                             - mean(LL_s(:,round([ll_w_t(1,n,2)-perdur]*sfx):round([ll_w_t(1,n,2)       ]*sfx)),2);

           
        elseif any(ll_w_t(:,n,3)) % electrode X semiology X CLONIC
            disp(['Clonic =', ll_w_t_labels{n}])
            disp(['n =', num2str(n)])
            disp(['ll_w_t(1,n,2) = ', num2str(ll_w_t(1,n,2))])
            LL_meandiff(:,n,3)=mean(LL_s(:,round([ll_w_t(1,n,3)       ]*sfx):round([ll_w_t(1,n,3)+perdur]*sfx)),2) ...
                             - mean(LL_s(:,round([ll_w_t(1,n,3)-perdur]*sfx):round([ll_w_t(1,n,3)       ]*sfx)),2);

        else
            LL_meandiff(:,n,:) = nan;
        end
    else
           LL_meandiff(:,n,:) = nan;
    end
end


%[row,col,page] = size(LL_diff);
[row,col,page] = size(LL_meandiff);

clean_LL_diff = []; %initialize 3d matrix of symptoms present


% loop to delete unused symptoms in all pages 
noNan = 0;
for c = 1:col
    for p = 1:page    
        if any(LL_meandiff(:,c,p))
            noNan = noNan + 1;
            clean_LL_diff(:,noNan,p) = LL_meandiff(:,c,p);
        end        
    end
end

        
% Intialize variables for plotting

[~,~,c_page] = size(clean_LL_diff);

figure; 

for i=1:c_page
    subplot(c_page,1,i); 
    pcolorjk(squeeze(clean_LL_diff(:,:,i))); 
    colorbar; 
    caxis_vec = [-1 1]*120;
    caxis(caxis_vec); 
    shading flat; 
    cmocean('balance');
    sgtitle('Electrode Weights')
    sgt.FontSize = 20;
    %y
    set(gca,'ytick',yt_LL,'yticklabel',ytl_LL,'ydir','reverse','yaxislocation','right','fontsize',8)
    %x
    set(gca,'xtick',[1.5:1:length(clean_ll_w_t_l)+.5],'xticklabel',clean_ll_w_t_l)
    hold on;
    plot(xlim,[u2_s u2_s],'k-')
    yline([min(xlim):1:max(xlim)],'k-',.01)
    if i == 1
        title('Automatism')
    elseif i == 2
        title('Tonic')
    elseif i == 3
        title('Clonic')
    end
end 



%-------------------------------------
%ANALYZE LL_diff 
%MAX, MIN, AVG

%for entire semiology

% 
% high_diff = [];
% low_diff = [];
% for n = 1:length(diff)
%     if diff(n) > mean(diff)
%         %need to find ytl_LL that corresponds with electrode #
%         high_diff(n) = diff(n);
%         %disp(diff(n))
%     elseif diff(n) < mean(diff)
%        low_diff(n) = diff(n);
%     end
% end
% %LL_vec(:,n,3) = ;



%NEED TO FIND A WAY TO MERGE DIFF AND DISP OF ELECTRODE NEUROANAT

% cut_or_not or all takes 0 or 1 
% 0 = do NOT cut the LL plot and ICEEG to only symptoms
% 1 = YES cut
%     Created by Natalia Sucher May 10 2022
% 1. align LL_s to sem times?
% 2. identify start and stop points in ll_w_t times