function sem_w8s(cut_or_not,ll_w_t,ll_w_t_labels,SEMperiod,LL_s,ytl_LL,yt_LL,u2_s,sfx)
format longG

%for entirety of sem period
LL_start_col = round(SEMperiod(1) * sfx);
LL_end_col = round(SEMperiod(2) * sfx);

%for only duration of each symptom, each element = 1 time point 
ll_true = nan(1,length(ll_w_t_labels));
%ll_true = []; %actual starting point for symptom
ll_start = nan(1,length(ll_w_t_labels));
%ll_start = []; %2 sec before
ll_end = nan(1,length(ll_w_t_labels));
% ll_end = [];

% LL_vec = [];
LL_vec=nan(size(LL_s,1),length(ll_w_t_labels),9); % Initialize 3D Matrix

% LL_diff = [];
LL_diff=nan(size(LL_s,1),length(ll_w_t_labels),3); % Initialize 3D Matrix 

LL_mean = nan(size(LL_vec,1),length(ll_w_t_labels),6);

clean_ll_w_t_l = {};
%   Page 1: auto diff
%   Page 2: tonic diff
%   Page 3: clonic diff
label_count = 0;
%align ll_w_t with ll_w_t_labels
for n = 1:length(ll_w_t_labels) % loop for symptom
    if isstring(ll_w_t_labels{n})
        label_count = label_count + 1; %collect number of present labels
        clean_ll_w_t_l{label_count} = ll_w_t_labels{n}; % collect present symptoms for x labels
        %disp(' ')
        %disp(ll_w_t_labels{n})
        %disp(' ')
        if any(ll_w_t(:,n,1)) % electrode X semiology X AUTOMATISM
            ll_true(n) = round(ll_w_t(1,n,1) * sfx);
            ll_start(n) = round(ll_w_t(2,n,1) * sfx);
            ll_end(n) = round(ll_w_t(3,n,1) * sfx);
% 
% 
%             ll_true(n) = round(ll_w_t(1,n,1));% * sfx;
%             ll_start(n) = round(ll_w_t(2,n,1));% * sfx;
%             ll_end(n) = round(ll_w_t(3,n,1));% * sfx;

%             disp('Automatism:')
%             disp(['start time: ', num2str(ll_w_t(1,n,1))])
%             disp(['2 seconds before: ', num2str(ll_w_t(2,n,1)),', ', num2str(ll_start(n))])
%             disp(['2 seconds after: ', num2str(ll_w_t(3,n,1)),', ', num2str(ll_end(n))])
            
            % rows = electrode number 
            % columns = n (symptom)
            % pages 
            % 1 = elec weight @ auto true start
            % 2 = elec weight @ auto start (x seconds before true start)
            % 3 = elec weight @ auto end (x seconds after true start)
            LL_vec(:,n,1)= LL_s(:,ll_true(n));
            LL_vec(:,n,2) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
            LL_vec(:,n,3) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix

            
%             LL_vec(:,n,1) = LL_s(:,ll_true(n));
%             LL_vec(:,n,2) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
%             LL_vec(:,n,3) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix

            for m = 1:size(LL_vec,1)
                % 3D matrix of electrode weights
                % page 1 = auto (avg true & BEFORE x sec)
                % page 2 = auto (avg true & AFTER x sec)
                LL_mean(m,n,1) = mean([LL_vec(m,n,1),LL_vec(m,n,2)]); 
                LL_mean(m,n,2) = mean([LL_vec(m,n,1),LL_vec(m,n,3)]);  
            end

            LL_diff(:,n,1) = LL_mean(:,n,2)-LL_mean(:,n,1); %(avg true & before) - (avg true & after)

        elseif any(ll_w_t(:,n,2)) % electrode X semiology X TONIC
            disp(num2str(n))
            ll_true(n) = round(ll_w_t(1,n,2) * sfx);
            ll_start(n) = round(ll_w_t(2,n,2) * sfx);
            ll_end(n) = round(ll_w_t(3,n,2) * sfx);    



            %ll_true(n)
            %ll_true = [ll_true, round(ll_w_t(1,n,2))]
%             ll_true(n) = round(ll_w_t(1,n,2));% * sfx
%             ll_start(n) = round(ll_w_t(2,n,2));% * sfx;
%             ll_end(n) = round(ll_w_t(3,n,2));% * sfx;    
%             disp(n)
            %disp(ll_start)
            %disp(ll_start(n))
%             disp('Tonic:')
%             disp(['start time: ', num2str(ll_w_t(1,n,2))])
%             disp(['2 seconds before: ', num2str(ll_w_t(2,n,2)), ', ' num2str(ll_start(n))])
%             disp(['2 seconds after: ', num2str(ll_w_t(3,n,2)), ', ' num2str(ll_end(n))])


            % rows = electrode number 
            % columns = n
            % pages 
            % 4 = elec weight @ tonic true start
            % 5 = elec weight @ tonic start
            % 6 = elec weight @ tonic end
%             disp(num2str(n))
%             disp(num2str(ll_true(n)))
            %disp(num2str(ll_true(n)*sfx))
            %disp(num2str(LL_s(:,ll_true(n))))
            LL_vec(:,n,4)= LL_s(:,ll_true(n));
            LL_vec(:,n,5) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
            LL_vec(:,n,6) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix

%             LL_vec(:,n,4) = LL_s(:,ll_true(n));
%             LL_vec(:,n,5) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
%             LL_vec(:,n,6) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix


            %LL_diff
            %LL_diff(:,n,2) = (LL_vec(:,n,4) + LL_vec(:,n,6))/2;

            %differences of the average!!! -- 
            %mean of true start and end start - mean of true start and before start
           
            for m = 1:size(LL_vec,1)
                % 3D matrix of electrode weights
                % page 3 = tonic (avg true & BEFORE x sec)
                % page 4 = tonic (avg true & AFTER x sec)

                LL_mean(m,n,3) = mean([LL_vec(m,n,4),LL_vec(m,n,5)]); 
                LL_mean(m,n,4) = mean([LL_vec(m,n,4),LL_vec(m,n,6)]);  
            end

            LL_diff(:,n,2) = LL_mean(:,n,4)-LL_mean(:,n,3); %(avg true & before) - (avg true & after)

        elseif any(ll_w_t(:,n,3)) % electrode X semiology X CLONIC
            ll_true(n) = round(ll_w_t(1,n,3) * sfx);
            ll_start(n) = round(ll_w_t(2,n,3) * sfx);
            ll_end(n) = round(ll_w_t(3,n,3) * sfx);


%             ll_true(n) = round(ll_w_t(1,n,3));% * sfx;
%             ll_start(n) = round(ll_w_t(2,n,3));% * sfx;
%             ll_end(n) = round(ll_w_t(3,n,3));% * sfx;
% 
%             disp('Clonic:')
%             disp(['start time: ', num2str(ll_w_t(1,n,3))])
%             disp(['2 seconds before: ', num2str(ll_w_t(2,n,3)),', ', num2str(ll_start(n))])
%             disp(['2 seconds after: ', num2str(ll_w_t(3,n,3)),', ', num2str(ll_end(n))])  

            % rows = electrode number 
            % columns = symptom (n)
            % pages 
            % 7 = elec weight @ clonic true start
            % 8 = elec weight @ clonic start
            % 9 = elec weight @ clonic end

            LL_vec(:,n,7)= LL_s(:,ll_true(n));
            LL_vec(:,n,8) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
            LL_vec(:,n,9) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix


%             LL_vec(:,n,7) = LL_s(:,ll_true(n));
%             LL_vec(:,n,8) = LL_s(:,ll_start(n)); %x second before values on page 1 of 3d matrix
%             LL_vec(:,n,9) = LL_s(:,ll_end(n)); %x second after values on page 1 of 3d matrix

                % page 5 = clonic (avg true & BEFORE x sec)
                % page 6 = clonic (avg true & AFTER x sec)

            for m = 1:size(LL_vec,1)                
                LL_mean(m,n,5) = mean([LL_vec(m,n,7),LL_vec(m,n,8)]); 
                LL_mean(m,n,6) = mean([LL_vec(m,n,7),LL_vec(m,n,9)]);  
            end

            LL_diff(:,n,3) = LL_mean(:,n,5)-LL_mean(:,n,6); %(avg true & before) - (avg true & after)

        else
            ll_true(n) = nan;
            ll_start(n) = nan;
            ll_end(n) = nan;
        end
    else
        ll_true(n) = nan;% * sfx;      
        ll_start(n) = nan;
        ll_end(n) = nan;
    end
end


[row,col,page] = size(LL_diff);
clean_LL_diff = []; %initialize 3d matrix of symptoms present


% loop to delete unused symptoms in all pages 
noNan = 0;
pageNum = 0;
for c = 1:col
    for p = 1:page    
        if any(LL_diff(:,c,p))
            noNan = noNan + 1;
            clean_LL_diff(:,noNan,p) = LL_diff(:,c,p);
        end        
    end
end


            
% Intialize variables for plotting

[~,c_col,c_page] = size(clean_LL_diff);

figure; 

for i=1:c_page
    subplot(3,1,i); 
%      for c = 1:c_col
%           if isnan(clean_LL_diff(:,c,i))
%               disp(['col: ',num2str(c)])
%               disp(['page: ',num2str(i)])
%               disp(' ')
%               clean_LL_diff(:,c,i) = 0;
%           end
%      end
    pcolorjk(squeeze(clean_LL_diff(:,:,i))); 
    pause
    colorbar; 
    pause
    caxis_vec = [-1 1]*120;
    caxis(caxis_vec); 
    pause
    shading flat; 
    cmocean('balance');
    pause
    %colormap(bone)
    %xlabel('Symptoms')
    sgtitle('Electrode Weights')
    sgt.FontSize = 20;
    %y
    set(gca,'ytick',yt_LL,'yticklabel',ytl_LL,'ydir','reverse','yaxislocation','right','fontsize',8)
    %x
    set(gca,'xtick',[1.5:1:length(clean_ll_w_t_l)+.5],'xticklabel',clean_ll_w_t_l)
    hold on;
    plot(xlim,[u2_s u2_s],'k-')
    yline([min(xlim):1:max(xlim)],'k-',.01)
    pause
    if i == 1
        title('Automatism')
    elseif i == 2
        title('Tonic')
    elseif i == 3
        title('Clonic')
    end
end 

disp(noNan)



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