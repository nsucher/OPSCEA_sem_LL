function [clean_ll_w_t_l, clean_LL_diff] = sem_w8s(cut_or_not,ll_w_t,ll_w_t_labels,LL_s,ytl_LL,yt_LL,u2_s,sfx,perdur,vid_period,S,u3_s)
format longG

%Created by NS & Dr. Jon Kleen, Updated 5/26 by NS


clean_ll_w_t_l = {};
label_count = 0;

ll_w_t(ll_w_t == 0) = NaN;

[LL_row, LL_col]= size(LL_s);

%subtract start of video from symptom time onset
ll_w_t = ll_w_t - vid_period(1);

LL_meandiff = nan(LL_row,1,3); % Initialize first column only 

%align ll_w_t with ll_w_t_labels
for n = 1:length(ll_w_t_labels) % loop for symptom
    if ischar(ll_w_t_labels{n})
        label_count = label_count + 1; %collect number of present labels
        clean_ll_w_t_l{label_count} = ll_w_t_labels{n}; % collect present symptoms for x labels
        if any(ll_w_t(:,n,1)) | any(ll_w_t(:,n,2)) | any(ll_w_t(:,n,3))
            if any(ll_w_t(:,n,1)) % electrode X semiology X AUTOMATISM
                if round([ll_w_t(1,n,1)+perdur]*sfx) <= LL_col %avoid exceeding array bounds   
                   %differences of the average!!! -- 
                   %mean of true start and end start - mean of true start and before start
                    LL_meandiff(:,n,1)=mean(LL_s(:,round([ll_w_t(1,n,1)       ]*sfx):round([ll_w_t(1,n,1)+perdur]*sfx)),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,1)-perdur]*sfx):round([ll_w_t(1,n,1)       ]*sfx)),2);
                elseif round([ll_w_t(1,n,1)+perdur]*sfx) > LL_col
                    LL_meandiff(:,n,1)=mean(LL_s(:,round([ll_w_t(1,n,1)       ]*sfx):LL_col                           ),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,1)-perdur]*sfx):round([ll_w_t(1,n,1)       ]*sfx)),2);
                end
            end
            if any(ll_w_t(:,n,2)) % electrode X semiology X TONIC
                  if round([ll_w_t(1,n,2)+perdur]*sfx) <= LL_col %avoid exceeding array bounds                       
                    LL_meandiff(:,n,2)=mean(LL_s(:,round([ll_w_t(1,n,2)       ]*sfx):round([ll_w_t(1,n,2)+perdur]*sfx)),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,2)-perdur]*sfx):round([ll_w_t(1,n,2)       ]*sfx)),2);
                  elseif round([ll_w_t(1,n,1)+perdur]*sfx) > LL_col
                    LL_meandiff(:,n,3)=mean(LL_s(:,round([ll_w_t(1,n,2)       ]*sfx):LL_col                           ),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,2)-perdur]*sfx):round([ll_w_t(1,n,2)       ]*sfx)),2);
                  end
            end              
            if any(ll_w_t(:,n,3)) % electrode X semiology X CLONIC
                if round([ll_w_t(1,n,3)+perdur]*sfx) <= LL_col %avoid exceeding array bounds               
                    LL_meandiff(:,n,3)=mean(LL_s(:,round([ll_w_t(1,n,3)       ]*sfx):round([ll_w_t(1,n,3)+perdur]*sfx)),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,3)-perdur]*sfx):round([ll_w_t(1,n,3)       ]*sfx)),2);
                elseif round([ll_w_t(1,n,3)+perdur]*sfx) > LL_col
                    LL_meandiff(:,n,3)=mean(LL_s(:,round([ll_w_t(1,n,3)       ]*sfx):LL_col                           ),2) ...
                                     - mean(LL_s(:,round([ll_w_t(1,n,3)-perdur]*sfx):round([ll_w_t(1,n,3)       ]*sfx)),2);
                end
            end
        else
              LL_meandiff(:,n,:) = nan;
        end
    else
           LL_meandiff(:,n,:) = nan;
    end
end


[~,col,page] = size(LL_meandiff);

clean_LL_diff = nan(LL_row,label_count,3); %initialize 3d matrix of symptoms present

clean_col = 0;

for c = 1:col
    if any(LL_meandiff(:,c,1)) | any(LL_meandiff(:,c,2)) | any(LL_meandiff(:,c,3))
        clean_col = clean_col + 1;   
        if any(LL_meandiff(:,c,1))
            clean_LL_diff(:,clean_col,1) = LL_meandiff(:,c,1);            
        end
        if any(LL_meandiff(:,c,2))
            clean_LL_diff(:,clean_col,2) = LL_meandiff(:,c,2);            
        end
        if any(LL_meandiff(:,c,3))
            clean_LL_diff(:,clean_col,3) = LL_meandiff(:,c,3);          
        end
    end
end

%EXPORT CLEAN_LL_DIFF TO BRAIN W8S TO PLOT W8 CHANGE FOR SPECIFIC SYMPTOM
%writematrix(clean_LL_diff,'q8_clean_LL_diff.csv') %need to make specific for file

 
%GROUP ELECTRODES BY YTL_LL NEUROANAT LABELS
% NEURO ANAT - overall decrease or increase 

[c_row, c_col,c_page] = size(clean_LL_diff);


% rows = # of neuroanat labels
% cols = # of symptoms
% pages = auto, tonic, clonic
% cell = weight change before and after first symptom occurance
anat_w8s = zeros(length(ytl_LL),c_col,c_page);
elec_w8s = zeros(c_row,c_col,3);

y_count = 0;

%anat_elec: 
%4 dimensional matrix to collect electrode weight changes of all symptoms in auto, tonic,& clonic

%col 1              col 2               col  3
% label number     anatomy string       electrode weight change

% 3D page 1 - # of symptoms
% each symptom

% 4D page 1 - 3 of auto, tonic, clonic

anat_elec = cell(c_row,3,c_col,c_page); %112x2 anat label number, anat label text for each electrode in row

for pages = 1:c_page
    for columns = 1:c_col
        for y_label = 1:length(ytl_LL)
            if y_label < length(ytl_LL)
               for elec = u2_s(y_label):u2_s(y_label+1)
                   y_count = y_count + 1;                  
                   anat_elec{elec,1,columns,pages} = u3_s(elec);
                   if u3_s(elec) == y_label
                       anat_elec{elec,2,columns,pages} = ytl_LL{y_label};
                       anat_elec{elec,3,columns,pages} = clean_LL_diff(elec,columns,pages);
                       anat_w8s(y_label,columns,pages) = anat_w8s(y_label,columns,pages) + anat_elec{elec,3,columns,pages};
%                        if anat_elec{elec,3,columns,pages} > 0
%                            %group into vector? 
%                        elseif anat_elec{elec,3,columns,pages} < 0
%                        else
%                        end
                   end
               end
           else
               for elec = u2_s(y_label):c_row
                   y_count = y_count + 1;                  
                   anat_elec{elec,1,columns,pages} = u3_s(elec);
                   anat_elec{elec,2,columns,pages} = ytl_LL{y_label};
                   anat_elec{elec,3,columns,pages} = clean_LL_diff(elec,columns,pages);
                   anat_w8s(y_label,columns,pages) = (anat_w8s(y_label,columns,pages) + anat_elec{elec,3,columns,pages})/y_count;
               end
            end
            y_count = 0;
        end
    end
end


        
% Plot

% figure; 
% 
% for i=1:c_page
%     subplot(c_page,1,i); 
%     pcolorjk(squeeze(clean_LL_diff(:,:,i))); 
%     colorbar; 
%     caxis(S.cax); 
%     shading flat; 
%     cmocean('balance');
%     sgtitle('Electrode Weights')
%     sgt.FontSize = 20;
%     %y
%     set(gca,'ytick',yt_LL,'yticklabel',ytl_LL,'ydir','reverse','yaxislocation','right','fontsize',8)
% 
%     %x
%     set(gca,'xtick',[1.5:1:length(clean_ll_w_t_l)+.5],'xticklabel',clean_ll_w_t_l)
%     hold on;
%     plot(xlim,[u2_s u2_s],'k-')
%     yline([min(xlim):1:max(xlim)],'k-',.01)
%     if i == 1
%         title('Automatism')
%     elseif i == 2
%         title('Tonic')
%     elseif i == 3
%         title('Clonic')
%     end
% end 
