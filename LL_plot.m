function LL_plot(anat,badch,LL,ts,i)
%     Created by Natalia Sucher and Jon Kleen May 10 2022 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Electrode Activity
    %Electrode Activity   
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
            
%             if islogical(badch) ~= true % for newer formatted badch doubles
%                 badch=logical(badch) %| all(isnan(LL),2); %update badch here, before final noneed update with badch
%             end

            noneed = noneed | badch; % now is either uncessary (noneed) or bad channels

            new_LL=LL;
            new_anat=anat;
            new_LL(noneed,:)=[];
            new_anat(noneed,:)=[];

           %Sort new anatomy
           [u1,~,u3] = unique(new_anat(:,4)); % u1 = new_anat(u2); u2 = index of new_anat; u3 = index of u1
           [~,si] = sort(u3); % sort u3 by electrode index low to high
           u3_s = u3(si); % store sorted u3
           [~,u2_s] = unique(u3_s); %sort u2 (entries of first mentions) by deleting  repetitions in u3_s  
        
           %Sort linelength data 
           LL_s = new_LL(si,:); % sort LL data in ascending electrodes to plot

           abv_u1 = {};
           %Abbreviate ytick labels
           for k = 1:length(u1)
               %contains(u1{k},'superiortemporal')
               switch u1{k} 
                   case {'parstriangularis'}; abv_u1(k) = {'pt'};
                   case {'parsopercularis'}; abv_u1(k) = {'pop'};
                   case {'parsorbitalis'}; abv_u1(k) = {'por'};
                   %frontal
                   case {'rostralmiddlefrontal'}; abv_u1(k) = {'rmf'};
                   case {'caudalmiddlefrontal'}; abv_u1(k) = {'cmf'};
                   case {'lateralorbitofrontal'}; abv_u1(k) = {'lof'};
                   %case {'superiorfrontal'}; abv_u1(k) = {'sf'};
                   case {'medialorbitofrontal'}; abv_u1(k) = {'mof'};
                   %central
                   case {'precentral'}; abv_u1(k) = {'prec'};
                   case {'postcentral'}; abv_u1(k) = {'postc'}; 
                   %temporal
                   case {'middletemporal'}; abv_u1(k) = {'mtg'};
                   case {'superiortemporal'}; abv_u1(k) = {'stg'};                   
                   case {'inferiortemporal'}; abv_u1(k) = {'itg'};
                   case {'parahippocampal'}; abv_u1(k) = {'ph'};                   
                   case {'Right-Hippocampus'}; abv_u1(k) = {'rhp'};
%                    case {'Left-Hippocampus'}; abv_u1(k) = {'lhp'};
                   case {'Right-Amygdala'}; abv_u1(k) = {'ram'};
%                    case {'Left-Amygdala'}; abv_u1(k) = {'lam'};
                   case {'entorhinal'}; abv_u1(k) = {'ent'};
                   case {'fusiform'}; abv_u1(k) = {'fus'};                       
                   %poles    
                   case {'temporalpole'}; abv_u1(k) = {'tp'};
                   %case {'frontalpole'}; abv_u1(k) = {'fp'};
                   %marginal
                   case {'supramarginal'}; abv_u1(k) = {'sm'};
                   %other                       
                   case {'lingual'}; abv_u1(k) = {'lg'};
    %                case 'Right-Inf-Lat-Vent'; 
    %                case 'Right-Cerebral-White-Matter'; 
    %                case 'ctx...'; 
               end
           end
    
           %Average space between yticks for spacious labels
%            avg_yt_LL = [];
%            for i = 1:length(u2_s)
%                if i < length(u2_s)
%                 avg_yt_LL(i) = u2_s(i) + (u2_s(i+1) - u2_s(i))/2;
%                else
%                 avg_yt_LL(i) = u2_s(i-1) + (u2_s(i) - u2_s(i-1)/2);
%                end
%            end

           avg_yt_LL = [];
           for i = 2:length(u2_s)
               avg_yt_LL(i-1) = u2_s(i-1) + (u2_s(i) - u2_s(i-1))/2;
           end

           avg_yt_LL = [];
           for i = 2:length(u2_s)
                avg_yt_LL(i-1) = u2_s(i-1) + (u2_s(i) - u2_s(i-1))/2;
           end
                               
           avg_yt_LL(length(u2_s)) = u2_s(i); % added to account for last average in for loop             

           
           % variables for y axis 
           ytl_LL = abv_u1; %ytick labels for LL_s plot                    
           yt_LL = avg_yt_LL; %yticks for LL_s plot 
           
           % set new_LL and neuroanatomy labels as global
           setGlobal_sem_w8s(new_LL,LL_s,ytl_LL)

           % display plot with pcolor               
           pcolor(ts,1:size(new_LL,1),LL_s);
           shading flat;
          
           % modify y axis
                %In case of y label overlap
                %set(gca,'yticklabel',[],'ydir','reverse','yaxislocation','right','fontsize',8)
           xlabel('Time (seconds)')           
           set(gca,'ytick',yt_LL,'yticklabel',ytl_LL,'ydir','reverse','yaxislocation','right','fontsize',8)

           title('Line Length')
           tx=diff(xlim);

           hold on; plot(xlim,[u2_s u2_s],'w-')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        %Time Stamp
    % Time Marking
            %hold on; 
            clear h; %replot red line (delete or clear?)
            h=plot([ts(i) ts(i)],ylim,'r-');
            

