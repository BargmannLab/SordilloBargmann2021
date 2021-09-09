% SordilloCode3 -- Sordillo A and Bargmann CI, eLife

%this function collects all the processed data from each trace of neural
%activity and gathers it together to be analyzed
function [names_all,initial_percent_timeon_ALL, initial_on_switches_ALL, initial_onstatemean_sec_ALL, initial_deltaFmean_ALL, initial_deltaFmax_ALL, initial_onstate_ALL_mean_sec_ALL,...
          initial_offstates_ALL,initial_full_offstates_ALL,initial_onstates_ALL,initial_full_onstates_ALL,...
          corrected_percent_timeon_ALL, corrected_on_switches_ALL, corrected_onstatemean_sec_ALL, corrected_deltaFmean_ALL, corrected_deltaFmax_ALL, corrected_onstate_ALL_mean_sec_ALL,...
          corrected_offstates_ALL,corrected_full_offstates_ALL,corrected_onstates_ALL,corrected_full_onstates_ALL...
          ]=get_all_binary_info(structure)
 
%generate arrays to be filled with data from all the traces
initial_percent_timeon_ALL = [];
initial_on_switches_ALL = [];
initial_onstatemean_sec_ALL = [];
initial_deltaFmean_ALL = [];
initial_deltaFmax_ALL = [];
initial_onstate_ALL_mean_sec_ALL = [];
initial_offstates_ALL = [];
initial_full_offstates_ALL = [];
initial_onstates_ALL = [];
initial_full_onstates_ALL = [];

corrected_percent_timeon_ALL = [];
corrected_on_switches_ALL = [];
corrected_onstatemean_sec_ALL = [];
corrected_deltaFmean_ALL = [];
corrected_deltaFmax_ALL = [];
corrected_onstate_ALL_mean_sec_ALL = [];
corrected_offstates_ALL = [];
corrected_full_offstates_ALL = [];
corrected_onstates_ALL = [];
corrected_full_onstates_ALL = [];

names_all=[];


fn=fieldnames(structure);
%loop through the structure to get the data. you must know how the
%structure is nested beforehand See Sub-Function 2 for details of data in structure.
for i=1: numel(fn)
    fn1=fieldnames(structure.(fn{i}));
        
    data=structure.(fn{i}).(fn1{3});        
    binarydata_initial=structure.(fn{i}).(fn1{6}); 
    binarydata_corrected=structure.(fn{i}).(fn1{7}); 
        
        % call function to further process data of each trace. 
        % See Sub-Function 1 below
        [initial_percent_timeon,initial_on_switches,initial_onstatemean_sec,initial_deltaFmean,initial_deltaFmax,initial_onstate_ALL_mean_sec]=get_binary_info( binarydata_initial,data);
        [corrected_percent_timeon,corrected_on_switches,corrected_onstatemean_sec,corrected_deltaFmean,corrected_deltaFmax,corrected_onstate_ALL_mean_sec]=get_binary_info( binarydata_corrected,data);

        % collect ON to OFF state transitions from each trace
        %(See Sub-Function 2)
        initial_offstates=structure.(fn{i}).(fn1{8}); 
        initial_full_offstates=structure.(fn{i}).(fn1{9}); 
        corrected_offstates=structure.(fn{i}).(fn1{10}); 
        corrected_full_offstates=structure.(fn{i}).(fn1{11}); 
        
        initial_onstates=structure.(fn{i}).(fn1{12}); 
        initial_full_onstates=structure.(fn{i}).(fn1{13}); 
        corrected_onstates=structure.(fn{i}).(fn1{14}); 
        corrected_full_onstates=structure.(fn{i}).(fn1{15}); 
        
        %gets name of each trace
        names_all=[names_all,fn{i}];
        
        % put all the data together for further analysis
        initial_percent_timeon_ALL = vertcat(initial_percent_timeon_ALL, initial_percent_timeon);
        initial_on_switches_ALL = vertcat(initial_on_switches_ALL,initial_on_switches);
        initial_onstatemean_sec_ALL = vertcat(initial_onstatemean_sec_ALL,initial_onstatemean_sec);
        initial_deltaFmean_ALL = vertcat(initial_deltaFmean_ALL,initial_deltaFmean);
        initial_deltaFmax_ALL = vertcat(initial_deltaFmax_ALL,initial_deltaFmax);
        initial_onstate_ALL_mean_sec_ALL = vertcat(initial_onstate_ALL_mean_sec_ALL,initial_onstate_ALL_mean_sec);
        
        if isempty(initial_offstates)==0
        initial_offstates_ALL = horzcat(initial_offstates_ALL,initial_offstates);
        else
        end
            
        if isempty(initial_full_offstates)==0
        initial_full_offstates_ALL = horzcat(initial_full_offstates_ALL,initial_full_offstates);
        else
        end
        
        if isempty(initial_onstates)==0
        initial_onstates_ALL = horzcat(initial_onstates_ALL,initial_onstates);
        else
        end
            
        if isempty(initial_full_onstates)==0
        initial_full_onstates_ALL = horzcat(initial_full_onstates_ALL,initial_full_onstates);
        else
        end
        
        corrected_percent_timeon_ALL = vertcat(corrected_percent_timeon_ALL, corrected_percent_timeon);
        corrected_on_switches_ALL = vertcat(corrected_on_switches_ALL,corrected_on_switches);
        corrected_onstatemean_sec_ALL = vertcat(corrected_onstatemean_sec_ALL,corrected_onstatemean_sec);
        corrected_deltaFmean_ALL = vertcat(corrected_deltaFmean_ALL,corrected_deltaFmean);
        corrected_deltaFmax_ALL = vertcat(corrected_deltaFmax_ALL,corrected_deltaFmax);
        corrected_onstate_ALL_mean_sec_ALL = vertcat(corrected_onstate_ALL_mean_sec_ALL,corrected_onstate_ALL_mean_sec);
        
        if isempty(corrected_offstates)==0
        corrected_offstates_ALL = horzcat(corrected_offstates_ALL,corrected_offstates);
        else
        end
        
        if isempty(corrected_full_offstates)==0
        corrected_full_offstates_ALL = horzcat(corrected_full_offstates_ALL,corrected_full_offstates);
        else
        end
        
        if isempty(corrected_onstates)==0
        corrected_onstates_ALL = horzcat(corrected_onstates_ALL,corrected_onstates);
        else
        end
        
        if isempty(corrected_full_onstates)==0
        corrected_full_onstates_ALL = horzcat(corrected_full_onstates_ALL,corrected_full_onstates);
        else
        end

end
end
%-----%
% END %
%-----%

%---------------%
% SUB-FUNCTION 1 %
%---------------%


function [percent_timeon,on_switches,onstatemean_sec,deltaFmean,deltaFmax,onstate_ALL_mean_sec]=get_binary_info(binarydata,data);

%generate arrays of appropriate size (our data is 2999-3000 frames)
on_switches=zeros(2999,size(binarydata,2));
off_switches=zeros(2999,size(binarydata,2));
onstate_lengths=[];
end_onstate_lengths=[];
beginning_onstate_lengths=[];

for h=1:size(binarydata,2)
for i=2:2999
    %finds transitions into the ON state 
    if binarydata((i-1),h)==0 && binarydata(i,h)==1
        on_switches(i,h)=1;
        length = 1;
        a=i;
        while a< 2999 && binarydata(a+1,h)==1
            length=length+1; %counts the amount of frames 
            a=a+1;
        end
        % special case if it is cut off at the end of our time bin
        if a==2999
            end_onstate_lengths=[end_onstate_lengths length]; 
        else
            onstate_lengths=[onstate_lengths length]; %FULL EVENTS
        end
       
       %special case for data cut off by beginning of time bin
       %find transition to on state
    elseif (i-1)==1 && binarydata((i-1),h)==1 && binarydata(i,h)==1
        length = 1;
        a=i;
        while a< 2999 && binarydata(a+1,h)==1 
            length=length+1; %counts the amount of frames 
            a=a+1;
        end
          % special case if it is cut off at the end of our time bin
        if a==2999
            end_onstate_lengths=[end_onstate_lengths length];
        else
            beginning_onstate_lengths=[beginning_onstate_lengths length];
        end
        
     %finds transitions into the OFF state 
    elseif binarydata((i-1),h)==1 && binarydata(i,h)==0
       off_switches(i,h)=1;
        
    end
end
end

time_on=sum(binarydata(1:2999,:)); %sum all ON frames
percent_timeon=time_on./size(on_switches,1); % calculate the fraction time in ON state over all frames
%switches=sum(on_switches+off_switches);
on_switches=sum(on_switches); % get all the transitions to the ON state
deltaFmean=mean(data(binarydata==1));% calculate the mean deltaF while ON

%not used in this manuscript
onstate_lengths_all=horzcat(onstate_lengths,end_onstate_lengths,beginning_onstate_lengths).';
onstatemean_sec=nanmean(onstate_lengths)/10; % calculate the mean length of FULL ON events in seconds (fps = 10)
onstate_ALL_mean_sec=nanmean(onstate_lengths_all)/10;   % calculate the mean length of ALL ON events in seconds (fps = 10)
deltaFmax=max(data(binarydata==1));% calculate the max deltaF while ON

end

%-----%
% END %
%-----%

%---------------%
% SUB-FUNCTION 2 %
%---------------%

% gets data from each trace
function Output = arenachip_gcamp_aug2021()

%load data
addpath(pwd);
PathName = uigetdir;
cd(PathName);
addpath(pwd);

fps = 10;
filenames = dir('*.csv'); %csv files generated by python script.  See SordilloCode4

% Calculate deltaF (also returns other versions of the data)
for d=1:length(filenames)
    [raw,deltaF,deltaFn,deltaF_true,time] = arenachip_getdeltas(filenames(d).name,fps,1);
    
    %parse file names to get names for plots and data structures
   
    USLocations = find(filenames(d).name == '_');

    name1 =filenames(d).name((USLocations(1)+1):(USLocations(2)-1));
    name2 =filenames(d).name((USLocations(2)+1):(USLocations(3)-1));
    

    % get time derivative
    Vtotal = bin_data_thirty(tderiv(deltaF,10));
    
    % get ON and OFF states; clean up binning to remove small deltaFs and
    % very short states called due to noise, etc.
    % see each subfunction below for explanation of parameters
    [Binary_Initial,param1,param2] = binneuron(bin_data_ten(deltaF),Vtotal);
    Binary_Corrected=filter_short_states(filter_small_changes(bin_data_ten(deltaF),Binary_Initial,0.3,5,.9),20);
    
    %initial binning for plots
    Binary_Max=filter_short_states(bin_by_max(bin_data_ten(deltaF),.5),20);
    
    %returns on to off state transitions for study of dynamics (not used in this manuscript)
    [off_states_Initial,off_states_full_Initial]=get_off_states(bin_data_ten(deltaF), Binary_Initial);
    [off_states_Corrected,off_states_full_Corrected]=get_off_states(bin_data_ten(deltaF), Binary_Corrected);
    
     %returns off to on states for study of dynamics (not used in this manuscript)
    [on_states_Initial,on_states_full_Initial]=get_on_states(bin_data_ten(deltaF), Binary_Initial);
    [on_states_Corrected,on_states_full_Corrected]=get_on_states(bin_data_ten(deltaF), Binary_Corrected);
    
   
    maxF = max(bin_data_ten(deltaF));
     %Ignore low amplitude data; filter data missing more than 1 frame
     %3000 frames (5 min) were analyzed for each time point in this figure
     if maxF > 10 && length(time)>=2999 
         
        %generate Output data structure
        Output.(strcat('LS_worm ',name1)).time = time;
        Output.(strcat('LS_worm',name1)).raw = raw;
        Output.(strcat('LS_worm',name1)).deltaF = deltaF;
        Output.(strcat('LS_worm',name1)).deltaF = deltaFn;
        Output.(strcat('LS_worm',name1)).deltaF = deltaF_true;
        Output.(strcat('LS_worm',name1)).binned = bin_data_ten(raw);
        Output.(strcat('LS_worm',name1)).Binary_Max= Binary_Max;
        Output.(strcat('LS_worm',name1)).Binary_Initial= Binary_Initial;
        Output.(strcat('LS_worm',name1)).Binary_Corrected= Binary_Corrected;
        Output.(strcat('LS_worm ',name1)).off_states_Initial=off_states_Initial; 
        Output.(strcat('LS_worm ',name1)).off_states_full_Initial=off_states_full_Initial; 
        Output.(strcat('LS_worm ',name1)).off_states_Corrected=off_states_Corrected; 
        Output.(strcat('LS_worm ',name1)).off_states_full_Corrected=off_states_full_Corrected; 
        Output.(strcat('LS_worm ',name1)).on_states_Initial=on_states_Initial; 
        Output.(strcat('LS_worm ',name1)).on_states_full_Initial=on_states_full_Initial; 
        Output.(strcat('LS_worm ',name1)).on_states_Corrected=on_states_Corrected; 
        Output.(strcat('LS_worm ',name1)).on_states_full_Corrected=on_states_full_Corrected; 
        
        % plot data for troubleshooting, etc.
        figure; 
        subplot(5,1,1)
        plot(time, bin_data_ten(raw),'LineWidth',2,'color','k');
        title(join(['raw trace: worm',string(name1),string(name2)],'-'))
        %ylim([0 40])
        set(gca,'FontSize',20)
        hold off;
        
        subplot(5,1,2)
        plot(time, bin_data_ten(deltaF),'LineWidth',2,'color','k');
        hold on;plot(time,Binary_Max*max(bin_data_ten(deltaF)),'LineWidth',2,'color','b')
        title('Threshold Binning')
        set(gca,'FontSize',20)
        
        subplot(5,1,3)
        plot(time,Vtotal,'LineWidth',2); hold on;
        yline(0);
        yline(param1);
        yline(param2);
        ylim([min(Vtotal), max(Vtotal)]);
        yticks([min(Vtotal), max(Vtotal)]);
        title('Derviative')
        %ylim([0 125])
        set(gca,'FontSize',20)
        hold off;
        
        subplot(5,1,4)
        plot(time, bin_data_ten(deltaF),'LineWidth',2,'color','k');
        hold on;plot(time,Binary_Initial*max(bin_data_ten(deltaF)),'LineWidth',2,'color','r');
        title('Derviative Binning')
        %ylim([0 125])
        set(gca,'FontSize',20)
        hold off;
        
        subplot(5,1,5)
        plot(time, bin_data_ten(deltaF),'LineWidth',2,'color','k');
        hold on;plot(time,Binary_Corrected*max(bin_data_ten(deltaF)),'LineWidth',2,'color','r');
        title('Derviative Binning:  Corrected')
        %ylim([0 125])
        set(gca,'FontSize',20)
        hold off;
        
       
%         if isempty(off_states_Corrected)==0
%         figure;
%         subplot(2,1,1)
%         plot(time, bin_data_ten(deltaF),'LineWidth',2,'color','k');
%         hold on;plot(time,Binary_Corrected*max(bin_data_ten(deltaF)),'LineWidth',2,'color','r');
%         title('Derviative Binning:  Corrected')
%         %ylim([0 125])
%         set(gca,'FontSize',20)
%         hold off;
%         
%         subplot(2,1,2)
%         plot_seconds(off_states_Corrected); [row,col]=size(off_states_Corrected);h=1:col;legend(num2str(h.'));
%         set(gca,'FontSize',20)
%         hold off;
%         else
%         end
        
        
        
        
        
        

    else
    end


   
end

%-----%
% END %
%-----%
end



%---------------%
% SUB-FUNCTIONS %
%---------------%

function [raw,deltaF,deltaFn,deltaF_true, time] = arenachip_getdeltas(name,fps,col)

data = importdata(name);

d = length(data);
time=(1:d)/(fps*60);


raw=data(:,col);
F01=min(data(:,col));
deltaF=((data(:,col)-F01)/F01)*100;
NP_sort=sort(data(:,col)); 
NP(1,col)=median(NP_sort(1:round(0.1*size(NP_sort,1)),1),1); 
%median of 10% lowest data
Initial=repmat(NP(1,col),size(data(:,col),1),1);
deltaF(:,col)=(data(:,col)-Initial)./Initial*100;

deltaF_true=deltaF ;
deltaF(deltaF<0)=0;

deltaFn = deltaF./max(deltaF);


end

function Vtotal = tderiv(deltaFtotal,frames) %Adapted from Gordus et al., 2015
                                             % original code not written by
                                             % A Sordillo
            
[d1,d2] = size(deltaFtotal);

if d2>d1
    deltaFtotal = deltaFtotal'; %transpose data if in wrong orientation
end

d = size(deltaFtotal,2);

Vtotal = NaN(length(deltaFtotal),d); %generate matrix to fill
    % calculate time derivative
for m=frames+1:length(deltaFtotal)-frames
    Vtotal(m,:) = (deltaFtotal(m+frames,:)-deltaFtotal(m-frames,:))/(2*frames);
end

end



function [BinaryTotal,param1,param2] = binneuron(data,Vtotal); %Adapted from Gordus et al., 2015
                                                               % original code not written by
                                                               % A Sordillo

thresh=.5;
BinaryTotal=filter_short_states(bin_by_max(bin_data_ten(data),thresh),20);

% param1 = positive dF/dt threshold
% param2 = negative dF/dt threshold
% ucountlimit = # of frames threshold for dF/dt to be above param1
% dcountlimit = # of frames threshold for dF/dt to be below param2

d = size(data,2);

%AVA WT PARAMS  

ucountlimit = 20; 
dcountlimit = 40; 

%adjust dF/dt that must be reached to call each state base on the 
% max(Vtotal) and min(Vtotal) that makes sense for your data
% this allows different genotypes to be parsed by the same script
if max(Vtotal) >0.35 
    param1 = max(Vtotal)*.075; %param1 - dF/dt threshold for ON states
else
    param1 = max(Vtotal)*.3; 
end
if min(Vtotal) <= -0.1907 && min(Vtotal)>=-0.3007  
    param2 =min(Vtotal)*.25; %param2 - dF/dt threshold for OFF states
elseif min(Vtotal) < -0.3007
    param2 =min(Vtotal)*.25; 
    dcountlimit = 30; 
else
    param2 = min(Vtotal)*.15; 
    dcountlimit=100;
end
   

 
 for m=1:d

    dcount = 0;
    ucount = 0;
    for n=2:length(data)
        
        % determine OFF and ON states based on parameters determined above
        if BinaryTotal(n-1,m) == 0 && Vtotal(n,m) > param1
            if ucount == ucountlimit
                BinaryTotal(n-ucountlimit:n,m) = ones(ucountlimit+1,1);
                ucount = 0;
            else
                ucount = ucount+1;
                BinaryTotal(n,m) = BinaryTotal(n-1,m);
            end
        % special case made for beginning of trace, make parameters more
        % lenient 26 is the first frame that has a nonNaN Vtotal in my data
        elseif n==26 && BinaryTotal(n-1,m) == 1 &&  Vtotal(n,m) < param2*.75 
                  o=n;
                  while Vtotal(o,m) < param2*.75
                        dcount = dcount+1;
                        o=o+1;
                  end
                  if dcount>=(dcountlimit*.3) 
                     BinaryTotal(o-dcount:o,m) = zeros(dcount+1,1);
                     n=o;
                     dcount=0;
                  else
                     n=o;
                     dcount=0;
                  end
                                        
       elseif BinaryTotal(n-1,m) == 1 && Vtotal(n,m) < param2 
            if dcount == dcountlimit
                BinaryTotal(n-dcountlimit:n,m) = zeros(dcountlimit+1,1);
                dcount = 0;
            else
                dcount = dcount+1;
                BinaryTotal(n,m) = BinaryTotal(n-1,m);
                
            end      
       
      else
            BinaryTotal(n,m) = BinaryTotal(n-1,m);
            ucount = 0;
            dcount = 0;
            
        end
        
    end
 end


 end



function binned_data_ten = bin_data_ten(data)
% bin over ten frames to smooth noise in data
[a,b]=size(data);
binned_data_ten = NaN(a,b);

    

for h=1:b
    for k=1:a
        if 5<k && k<=(a-5)
            binned_data_ten(k,h) = mean(data(k-5:k+5,h));
        elseif k<=5
            binned_data_ten(k,h) = mean(data(k-(k-1):k+5,h));
        else 
            binned_data_ten(k,h) = mean(data(k-5:end,h));
        end
    end
end

end

function binned_data_twenty = bin_data_twenty(data)
% bin over twenty frames to smooth noise in data (not used in this
% manuscript)
[a,b]=size(data);
binned_data_twenty = NaN(a,b);

    

for h=1:b
    for k=1:a
        if 10<k && k<=(a-10)
            binned_data_twenty(k,h) = mean(data(k-10:k+10,h));
        elseif k<=10
            binned_data_twenty(k,h) = mean(data(k-(k-1):k+10,h));
        else 
            binned_data_twenty(k,h) = mean(data(k-10:end,h));
        end
    end
end

end

function binned_data_thirty = bin_data_thirty(data)
% bin over thirty frames to smooth noise in data 
[a,b]=size(data);
binned_data_thirty = NaN(a,b);

for h=1:b
    for k=1:a
        if 15<k && k<=(a-15)
            binned_data_thirty(k,h) = mean(data(k-15:k+15,h));
        elseif k<=15
            binned_data_thirty(k,h) = mean(data(k-(k-1):k+15,h));
        else 
            binned_data_thirty(k,h) = mean(data(k-15:end,h));
        end
    end
end

end


function hist_binned_data=bin_by_hist(data)
%can be used to call ON and OFF states using the most common bins in a
%histogram (not used in this manuscript)

hist_binned_data=zeros(size(data));
[V,E]=histcounts(data);
% Retrieve some properties from the histogram
    
peaks=maxk(V,2);
bin1=max(E(V==peaks(1)));
%bin2=min(E(V==peaks(2)));
bin2=max(E(V==peaks(2)));

bigbin=max([bin1,bin2]);
% Find the centers of the bins that maxk identified as peaks
%center = (bin1 + bin2)/2;
thresh=bigbin*.8;

%on_vals=data>center;
on_vals=data>thresh;
hist_binned_data(on_vals)=1;

end

function on_off_states=bin_by_max(data,thresh)
% bins data based on if it surpasses a threshold of the max data.  we used
% 0.5 for thresh
on_off_states=zeros(size(data));
maxdata=max(data);
thresh=maxdata*thresh;

on_vals=data>=thresh;
on_off_states(on_vals)=1;

end

function data=filter_short_states(data,thresh)
%filters out short states due to noise,etc. we used a threshold of 20
%frames

% find all points where switches between states happen
split_starts = find(diff(data)==-1);
split_ends = find(diff(data)==1);
splits=sort(vertcat(split_starts,split_ends));

% filter short states
if isempty(splits)==0
    for i=1:length(splits)-1
        if abs(splits(i)-splits(i+1))<thresh 
            data((splits(i)):(splits(i+1)))=data(splits(i)); %makes short state same as previous state
        else
        end
    end
end

end

function bin_data=filter_small_changes(data,bin_data,thresh1,thresh2,thresh3)

% find all points where switches between states happen
split_ends = find(diff(bin_data)==-1);
split_starts = find(diff(bin_data)==1);
splits=sort(vertcat(split_starts,split_ends));

%clean binary states
%thresh1 = the smallest delta F you want to consider for each event
%thresh 2= the smallest max delta F you want to consider
%thresh 2= the largest min delta F you want to consider
if length(splits)>=2
            event=data(splits(1):splits(2));
          
            change=abs(max(event)-min(event));
            if change<(max(event)*thresh1) || max(event)<thresh2 || (min(event)>thresh3*max(data) && bin_data(splits(1))==1)...
                    
                bin_data((splits(1)):(splits(2)))=bin_data(splits(1));
            else
            end
      
        for i=2:length(splits)-1
            
            prev_event=data(splits(i-1):splits(i));
            event=data(splits(i):splits(i+1));
          
            change=abs(max(event)-min(event));
            
            % attempts to filters slow decreases that occur during ON
            % states without switching...can change parameters (0.3 and
            % 100) to fit your needs
            if change<(max(event)*thresh1) || max(event)<thresh2 || (min(event)>thresh3*max(data) && bin_data(splits(i))==1)...
                    || (min(event)>0.3*max(prev_event) && bin_data(splits(i))==1 && length(event) > 100)
                bin_data((splits(i)):(splits(i+1)))=bin_data(splits(i));
            else
            end
        end
        i=length(splits);
       
        %special case for final events
            event=data(splits(i):end);
            change=abs(max(event)-min(event));
                if change<(max(event)*thresh1) || max(event)<thresh2 || (min(event)>thresh3*max(data) && bin_data(splits(i))==1)
                     
                    bin_data(splits(i):end)=bin_data(splits(i));
                else
                end
         
       

            
            
else
end
end

function [off_states,full_off_states]=get_off_states(data,bin_data)
%puts all ON to OFF transitions into a matrix for analysis (not used in
%this manuscript)
% events that had a minimum within 5% of the max were considered full
% events

split_ends = find(diff(bin_data)==-1);
split_starts = find(diff(bin_data)==1);
splits=sort(vertcat(split_starts,split_ends));
off_states=NaN(3000,size(split_ends,1));
full_off_states=NaN(3000,size(split_ends,1));


for i=1:size(split_ends,1)
    index=find(splits==split_ends(i));
    if index+1<=length(splits)
        event=data(splits(index):splits(index+1));
    else
        event=data(splits(index):end);
    end
        event_n=event./max(event);
        off_states(1:length(event_n),i)=event_n;
        if any(event_n(:) <=0.05)
            full_off_states(1:length(event_n),i)=event_n;
        else
        end
end
end



function [on_states,full_on_states]=get_on_states(data,bin_data)
%puts all OFF to ON transitions into a matrix for analysis (not used in
%this manuscript)
% events that had a minimum within 5% of the max were considered full
% events
split_ends = find(diff(bin_data)==-1);
split_starts = find(diff(bin_data)==1);
splits=sort(vertcat(split_starts,split_ends));
on_states=NaN(3000,size(split_starts,1));
full_on_states=NaN(3000,size(split_starts,1));


for i=1:size(split_starts,1)
    index=find(splits==split_starts(i));
    if index+1<=length(splits)
        event=data(splits(index):splits(index+1));
    else
        event=data(splits(index):end);
    end
        event_n=event./max(event);
        on_states(1:length(event_n),i)=event_n;
        if any(event_n(:) <=0.05)
            full_on_states(1:length(event_n),i)=event_n;
        else
        end
end
end