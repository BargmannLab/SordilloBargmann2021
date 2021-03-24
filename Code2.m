% SordilloCode2 -- Sordillo A and Bargmann CI

function[RevLengths,RevSpeeds,RevDurations,...
    FwdDurationsFiltered, FwdSpeedsFiltered,...
    ROLengths,ROSpeeds,RODurations,...
    PureRevLengths,PureRevSpeeds, PureRevDurations]=get_parameters(Tracks,t1,t2,FrameRate)

% Inputs:
% Tracks - an array that contains the output of .linkedTracks, 
% generated using join_multiple_sansfood_track_arrays on data structures
% generated by BargmannWormTracker
% see Lopez-Cruz et al., 2019
%t1,t2- beggining and end of time window of interest (in minutes)
% FrameRate - framerate of acquisition

%Outputs:
%Vectors containing Reveral and Forward Run Parameters (length, speed, duration)
% in a given time window

% 8100 - total number of frames recorded
RevLengths=NaN(size(Tracks,2),8100); %create empty array for all reversal lengths
RevSpeeds=NaN(size(Tracks,2),8100); %create empty array for all reversal speeds
RevDurations=NaN(size(Tracks,2),8100); %create empty array for all reversal durations

FwdSpeeds=NaN(size(Tracks,2),8100); %create empty array for all forward run speeds
FwdDurations=NaN(size(Tracks,2),8100); %create empty array for all forward run durations

ROLengths=NaN(size(Tracks,2),8100); %create empty array for reversal omega lengths
ROSpeeds=NaN(size(Tracks,2),8100); %create empty array for reversal omega speeds
RODurations=NaN(size(Tracks,2),8100); %create empty array for all reversal omega durations

PureRevLengths=NaN(size(Tracks,2),8100); %create empty array for pure reversal lengths
PureRevSpeeds=NaN(size(Tracks,2),8100); %create empty array for pure reversal speeds
PureRevDurations=NaN(size(Tracks,2),8100); %create empty array for all pure reversal durations


% First go through Tracks and find reversals and forward runs for the time
% window.  Calculate parameters.
for h = 1:size(Tracks,2)
    c=1;
    d=1;
    if length(Tracks(h).Frames) > (5*60*FrameRate) % only consider tracks that are at least 5 minutes long
        for i=2:size(Tracks(h).Frames,2)
            if Tracks(h).Frames(i) >= (FrameRate*60*t1) %track must being at or after t1
                % considers all reversals based on state numbers defined by the
                % tracker. finds the beginning of the reversal
            if  (abs(Tracks(h).State(i-1)-4.7)>=0.1  && abs(Tracks(h).State(i)-4.7)<0.1 )|| (abs(Tracks(h).State(i-1)-5.7)>=0.1  && abs(Tracks(h).State(i)-5.7)<0.1 )...
                ||( abs(Tracks(h).State(i-1)-4.3)>=0.1  && abs(Tracks(h).State(i)-4.3)<0.1 )|| (abs(Tracks(h).State(i-1)-5.3)>=0.1  && abs(Tracks(h).State(i)-5.3)<0.1) ...
                ||(abs(Tracks(h).State(i-1)-4.1)>=0.1  && abs(Tracks(h).State(i)-4.1)<0.1) ||( abs(Tracks(h).State(i-1)-5.1)>=0.1  && abs(Tracks(h).State(i)-5.1)<0.1)...
                && i+1 < size(Tracks(h).State,2)
             j=i;
            while abs(Tracks(h).State(j)-99)>=0.1 && abs(Tracks(h).State(j)-100)>=0.1 && (j+1)<=size(Tracks(h).State,2)
                j=j+1; % these state definitions refer to gaps in the tracks or tracks approaching the copper ring
           end
           if  Tracks(h).Frames(j)>(FrameRate*60*t2) %only considers clean tracks for entire time frame considered
                a=i;
                b=i+1;
                while abs(Tracks(h).State(b))==abs(Tracks(h).State(b-1)) && b+1 < size(Tracks(h).State,2)
                    b=b+1; % get frames of the reversal
                end
                % verify reversal is complete and contained within the
                % timepoints specified. 
                % calculate length,speed, duration of reversal
                if b ~= length(Tracks(h).State) &&  Tracks(h).Frames(a)>=(FrameRate*60*t1) && Tracks(h).Frames(a)<=(FrameRate*60*t2) && Tracks(h).Frames(b)>=(FrameRate*60*t1) && Tracks(h).Frames(b)<=(FrameRate*60*t2)
                    RevLengths(h,c)=((sqrt(((Tracks(h).SmoothX(b)-Tracks(h).SmoothX(a))^2)+((Tracks(h).SmoothY(b)-Tracks(h).SmoothY(a))^2)))*Tracks(h).PixelSize)/Tracks(h).Wormlength;
                    RevSpeeds(h,c)=((nanmean(Tracks(h).Speed(a:b)))+(nanmedian(Tracks(h).Speed(a:b))))/2;   
                    RevDurations(h,c)=(Tracks(h).Frames(b)/FrameRate)-(Tracks(h).Frames(a)/FrameRate);
                    c=c+1;                   
                else 
                end
           else
           end
           elseif abs(Tracks(h).State(i-1)-1)>=0.1  && abs(Tracks(h).State(i)-1)<0.1 && i+1 < size(Tracks(h).State,2)
                 j=i;
            while abs(Tracks(h).State(j)-99)>=0.1 && abs(Tracks(h).State(j)-100)>=0.1 && (j+1)<=size(Tracks(h).State,2)
                j=j+1; % these state definitions refer to gaps in the tracks or tracks approaching the copper ring
           end
           if Tracks(h).Frames(j)>(FrameRate*60*t2) %only considers clean tracks for entire time frame considered
                a=i;
                b=i+1;
                while abs(Tracks(h).State(b))==abs(Tracks(h).State(b-1)) && b+1 < size(Tracks(h).State,2)
                    b=b+1; % get frames of forward run
                end
                  if b ~= length(Tracks(h).State) && Tracks(h).Frames(a)>=(FrameRate*60*t1) && Tracks(h).Frames(a)<=(FrameRate*60*t2) && Tracks(h).Frames(b)>=(FrameRate*60*t1) && Tracks(h).Frames(b)<=(FrameRate*60*t2)
                  % verify reversal is complete and contained within the
                  % timepoints specified. 
                  % calculate duration and speed of reversal forward run
                    FwdDurations(h,d)=(Tracks(h).Frames(b)/FrameRate)-(Tracks(h).Frames(a)/FrameRate);
                    FwdSpeeds(h,d)=((nanmean(Tracks(h).Speed(a:b)))+(nanmedian(Tracks(h).Speed(a:b))))/2;  
                    d=d+1;
                  else
                  end
           else
            end
            else
            end
            else
            end
        end
    else
    end
end

% Second loop is to parse the different reversal types for the same time
% frame.  Calculate parameters.
for h = 1:size(Tracks,2)
    c=1;

    if length(Tracks(h).Frames) > (5*60*FrameRate) % only consider tracks that are at least 5 minutes long
        for i=2:size(Tracks(h).Frames,2)
            if Tracks(h).Frames(i) >= (FrameRate*60*t1) %track must being at or after t1
                % parses reversal omegas based on state numbers defined by the
                % tracker. finds the beginning of the reversal
            if  (abs(Tracks(h).State(i-1)-4.7)>=0.1  && abs(Tracks(h).State(i)-4.7)<0.1 )|| (abs(Tracks(h).State(i-1)-5.7)>=0.1  && abs(Tracks(h).State(i)-5.7)<0.1 )...
                && i+1 < size(Tracks(h).State,2)
             j=i;
            while abs(Tracks(h).State(j)-99)>=0.1 && abs(Tracks(h).State(j)-100)>=0.1 && (j+1)<=size(Tracks(h).State,2)
                j=j+1;  % these state definitions refer to gaps in the tracks or tracks approaching the copper ring
           end
           if  Tracks(h).Frames(j)>(FrameRate*60*t2) %only consider clean tracks
                a=i;
                b=i+1;
                while abs(Tracks(h).State(b))==abs(Tracks(h).State(b-1)) && b+1 < size(Tracks(h).State,2)
                    b=b+1; % get frames of the reversal
                end
                % verify reversal is complete and contained within the
                % timepoints specified. 
                %calculate length and speed of reversal
                if b ~= length(Tracks(h).State) &&  Tracks(h).Frames(a)>=(FrameRate*60*t1) && Tracks(h).Frames(a)<=(FrameRate*60*t2) && Tracks(h).Frames(b)>=(FrameRate*60*t1) && Tracks(h).Frames(b)<=(FrameRate*60*t2)
                    ROLengths(h,c)=((sqrt(((Tracks(h).SmoothX(b)-Tracks(h).SmoothX(a))^2)+((Tracks(h).SmoothY(b)-Tracks(h).SmoothY(a))^2)))*Tracks(h).PixelSize)/Tracks(h).Wormlength;
                    ROSpeeds(h,c)=((nanmean(Tracks(h).Speed(a:b)))+(nanmedian(Tracks(h).Speed(a:b))))/2; 
                    RODurations(h,c)=(Tracks(h).Frames(b)/FrameRate)-(Tracks(h).Frames(a)/FrameRate);
                    c=c+1;
                else 
                end
           else
           end 
            % parses pure reversals based on state numbers defined by the
                % tracker. finds the beginning of the reversal
           elseif  (abs(Tracks(h).State(i-1)-4.1)>=0.1  && abs(Tracks(h).State(i)-4.1)<0.1 )|| (abs(Tracks(h).State(i-1)-5.1)>=0.1  && abs(Tracks(h).State(i)-5.1)<0.1 )...
                && i+1 < size(Tracks(h).State,2)
             j=i;
            while abs(Tracks(h).State(j)-99)>=0.1 && abs(Tracks(h).State(j)-100)>=0.1 && (j+1)<=size(Tracks(h).State,2)
                j=j+1; % these state definitions refer to gaps in the tracks or tracks approaching the copper ring
           end
           if  Tracks(h).Frames(j)>(FrameRate*60*t2) %only consider clean tracks
                a=i;
                b=i+1;
                while abs(Tracks(h).State(b))==abs(Tracks(h).State(b-1)) && b+1 < size(Tracks(h).State,2)
                    b=b+1;
                end
                % verify reversal is complete and contained within the
                % timepoints specified. 
                %calculate length and speed of reversal
                if b ~= length(Tracks(h).State) &&  Tracks(h).Frames(a)>=(FrameRate*60*t1) && Tracks(h).Frames(a)<=(FrameRate*60*t2) && Tracks(h).Frames(b)>=(FrameRate*60*t1) && Tracks(h).Frames(b)<=(FrameRate*60*t2)
                    PureRevLengths(h,c)=((sqrt(((Tracks(h).SmoothX(b)-Tracks(h).SmoothX(a))^2)+((Tracks(h).SmoothY(b)-Tracks(h).SmoothY(a))^2)))*Tracks(h).PixelSize)/Tracks(h).Wormlength;
                    PureRevSpeeds(h,c)=((nanmean(Tracks(h).Speed(a:b)))+(nanmedian(Tracks(h).Speed(a:b))))/2; 
                    PureRevDurations(h,c)=(Tracks(h).Frames(b)/FrameRate)-(Tracks(h).Frames(a)/FrameRate);
                    c=c+1;
                    
                else 
                end
           else
           end
           else
           end
           else
            end
        end
    else
    end
end


% put all data into one dimensional arrays and remove nans
RevSpeeds=RevSpeeds(:);
RevLengths=RevLengths(:);
RevDurations=RevDurations(:);
RevSpeeds=RevSpeeds(~isnan(RevSpeeds));
RevLengths=RevLengths(~isnan(RevLengths));
RevDurations=RevDurations(~isnan(RevDurations));

ROSpeeds=ROSpeeds(:);
ROLengths=ROLengths(:);
RODurations=RODurations(:);
ROSpeeds=ROSpeeds(~isnan(ROSpeeds));
ROLengths=ROLengths(~isnan(ROLengths));
RODurations=RODurations(~isnan(RODurations));

PureRevSpeeds=PureRevSpeeds(:);
PureRevLengths=PureRevLengths(:);
PureRevSpeeds =PureRevSpeeds(~isnan(PureRevSpeeds));
PureRevLengths =PureRevLengths(~isnan(PureRevLengths));
PureRevDurations =PureRevDurations(~isnan(PureRevDurations));

FwdSpeeds=FwdSpeeds(:);
FwdSpeeds =FwdSpeeds(~isnan(FwdSpeeds));
FwdDurations=FwdDurations(:);
FwdDurations =FwdDurations(~isnan(FwdDurations));

% copy forward run data to filter runs from
FwdDurationsFiltered=FwdDurations;
FwdSpeedsFiltered=FwdSpeeds;
% remove forward runs less than 2 s in duration from data sets
FwdDurationsFiltered(FwdDurationsFiltered<2) = NaN;
FwdDurationsFiltered=FwdDurationsFiltered(~isnan(FwdDurationsFiltered));
FwdSpeedsFiltered=FwdSpeedsFiltered(~isnan(FwdDurationsFiltered));

end