%% Return Activity Locations for labeled data (in-lab)
% Inputs: 
%   dataTS - nx1 array of data timestamps (only for activities to label)
%   locations - nx1 cell array of activity locations (from labels)
%   activityStartEnd - nx2 array of start and end times of activities
% Output:
%   dataLocations - nx1 cell array of locations assigned to activities


function dataLocations=getActivityLocations(dataTS, locations, activityStartEnd)

    dataLocations = cell(length(dataTS),1);

    starts=activityStartEnd(:,1);
    ends=activityStartEnd(:,2);
    
    pouchinds=false(length(dataTS),1);
    pocketinds=false(length(dataTS),1);
    handinds=false(length(dataTS),1);
    
    pclips=strcmp('pouch',locations);
    pstart=find(diff([0; pclips])==1);
    pend=find(diff([0; pclips])==-1)-1;
    for i=1:length(pstart)
        if i>length(pend)
            startt=starts(pstart(i));
            endt=ends(pstart(i));
        else
            startt=starts(pstart(i));
            endt=ends(pend(i));
        end
        
        pouchinds=pouchinds | all(dataTS>startt & dataTS<endt,2);
    end
        
    pclips=strcmp('pocket',locations);
    pstart=find(diff([0; pclips])==1);
    pend=find(diff([0; pclips])==-1)-1;
    for i=1:length(pstart)
        if i>length(pend)
            startt=starts(pstart(i));
            endt=ends(pstart(i));
        else
            startt=starts(pstart(i));
            endt=ends(pend(i));
        end
        
        pocketinds=pocketinds | all(dataTS>startt & dataTS<endt,2);
    end
    
    handinds(~(pocketinds | pouchinds))=true;
    
    dataLocations(pouchinds)={'pouch'};
    dataLocations(pocketinds)={'pocket'};
    dataLocations(handinds)={'hand'};
    
end