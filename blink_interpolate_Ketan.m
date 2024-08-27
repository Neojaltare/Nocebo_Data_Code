
function [outdata,nanids] = blink_interpolate_Ketan(data,time)

% Identify 0s and linear interpolate +/-150 ms
blinks = double(data==0);
onoff = diff(blinks); % positive values are onset and negative are offset
onsets = time(onoff>.5);
offsets = time(onoff<-.5);

% Make sure everything is identified correctly
while onsets(1)>offsets(1)
    offsets(1) = [];
end

if length(onsets) - length(offsets) == 1
    onsets = onsets(1:end-1);
end
invalids = [onsets;offsets];

% now merge any periods that are less than 250ms apart
temp = invalids(:,1);
for d = 1:size(invalids,2)-1
    if abs(invalids(2,d) - invalids(1,d+1)) < 250
        temp(2,end) = invalids(2,d+1);
    else
        temp(:,end+1) = invalids(:,d+1);
    end
end
invalids = temp; clear temp

% pad the blinks and then check for times outside the data range
invalids(1,:) = invalids(1,:) - 150;
invalids(2,:) = invalids(2,:) + 250;
outofrangeend = find(invalids(2,:)>size(data,2));
if ~isempty(outofrangeend)
    invalids(:,outofrangeend) = [];
end


% Initial rough Interpolation of blinks
for i = 1:length(invalids)
    start = invalids(1,i);
    stop = invalids(2,i);
    try
    interpolated = interp1([start stop],[data(time==start) data(time==stop)],start:stop,'linear');
    data(find(time==start):find(time==stop)) = interpolated;
    catch
        disp(['Something happened with interp on iteration ',num2str(i)])
    end
end


% Vectorise all invalid time points
invalvec = [];
for i = 1:length(invalids)
    invalvec = [invalvec invalids(1,i):invalids(2,i)];
end

outdata = data;
nanids = invalvec;
