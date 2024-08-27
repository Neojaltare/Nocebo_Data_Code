
addpath('Path to fieldtrip /fieldtrip-20220426')
ft_defaults

%% Fix this code for extracting all the relevant file names

% % For loading all the files
filepath = 'Path to folder where all the raw pupil data is stored/';
files = dir(strcat(filepath,'*.asc'));
files(1:60) = [];

for i = 1:length(files)
    disp(i)
    
    filename = files(i).name;
    filename_eye = fullfile(filepath,filename);

    cfg = [];
    cfg.dataset = filename_eye;
    data_eye = ft_preprocessing(cfg);
    
    % Open the file
    fid = fopen(filename_eye, 'rt');  % replace 'filename.asc' with your actual file name
    if fid == -1
        error('Cannot open file');
    end
    
    % Initialize an array to store the numeric data
    Events = [];
    
    % Read the file line by line
    while true
        line = fgets(fid);  % Read one line
        if line == -1
            break;  % End of file
        end
    
        % Check if the line starts with 'BUTTON'
        if strncmp(line, 'BUTTON', length('BUTTON'))
            % Use textscan to read the numeric data from the rest of the line
            % This assumes that after 'BUTTON' there are space-separated values
            lineData = textscan(line(length('BUTTON')+1:end), '%f');  % Skip 'BUTTON' and read numbers
    
            % Concatenate the lineData to the numericData array
            % lineData{1} contains the numeric data from the current line as a cell array
            Events = [Events; lineData{1}'];  % Assuming horizontal concatenation is desired
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Events now contains all the extracted numeric data after 'BUTTON' from each line
    Events = Events(Events(:,3) == 1,:);
    
    % Adjust latency for the pinprick triggers
    Events(Events(:,2)==2,1) = Events(Events(:,2)==2,1) + 250;
    
    % Put the data back in the data_eye struct
    data_eye.trial{1}(5,Events(:,1)) = Events(:,2);
    
    
    %% Pre Processing
    data = data_eye.trial{1}(4,1:length(data_eye.trial{1}(1,:)));
    time = data_eye.trial{1}(1,:);
    data(time==0) = [];
    time(time==0) = [];
    
    temptime = 1:max(time);
    tempdata = zeros(1,max(time));
    tempdata(time) = data;
    
    data = tempdata;
    time = temptime;
    datacopy = data;
    
    plotme = true;
    dat = struct;
    dat.pupil = data;
    dat.fsample = 1000;
    dat.times = time;
    
    % The following two lines call the relevant function for the processing
    % of the pupil data. The variable names may be a bit confusing. Note
    % that the functions return the processed data with the indices of
    % interpolated samples
    [proc_data_ketan,knans] = blink_interpolate_Murphy(dat, false); 
    %[proc_data_ketan,knans] = blink_interpolate_Ketan(data, time);

    
    %% low pass filter at 6Hz
    
    % Example parameters
    fs = 1000; % Sampling frequency in Hz
    fc = 6;   % Cutoff frequency in Hz
    order = 3; % Order of the filter
    
    % Normalize the frequency
    Wn = fc/(fs/2);
    
    % Create the low-pass Butterworth filter
    [b, a] = butter(order, Wn, 'low');
    
    % Apply the filter
    filtered_data = filtfilt(b, a, proc_data_ketan);
    
    %% Set interpolated segments longer than 1 second to be NANs
    
    interpidx = zeros(1,size(filtered_data,2));
    interpidx(knans) = 1;
    Connected = bwconncomp(interpidx);
    numpoints = cellfun(@numel,Connected.PixelIdxList);
    
    longs = numpoints>1000;
    selectedCells = Connected.PixelIdxList(longs);
    invalvec = vertcat(selectedCells{:});
    datcopy = filtered_data;
    filtered_data(invalvec) = nan;
    
    
    %% Extract the relevant sections of the data for RS and the event related pupil
    
    notes = {};
    dat = struct;
    excludethresh = 15;
    
    % Check number of events
    numRS = sum(Events(:,2)==1);
    numPP = sum(Events(:,2)==2);
    restingidx = Events(Events(:,2)==1,1);
    Prickidx = find(Events(:,2)==1);
    
    
    % RS1
    try
        RS1idx = restingidx(1:2);
        beg = find(time==RS1idx(1));
        fin = find(time==RS1idx(2));
        overlap = ismember(beg:fin,invalvec); % see how much of the data was interpolated
    
        if sum(overlap)<=60000*.3
            RS1 = datcopy(beg:fin);
            dat.RS1 = RS1;
        elseif sum(overlap)>60000*.3 % If it was more than 30%, make a note
            notes{end+1} = 'More than 30 percent of RS1 is nans';
        end
    catch
        notes{end+1} = 'Problem with RS1. Check the plot. ';
    end
    
    
    % RS2
    try
        RS2idx = restingidx(3:4);
        beg = find(time==RS2idx(1));
        fin = find(time==RS2idx(2));
        overlap = ismember(beg:fin,invalvec); % see how much of the data was interpolated
        if sum(overlap)<=60000*.3
            RS2 = datcopy(beg:fin);
            dat.RS2 = RS2;
        elseif sum(overlap)>60000*.3 % If it was more than 30%, make a note
            notes{end+1} = 'More than 30 percent of RS2 is nans';
        end
    
    catch
        notes{end+1} = 'Problem with RS2. Check the plot.';
    end
    
    
    try
%         T0idx = Events(Prickidx(2)+1:Prickidx(3)-1,:);
        T0idx = Events(Prickidx(2)+1:end,:);

        T0idx(:,1) = T0idx(:,1) + 242;
        T0 = zeros(length(T0idx),6000);
        index = 1;
        for i = 1:length(T0idx)
            try
                beg = find(time==(T0idx(i,1)-1000));
                fin = find(time==(T0idx(i,1)+4999));
                temp = filtered_data(beg:fin);
            catch
                warning(['False Trigger: ',num2str(i)])
            end
    
            if sum(temp) == 0
                notes{end+1} = ['Trigger ',num2str(i),' at T1 was all zeros.'];
                continue
            end
    
            if ~isempty(temp)
                overlap = ismember(beg:fin,invalvec);
                missing = isnan(temp);
                Connected = bwconncomp(missing);
                numpoints = cellfun(@numel,Connected.PixelIdxList);
                if sum(overlap)<=3000 %&&  ~any(numpoints>1000) % Check that total number of nans are less than 50 percent
                    T0(index,:) = datcopy(beg:fin);
                    index = index+1;
%                 elseif any(numpoints>1000)
%                     notes{end+1} = ['Trial ',num2str(i),' at T0 excluded. More than one second missing.'];
                end
    
                if sum(overlap)>3000
                    notes{end+1} = ['Trial ',num2str(i),' at T0 excluded. More than 50% missing.'];
    
                end
            end
    
        end
    
        T0(index:end,:) = [];
    
        dat.T0 = T0;
    
        if size(T0,1) < excludethresh
            notes{end+1} = 'More than 30 percent of trials at T0 missing.';
        end
    
    catch
        notes{end+1} = 'Markers for T0 may be missing';
    end
    
    
    
    
    try
        T1idx = Events(Prickidx(4)+1:end,:);
        T1idx(:,1) = T1idx(:,1) + 242;
        T1 = zeros(length(T1idx),6000);
        index = 1;
        for i = 1:length(T1idx)
            try
                beg = find(time==(T1idx(i,1)-1000));
                fin = find(time==(T1idx(i,1)+4999));
                temp = filtered_data(beg:fin);
            catch
                warning(['False Trigger: ',num2str(i)])
            end
    
            if sum(temp) == 0
                notes{end+1} = ['Trigger ',num2str(i),' at T2 was all zeros.'];
                continue
            end
    
    
            if ~isempty(temp)
                overlap = ismember(beg:fin,invalvec);
                missing = isnan(temp);
                Connected = bwconncomp(missing);
                numpoints = cellfun(@numel,Connected.PixelIdxList);
                if sum(overlap)<=3000 %&&  ~any(numpoints>1000) % Check that total number of nans are less than 50 percent
                    T1(index,:) = datcopy(beg:fin);
                    index = index+1;
%                 elseif any(numpoints>1000)
%                     notes{end+1} = ['Trial ',num2str(i),' at T1 excluded. More than one second missing.'];
                end
    
                if sum(overlap)>3000
                    notes{end+1} = ['Trial ',num2str(i),' at T1 excluded. More than 50% missing.'];
    
                end
            end
        end
        T1(index:end,:) = [];
        dat.T1 = T1;
    
        if size(T1,1) < excludethresh
            notes{end+1} = 'More than 30 percent of trials at T1 missing.';
        end
    
    catch
        notes{end+1} = 'Markers for T1 may be missing';
    end
    
    dat.notes = notes;
    
    save(filename_eye(1:end-4),"dat")
    
    clear dat

end


%% Inspect the saved files and the figures and check which ones can be included

filepath = 'Path to the folder where the processed files are saved/';
files = dir(strcat(filepath,'*.mat'));
figs = dir(strcat(filepath,'*.fig'));
figs(1:3) = [];
ketanproc = [];
compareName = "KetanProc";
for i  = 1:length(files)
    tempfile = files(i).name;
    tempfig = figs(i).name;
    if contains(tempfile,compareName)
        ketanproc(i) = true;
    else
        ketanproc(i) = false;
    end
end
files = files(~ketanproc);

excludethresh = 15;
InclusionsPEP = {};
InclusionsRS = {};
includeRS = false;
includeT0 = false;
includeT1 = false;

for i = 1:length(files)

    filename = files(i).name;
    data = load( fullfile(filepath,filename) );
    disp(strcat('You are looking at participant : ',filename))


    if isfield(data.dat,'RS1') && isfield(data.dat,'RS2')
        includeRS = true;
    else
        includeRS = false;
    end

    if isfield(data.dat,'T0') && isfield(data.dat,'T1')
        if size(data.dat.T0,1) >= excludethresh
            includeT0 = true;
        else
            includeT0 = false;
        end

        if size(data.dat.T1,1) >= excludethresh
            includeT1 = true;
        else
            includeT1 = false;
        end

       
    else
        includeT0 = false;
        includeT1 = false;
    end
    

    if includeRS
        InclusionsRS{end+1} = string(filename);
    end

    if includeT1 && includeT0
        InclusionsPEP{end+1} = string(filename);
    end

end


save Nocebo_InclusionsPEP.mat InclusionsPEP
save Nocebo_InclusionsRS.mat InclusionsRS

%% Run the analysis for the included PEP data
addpath('Path to the bayesian analysis toolbox for matlab /klabhub-bayesFactor-c911f2e'); % Add the toolbox to the path
filepath = 'Path to where the processed data is stored/';
files = dir(strcat(filepath,'*.mat'));
ketanproc = [];
compareName = "KetanProc";
for i  = 1:length(files)
    tempfile = files(i).name;
    if contains(tempfile,compareName)
        ketanproc(i) = true;
    else
        ketanproc(i) = false;
    end
end

files = files(~ketanproc);

load('Nocebo_InclusionsPEP.mat', 'InclusionsPEP');
HFSdata = readtable("Path to the HFS data sheet containing group allocations/HFS_Data.xlsx");

T0 = zeros(length(InclusionsPEP),6000);
T1 = zeros(length(InclusionsPEP),6000);
[T0peak, T1peak] = deal(zeros(length(InclusionsPEP),2));
participants = zeros(length(InclusionsPEP),1);
grp = zeros(length(InclusionsPEP),1);

for i = 1:length(InclusionsPEP)
    for j = 1:length(files)

        if InclusionsPEP{i} == string(files(j).name)
            filename = char(InclusionsPEP{i});
            data = load( fullfile(filepath,filename) );
            nums = filename(isnumber(filename));
            participant = str2double(nums);
            grp(i) = HFSdata.group(HFSdata.participant == participant);
            participants(i) = participant;

            % Calculate ERPs and baseline correct
            T0erp = mean(data.dat.T0);
            T0erp = T0erp - mean(T0erp(1:1000));
            T0(i,:) = T0erp;

            T1erp = mean(data.dat.T1);
            T1erp = T1erp - mean(T1erp(1:1000));
            T1(i,:) = T1erp;

            T0base = mean(T0erp(1:1000));
            T0max = max(T0erp(1000:end));
            T0peak(i,1) = T0base;
            T0peak(i,2) = T0max;

            T1base = mean(T1erp(1:1000));
            T1max = max(T1erp(1000:end));
            T1peak(i,1) = T1base;
            T1peak(i,2) = T1max;

        end
    end
end

srate = 1000;
timevec = linspace(-1,5,srate*6);
% Test for T0 vs T1
T1long = T1(:,1001:end);
T0long = T0(:,1001:end);
timelong = timevec(1001:end);
[clusters, p_values, t_sums, permutation_distribution ] = permutest( T1long', T0long', true, .05, 10000, false);
sigtimes = timelong(clusters{1});

windowmeanT0 = mean(T0long(:,clusters{1}),2);
windowmeanT1 = mean(T1long(:,clusters{1}),2);

[bf10,p] = bf.ttest(windowmeanT1,windowmeanT0,'tail','right','scale', sqrt(2)/2);

%% Define line widths and colors
lineWidth = 2;
colorT0 = 'k';   % Black color for T0
colorT1 = 'r';   % Red color for T1
colorSig = [0.5, 0.5, 0.5];  % Transparent grey color for significance lines

% Plot T0 and T1
figure; set(gcf,'color','w');
plot(timevec, mean(T0), 'LineWidth', lineWidth, 'Color', colorT0)
hold on
plot(timevec, mean(T1), 'LineWidth', lineWidth, 'Color', colorT1)
ylim([-50 250])

% Add black dashed lines at 0 on the x and y axes
plot(get(gca, 'xlim'), [0, 0], 'k--', 'LineWidth', 1.5)
plot([0, 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5)

% Highlight significant times
x = [sigtimes(1), sigtimes(1), sigtimes(end), sigtimes(end)];
y = [get(gca, 'ylim'), fliplr(get(gca, 'ylim'))];
patch(x, y, colorSig, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

% Legend and labels
legend('T0', 'T1', '', '','Non-sig cluster', 'FontSize', 10) % Adjust legend font size
xlabel('Time (s)', 'FontSize', 14) % Adjust x-axis label font size
ylabel('Pupil Diameter (a.u)', 'FontSize', 14) % Adjust y-axis label font size
title('Pinprick Evoked Pupil Responses', 'FontSize', 16) % Adjust title font size

% Adjust figure properties
grid on
set(gca, 'FontName', 'Arial', 'FontSize', 20) % Adjust tick label font size
box off


%% With standard error patch

% Define colors and line width
lineWidth = 2;
colorT1 = [0.85, 0.33, 0.1]; % Adjust color as needed
colorT0 = [0, 0, 0]; % Adjust color as needed
% figSize = [100, 100, 800, 800]; % [left, bottom, width, height]
colorSig = [0.4, 0.6, 0.8];  % Soft blue with 50% transparency

% Calculate standard errors
SE_T0 = std(T0) / sqrt(size(T0, 1));
SE_T1 = std(T1) / sqrt(size(T1, 1));

% Plot T0 and T1 with standard error patches
figure('Color', 'w'); set(gcf,'color','w');
hold on;
fill([timevec, fliplr(timevec)], [mean(T0) + SE_T0, fliplr(mean(T0) - SE_T0)], colorT0, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
fill([timevec, fliplr(timevec)], [mean(T1) + SE_T1, fliplr(mean(T1) - SE_T1)], colorT1, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(timevec, mean(T0), 'LineWidth', lineWidth, 'Color', colorT0);
plot(timevec, mean(T1), 'LineWidth', lineWidth, 'Color', colorT1);

% Add black dashed lines at 0 on the x and y axes
plot(get(gca, 'xlim'), [0, 0], 'k--', 'LineWidth', 1.5);
plot([0, 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);

% Highlight significant times
x = [sigtimes(1), sigtimes(1), sigtimes(end), sigtimes(end)];
y = [get(gca, 'ylim'), fliplr(get(gca, 'ylim'))];
patch(x, y, colorSig, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Legend and labels
legend('','','T0', 'T1','','','Non-Sig Cluster', 'Location', 'best');
xlabel('Time (s)');
ylabel('Pupil Diameter (a.u)');
title('Pinprick Evoked Pupil Responses');

% Adjust figure properties
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 12);
box off;

%%

% Extract critical t value
critical_t = tinv(1 - .05, size(T0,1)-1);
ts = zeros(size(T0,2),1);
for i = 1:size(T0,2)
    [h, p, ci, stats] = ttest(T1(:,i), T0(:,i));
    ts(i) = stats.tstat;
end
siginds = ts>critical_t;
Connected = bwconncomp(siginds);
siginds = Connected.PixelIdxList{2};
sigtimes = timevec(siginds);

sigBFs = zeros(size(siginds));
for ind = 1:length(siginds)

    [sigBFs(ind),p] = bf.ttest('T',ts(siginds(ind)),'N',size(T1,1));

end
plot(sigBFs)

%%
% Define line widths and colors
lineWidth = 2;
colorObserved = 'k';    % Black color for observed t values
colorCritical = 'r';    % Red color for critical T (one sided)
colorT0 = 'k';          % Black color for T0
colorT1 = 'r';          % Red color for T1

% Font sizes
titleFontSize = 14;
axisLabelFontSize = 12;
legendFontSize = 10;

% Plot 1: Observed t values and Critical T (one sided)
figure; set(gcf,'color','w');
% subplot(211)
plot(timevec, ts, 'LineWidth', lineWidth, 'Color', colorObserved)
ylim([-5 5])
hold on
yline(critical_t, '--', 'Color', colorCritical, 'LineWidth', lineWidth);
plot([sigtimes(1); sigtimes(1)], get(gca, 'ylim'), 'r-', 'LineWidth', lineWidth)
plot([sigtimes(end); sigtimes(end)], get(gca, 'ylim'), 'r-', 'LineWidth', lineWidth)
yline(0, '--k');
xline(0, '--k');
legend('Observed t values', 'Critical T (one sided)', '', '', 'FontSize', legendFontSize)
grid on
title('A', 'FontSize', titleFontSize)
xlabel('Time (s)', 'FontSize', axisLabelFontSize)
ylabel('t Values', 'FontSize', axisLabelFontSize)
grid on
set(gca, 'FontName', 'Arial', 'FontSize', 20) % Adjust tick label font size
box off

% Plot 2: T0 and T1
subplot(212)
plot(timevec, mean(T0), 'LineWidth', lineWidth, 'Color', colorT0)
hold on
plot(timevec, mean(T1), 'LineWidth', lineWidth, 'Color', colorT1)
ylim([-50 250])
yline(0, '--k');
xline(0, '--k');
plot([sigtimes(1); sigtimes(1)], get(gca, 'ylim'), 'r-', 'LineWidth', lineWidth)
plot([sigtimes(end); sigtimes(end)], get(gca, 'ylim'), 'r-', 'LineWidth', lineWidth)
legend('T0', 'T1', '', '', 'FontSize', legendFontSize)
grid on
title('B', 'FontSize', titleFontSize)
xlabel('Time (s)', 'FontSize', axisLabelFontSize)
ylabel('Pupil Diameter (a.u.)', 'FontSize', axisLabelFontSize)

grid on
set(gca, 'FontName', 'Arial', 'FontSize', 20) % Adjust tick label font size
box off

%%
% Test for Control vs Nocebo
noceboT0 = T0(grp == 0,:);
noceboT1 = T1(grp == 0,:);
controlT0 = T0(grp == 1,:);
controlT1 = T1(grp == 1,:);

nocebo = noceboT1 - noceboT0;
control = controlT1 - controlT0;

nocebolong = nocebo(:,1001:end);
controllong = control(:,1001:end);
[clusters, p_values, t_sums, permutation_distribution ] = permutest( nocebolong', controllong', false, .05, 10000, false);
sigtimes = timelong(clusters{1});



%% Plots for each group 
% Define the figure properties
% figSize = [100, 100, 800, 800]; % Large square figure
figure('Color', 'w');

% Define line widths, colors, and data for the Control group
lineWidth = 2;
colorT0 = 'k'; % Black for T0
colorT1 = 'r'; % Red for T1
SE_T0_control = std(controlT0) / sqrt(size(controlT0, 1));
SE_T1_control = std(controlT1) / sqrt(size(controlT1, 1));

% Plot Control group data
subplot(2,1,1)
hold on;
fill([timevec, fliplr(timevec)], [mean(controlT0) + SE_T0_control, fliplr(mean(controlT0) - SE_T0_control)], colorT0, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
fill([timevec, fliplr(timevec)], [mean(controlT1) + SE_T1_control, fliplr(mean(controlT1) - SE_T1_control)], colorT1, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(timevec, mean(controlT0), 'LineWidth', lineWidth, 'Color', colorT0);
plot(timevec, mean(controlT1), 'LineWidth', lineWidth, 'Color', colorT1);
ylim([-50 250]);
title('Pupil Responses: Control Group', 'FontSize', 16);

% Define data for the Nocebo group
SE_T0_nocebo = std(noceboT0) / sqrt(size(noceboT0, 1));
SE_T1_nocebo = std(noceboT1) / sqrt(size(noceboT1, 1));

% Plot Nocebo group data
subplot(2,1,2)
hold on;
fill([timevec, fliplr(timevec)], [mean(noceboT0) + SE_T0_nocebo, fliplr(mean(noceboT0) - SE_T0_nocebo)], colorT0, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
fill([timevec, fliplr(timevec)], [mean(noceboT1) + SE_T1_nocebo, fliplr(mean(noceboT1) - SE_T1_nocebo)], colorT1, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(timevec, mean(noceboT0), 'LineWidth', lineWidth, 'Color', colorT0);
plot(timevec, mean(noceboT1), 'LineWidth', lineWidth, 'Color', colorT1);
ylim([-50 250]);
title('Pupil Responses: Nocebo Group', 'FontSize', 16);

% Shared Figure Adjustments
for i = 1:2
    subplot(2,1,i);
    plot(get(gca, 'xlim'), [0, 0], 'k--', 'LineWidth', 1.5);
    plot([0, 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Pupil Diameter (a.u)', 'FontSize', 14);
    grid on;
    legend('','','T0', 'T1', 'Location', 'best');
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    box off;
end

%% Now plot for the differences

% Figure for the differences
% Define line widths, colors, and data for standard errors
lineWidth = 2;
colorControl = 'k';   % Black color for Control
colorNocebo = 'r';   % Red color for Nocebo
colorSig = [0.4, 0.6, 0.8];  % Soft blue with 50% transparency
SE_Control = std(control) / sqrt(size(control, 1));
SE_Nocebo = std(nocebo) / sqrt(size(nocebo, 1));


% Plot Control and Nocebo differences with error patches
figure; set(gcf, 'color', 'w'); % Large square figure
hold on;
fill([timevec, fliplr(timevec)], [mean(control) + SE_Control, fliplr(mean(control) - SE_Control)], colorControl, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
fill([timevec, fliplr(timevec)], [mean(nocebo) + SE_Nocebo, fliplr(mean(nocebo) - SE_Nocebo)], colorNocebo, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
plot(timevec, mean(control), 'LineWidth', lineWidth, 'Color', colorControl);
plot(timevec, mean(nocebo), 'LineWidth', lineWidth, 'Color', colorNocebo);
ylim([-70 70]);

% Add black dashed lines at 0 on the x and y axes
plot(get(gca, 'xlim'), [0, 0], 'k--', 'LineWidth', 1.5);
plot([0, 0], get(gca, 'ylim'), 'k--', 'LineWidth', 1.5);

% Highlight significant times (if any)
x = [sigtimes(1), sigtimes(1), sigtimes(end), sigtimes(end)];
y = [get(gca, 'ylim'), fliplr(get(gca, 'ylim'))];
patch(x, y, colorSig, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Legend and labels
legend('','','Control', 'Nocebo','','','Non-Sig Cluster', 'Location', 'best', 'FontSize', 10);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Pupil Diameter Differences (a.u)', 'FontSize', 14);
title('Difference Waves (T1 - T0)', 'FontSize', 16);

% Adjust figure properties
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 20);
box off;
%%

% Extract critical t value
ts = zeros(size(control,2),1);
ps = zeros(size(control,2),1);
for i = 1:size(T0,2)
    [h, p, ci, stats] = ttest2(control(:,i), nocebo(:,i));
    ts(i) = stats.tstat;
    ps(i) = p;
end
critical_t = tinv(1 - .05, stats.df);
siginds = ts>critical_t;
Connected = bwconncomp(siginds);
siginds = Connected.PixelIdxList{1};
sigtimes = timevec(siginds);

figure; subplot(211)
plot(timevec,mean(control),'k-',LineWidth=2)
hold on
plot(timevec,mean(nocebo),'r-',LineWidth=2)
ylim([-50 50])
yline(0, '--k');
xline(0, '--k');
plot([sigtimes(1);sigtimes(1)],get(gca,'ylim'),'r-')
plot([sigtimes(end);sigtimes(end)],get(gca,'ylim'),'r-')
legend('Control','Nocebo')

subplot(212)
plot(timevec,ts,'k-',LineWidth=2)
hold on
yline(critical_t, '--r', LineWidth=2);
ylim([-5 5])
% yline(-critical_t, '--r', LineWidth=2);
plot([sigtimes(1);sigtimes(1)],get(gca,'ylim'),'r-')
plot([sigtimes(end);sigtimes(end)],get(gca,'ylim'),'r-')
legend('Observed t values','Critical T (one sided)')



% Analyse the peak data
all_vals = [T0peak;T1peak];
all_vals = array2table(all_vals,"VariableNames",{'Baseline','Max'});
all_vals.Participant = [participants;participants];
all_vals.Time = [repmat("T0",size(T0peak,1),1); repmat("T1",size(T1peak,1),1)];
all_vals.Group = [grp;grp];
all_vals = stack(all_vals,{'Baseline','Max'},'IndexVariableName','Value');

filename = 'Peak_vals.xlsx';
writetable(all_vals, filename);

%% Run the analysis for the included RS data

load('Nocebo_InclusionsRS.mat', 'InclusionsRS');
filepath = 'Path to where the processed pupil data is stored/';
files = dir(strcat(filepath,'*.mat'));
ketanproc = [];
compareName = "KetanProc";
for i  = 1:length(files)
    tempfile = files(i).name;
    if contains(tempfile,compareName)
        ketanproc(i) = true;
    else
        ketanproc(i) = false;
    end
end
files = files(~ketanproc);

HFSdata = readtable("Path to HFS data sheet containing group allocations/HFS_Data.xlsx");

RS1 = zeros(length(InclusionsRS),1);
RS2 = zeros(length(InclusionsRS),1);
grp = [];
participants = zeros(length(InclusionsRS),1);

for i = 1:length(InclusionsRS)
    for j = 1:length(files)

        if InclusionsRS{i} == string(files(j).name)
            filename = char(InclusionsRS{i});
            data = load( fullfile(filepath,filename) );
            nums = filename(isnumber(filename));
            participant = str2double(nums);
            participants(i) = participant;
            grp(i) = HFSdata.group(HFSdata.participant == participant);

            % Calculate RS1 and RS2 means
            RS1(i) = mean(data.dat.RS1);
            RS2(i) = mean(data.dat.RS2);

        end
    end
end

RS = [RS1,RS2];
RS = array2table(RS,"VariableNames",{'RS1','RS2'});
RS.Group = grp';
RS.Participant = participants;


filename = 'RS_vals.xlsx';
writetable(RS, filename);
