function LAVES_protocol_update(monkey, dates, excel_file, protocol, entry_person)
% monkey = 'Cornelius';
% dates = [20170329; 20170331];
% excel_file = 'Y:\Protocols\Cornelius\Cornelius_protocol.xls';
% protocol = 'Y:\Personal\Uwe\MATLAB\uz\cor_test_protocol.txt';
% entry_person = 'uwz';

sheet = 'CNL protocol';
%protocol_new = 'Y:\Personal\Uwe\MATLAB\uz\cor_test_protocol_new.txt';

% TXT_PROTOCOL_UPDATE
% This function updates the txt-file of the LAVES protocol by the specified
% dates using the following inputs:

% This function will search the specified excel file for the following 
% columns to find the data necessary for protocol update:
column_session = 'Session';
column_start_chair_time = 'chair in';
column_end_chair_time = 'chair out';
column_method = 'method';
column_method_comment = 'method_comments';
column_fluid_reward = 'Reward [ml]';
column_fluid_reward_code = 'Reward code';
column_non_fluid_reward = 'Fruits in Experiment [ml]';
column_additional_fluid = 'additional water [ml]';
column_additional_non_fluid = 'Fruits [ml]';
column_weight = 'Weight';
column_aim = 'Task';
column_daily_trials = 'initiated';
column_hitrate = 'hits_per_initiated';

% open file for reading

fID = fopen(protocol,'a+');
[num,txt,data_cell] = xlsread(excel_file,sheet);

% identify the dates to update
all_dates = dates_from_range(dates, 'yyyyMMdd', 'yyyyMMdd', 'type_double');

% For each date, search for the data in the excel sheet
% find the column indeces of the column in data_cell
column_idx_sessions  = DAG_find_column_index(data_cell,column_session);
column_idx_start_chair_time  = DAG_find_column_index(data_cell,column_start_chair_time);
column_idx_end_chair_time  = DAG_find_column_index(data_cell,column_end_chair_time);
column_idx_method  = DAG_find_column_index(data_cell,column_method);
column_idx_method_comment  = DAG_find_column_index(data_cell,column_method_comment);
column_idx_fluid_reward  = DAG_find_column_index(data_cell,column_fluid_reward);
column_idx_fluid_reward_code = DAG_find_column_index(data_cell,column_fluid_reward_code);
column_idx_non_fluid_reward  = DAG_find_column_index(data_cell,column_non_fluid_reward);
column_idx_additional_fluid  = DAG_find_column_index(data_cell,column_additional_fluid);
column_idx_additional_non_fluid  = DAG_find_column_index(data_cell,column_additional_non_fluid);
column_idx_weight  = DAG_find_column_index(data_cell,column_weight);
column_idx_aim  = DAG_find_column_index(data_cell,column_aim);
column_idx_daily_trials  = DAG_find_column_index(data_cell,column_daily_trials);
column_idx_hitrate  = DAG_find_column_index(data_cell,column_hitrate);

monkey_upper_case = monkey(1:3);
monkey_lower_case = lower(monkey_upper_case);

for dat = 1:numel(all_dates)
    session = all_dates(dat);
    session_dt = datetime(num2str(session),'InputFormat','yyyyMMdd','format','yyyyMMdd');
    session_minus = datestr(datetime(num2str(session),'InputFormat','yyyyMMdd','format','yyyy-MM-dd'),'yyyy-mm-dd');
    session_str = datestr(session_dt,'yymmdd'); 
    
    % find row index of current session in data_cell
    session_column = {data_cell{:,column_idx_sessions}}';
    row_idx = DAG_find_row_index(session_column,session);
    
    % find row index of previous day in data_cell
    prev_day = str2num(datestr(session_dt-1,'yyyymmdd')); %,'InputFormat','yyyyMMdd'
    row_idx_prev_day = DAG_find_row_index(session_column,prev_day);
    
    if isempty(row_idx_prev_day)
        no_session_on_prev_day = true;
        prev_day_fluid = -1;
        sprintf('No session previous to %d was fount so previous day fluid of "-1" was assumed!', session);
    else
        prev_day_fluid = nansum([data_cell{row_idx_prev_day,column_idx_fluid_reward},...
            data_cell{row_idx_prev_day,column_idx_non_fluid_reward},...
            data_cell{row_idx_prev_day,column_idx_additional_fluid},...
            data_cell{row_idx_prev_day,column_idx_additional_non_fluid}]);
    end
    
    % read out data from data_cell (for an entire training/experimental entry)
    if ~isempty(row_idx)
        start_chair_time = datetime(data_cell{row_idx,column_idx_start_chair_time},'ConvertFrom','excel','format','MM/dd/yyyy HH:mm:ss'); % ,'InputFormat','h:mm'
        start_chair_time = datestr(start_chair_time, 'HH:MM');
        
        end_chair_time = datetime(data_cell{row_idx,column_idx_end_chair_time},'ConvertFrom','excel','format','MM/dd/yyyy HH:mm:ss'); % ,'InputFormat','h:mm'
        end_chair_time = datestr(end_chair_time, 'HH:MM');
        
        method = data_cell{row_idx,column_idx_method};
        
        method_comment = data_cell{row_idx,column_idx_method_comment};
        if isnan(method_comment)
            method_comment_out = '""';
        else
            method_comment_out = ['"',method_comment,'"'];
        end
        
        fluid_reward = data_cell{row_idx,column_idx_fluid_reward};
        
        fluid_reward_code = data_cell{row_idx,column_idx_fluid_reward_code};
        if ~isnan(fluid_reward_code) % add fluid reward code if there is one
            fluid_reward_out = ['"',num2str(fluid_reward),' ', fluid_reward_code,'"'];
        else
            fluid_reward_out = ['"',num2str(fluid_reward),'"'];
        end
        
        non_fluid_reward = data_cell{row_idx,column_idx_non_fluid_reward};
        if isnan(non_fluid_reward)
            non_fluid_reward = 0;
        end
        non_fluid_reward_out = ['"',num2str(non_fluid_reward),'"'];
        
        additional_fluid = data_cell{row_idx,column_idx_additional_fluid};
        if isnan(additional_fluid)
            additional_fluid = 0;
        end
        additional_fluid_out = ['"',num2str(additional_fluid),'"'];
        
        additional_non_fluid = data_cell{row_idx,column_idx_additional_non_fluid};
        if isnan(additional_non_fluid)
            additional_non_fluid = 0;
        end
        additional_non_fluidout = ['"',num2str(additional_non_fluid),'"'];
        
        total_daily_fluid = sum([fluid_reward,non_fluid_reward,additional_fluid,additional_non_fluid]);
        
        weight = data_cell{row_idx,column_idx_weight}; % fprintf('%.2f', weight)
        
        aim = ['"',data_cell{row_idx,column_idx_aim},'"'];
        
        daily_trials = data_cell{row_idx,column_idx_daily_trials};
        
        hitrate = data_cell{row_idx,column_idx_hitrate};
        
        data_File = [monkey_upper_case,session_minus,'.mat'];
        
        
        % round chair time
        %     start_chair_time_str = datestr(start_chair_time, 'HH:MM');
        %     start_chair_time_vec = datevec(start_chair_time_str, 'HH:MM');            % Generate ‘datevec’ Date Vectors
        %     start_chair_time_rounded = datevec(datenum([start_chair_time_vec(:,4) [30*(start_chair_time_vec(:,5)<30) + 60*(start_chair_time_vec(:,5)>=30)] ]));    % Rounded Date Vectors (not working)
        
        % create output for txt file
        
        % create output for one entire day
        out.method = {1 session_str '24:00' entry_person monkey_lower_case 1 'method' method method_comment_out};
        format.method = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t\t%s\t\t%s';
        out.prevDayFluid = {1 session_str '24:00' entry_person monkey_lower_case 1 'prevDayFluid' prev_day_fluid};
        format.prevDayFluid = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%d';
        out.dailyWater = {1 session_str '24:00' entry_person monkey_lower_case 2 'dailyWater' '""' total_daily_fluid fluid_reward_out non_fluid_reward_out additional_fluid_out additional_non_fluidout};
        format.dailyWater = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s';
        out.weight = {1 session_str '24:00' entry_person monkey_lower_case 1 'weight' weight};
        format.weight = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t\t%1.1f';
        out.chairStateIn = {1 session_str start_chair_time entry_person monkey_lower_case 1 'chairState' 'in' '""'};
        format.chairStateIn = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%s %s';
        out.aim = {1 session_str '24:00' entry_person monkey_lower_case 1 'aim' aim};
        format.aim = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t\t%s';
        out.dayNumTrials = {1 session_str '24:00' entry_person monkey_lower_case 1 'dayNumTrials' daily_trials hitrate fluid_reward};
        format.dayNumTrials = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%d %d %d';
        out.dataFile = {1 session_str '24:00' entry_person monkey_lower_case 1 'dataFile' data_File};
        format.dataFile = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%s';
        out.chairStateOut = {1 session_str end_chair_time entry_person monkey_lower_case 1 'chairState' 'out' '""'};
        format.chairStateOut = '\n%d\t%s\t%s\t\t\t\t%s\t%s\t%d\t%s\t%s %s';
        
        % print output into txt file
        fields = fieldnames(out);
        for fn = 1:numel(fields)
            fprintf(fID, format.(fields{fn}), out.(fields{fn}){1,:});
        end
        fprintf(fID, '\n');
    end
end
% fID = fopen(protocol,'a+');
fclose(fID);

disp('Done! :)')





%%%%%%% Non-functioning code lines trying to insert updates lines at the
%%%%%%% correct position in the .txt file
% % open file for reading
% 
% fIDr = fopen(protocol,'a+');
% 
% % open file for writing
% 
% fIDw = fopen(protocol_new,'w') ;
% 
% % while end of file has not been reached
% 
% while ( ~feof(fidr) )
% 
%       % read line from reading file
% 
%       str = fgets(fIDr) ;
% 
%       % match line to regular expression to determine if replacement needed
% 
%       match = regexp(str,'(?<=Line )\d(?= is this.)', 'match' ) ;
% 
%       % if line is to be replaced
% 
%       if ( ~isempty(match) )
% 
%           % define replacement line
% 
%           str = ['This is line ',match{1},'.',char(10)] ;
% 
%       end
% 
%       % write line to writing file
% 
%       fwrite(fIDw,str) ;
% 
% end
% 
% C = {1 170404 '17:00' 'uwz' 'cor' 1 'chairState' 'out' '""'};
% formatSpec = '\n%d\t%d\t%s\t\t\t\t%s\t%s\t%d\t%s\t%s %s';
% fprintf(fIDr, formatSpec, C{1,:});
% while ~feof(fid)
%     inData = textscan(fIDr, '%s', 1);
%     while isempty(inData{1})
%             
%             inData = textscan(fid, '%s', 1);
%             if feof(fid) && isempty(inData{1})
%                 exitInd = true;
%                 break;
%             end %if
%             
%         end %while
%         if exitInd
%             break;
%         end %if
%     
%     fileInformation = dir(protocol);
%     fileSize = fileInformation.bytes;
%     positionInFile = ftell(fIDr);
%     
% end
% fclose(fIDr);
% % Search for the inputs
% % Search for the column containing sessions
% % str = sprintf('Control Chart For \t %d water flow \t%d oil flow',5,3);
% % fprintf(fileID, formatSpec, str);
% del = '\t';
% A = fprintf(['a' del 'b']);
% DAG_find_column_index
% DAG_update_mastertable_cell
end

%% Subfunctions
function [ all_dates ] = dates_from_range( dates, input_format, output_format, output_type)
%DATES_FROM_RANGE Returns all dates between a start date and an end date
%   Input: 
%       dates - Array containing start date and end date in format yyyymmdd
%           Example: [20170329; 20170405];
%       input_format - format of the input dates as recognized by the
%           datetime function (type help datetime and check InputFormat)
%           Default: 'yyyyMMdd'. 
%       output_format - format of the output dates as recognized by the
%           datetime function (type help datetime and check InputFormat)
%           Default: 'yyyyMMdd'. 
%       output_type - data type of the output dates. Possible values:
%           'type_double', 'type_datetime'. 
%           Default: 'type_datetime'
%   Output:
%       all_dates - A row vector of all dates as datetime in the defined format and data type;

% Created by Uwe Zimmermann (20170420)
if nargin < 2
    input_format = 'yyyyMMdd';
end
if nargin < 3
    output_format = 'yyyyMMdd';
end
if nargin < 4
    output_type = 'type_datetime';
end

dates = datetime(num2str(dates),'InputFormat',input_format,'format',output_format);
if numel(dates) > 1
    
    start_date = dates(1);
    end_date = dates(2);
    interval = 1 ; % 1 day interval
    all_dates = [start_date:interval:end_date];
    
else
    all_dates = dates;
end

switch output_type
    case 'type_double'
    all_dates = str2num(datestr(all_dates,'yyyymmdd')); 
    case 'type_datetime'
        all_dates = all_dates;
end
end

function column_index=DAG_find_column_index(inputcell,title)
column_index=[];
for m=1:size(inputcell,2)
    if strcmp(inputcell{1,m},title)
        column_index=m;
    end
end
end

function row_index=DAG_find_row_index(inputcell,title)
row_index=[];
if ischar(title) % search for a string
for m=1:size(inputcell,1)
    if strcmp(inputcell{m,1},title)
        row_index=m;
    end
end
end
if isnumeric(title) % search for a number
for m=1:size(inputcell,1)
    if isequal(inputcell{m,1},title)
        row_index=m;
    end
end
end
end

