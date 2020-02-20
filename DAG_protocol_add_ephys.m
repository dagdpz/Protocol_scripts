function DAG_protocol_add_ephys(monkey)
%monkey='Linus_phys';

dag_drive_IP=get_dag_drive_IP;
protocol_folder=[dag_drive_IP, 'Protocols' filesep monkey];
protocol_xls_file= most_recent_version(protocol_folder,[monkey '_protocol.xls']);

[~, ~, protocol_runs_table]=xlsread(protocol_xls_file,'Runs');
[~, ~, protocol_sessions_table]=xlsread(protocol_xls_file,'Sessions');

dropboxpath=['C:\Users\' getUserName '\Dropbox' filesep 'DAG' filesep 'phys' filesep monkey '_dpz'];
[~, sheets_available]=xlsfinfo([dropboxpath filesep monkey(1:3) '_sorted_neurons.xlsx']);
if ismember('final_sorting',sheets_available)
    [~, ~, sorting_table]=xlsread([dropboxpath filesep monkey(1:3) '_sorted_neurons.xlsx'],'final_sorting');
elseif ismember('automatic_sorting',sheets_available)
    [~, ~, sorting_table]=xlsread([dropboxpath filesep monkey(1:3) '_sorted_neurons.xlsx'],'automatic_sorting');
else
    sorting_table={'Monkey','Session','Date','Run','Block','Chan','z','Unit','N_spk','Neuron_ID','Times_same_unit','Site_ID'};
end

sort_session_idx=find_column_index(sorting_table,'Date');
sort_run_idx=find_column_index(sorting_table,'Run');
sort_neuron_idx=find_column_index(sorting_table,'Neuron_ID');
sort_site_idx=find_column_index(sorting_table,'Site_ID');
sort_target_idx=find_column_index(sorting_table,'Target');

for x=1:size(protocol_runs_table,2)
fname=strrep(protocol_runs_table{1,x},' ','_');
fname=strrep(fname,'?','');
fname=strrep(fname,'@','');
fname=strrep(fname,'!','');
fname=strrep(fname,'(','');
fname=strrep(fname,')','');
if ~ischar(fname)
    continue;
end
idx.(fname)=x;
end

sort_table_extension={'Effector','Type','choice_fraction','demanded_hand_fraction'};
%sort_table_extension={'Effector','Type','choice_fraction','demanded_hand_fraction','x_positions','y_positions'};
for r=2:size(sorting_table,1)
    
    session=sorting_table{r,sort_session_idx};
    run=sorting_table{r,sort_run_idx};
    
    current_entry=[false [protocol_runs_table{2:end,idx.Session}]==session & [protocol_runs_table{2:end,idx.Run}]==run];
    
    for x=1:size(sort_table_extension,2)
        if any(current_entry)
            fname=sort_table_extension{1,x};
            sort_table_extension(r,x)=protocol_runs_table(current_entry,idx.(fname));
        else
            sort_table_extension(r,x)={'-1 -1'};
        end
    end
    
%     current_entries=[false [sorting_table{2:end,idx.Date}]==session & [sorting_table{2:end,idx.Run}]==run];
%     Sites=unique(sorting_table(current_entries,idx.Site_ID));
%     Cells=unique(sorting_table(current_entries,idx.Neuron_ID));    
%     Cells(~cellfun(@isempty,strfind(Cells,'_00')))=[]; % rmoving "no cell", only LFP site entries
%     
%     ephys(r).Sites=Sites;
%     ephys(r).Cells=Cells;
    
end

for x=1:size(sort_table_extension,2)
    col_idx=find(~cellfun(@isempty,strfind(sorting_table(1,:),sort_table_extension{1,x})));
    if isempty(col_idx)
    col_idx=size(sorting_table,2)+1;
    end
    sorting_table(:,col_idx)=sort_table_extension(:,x);
end

%Adjusting in order to have it the way we want
    idx_1=find_column_index(sort_table_extension,'choice_fraction');sort_table_extension{1,idx_1}='choices_present';
    idx_2=find_column_index(sort_table_extension,'demanded_hand_fraction');sort_table_extension{1,idx_2}='hands_interleaved';
for x=2:size(sort_table_extension,1)
    sort_table_extension{x,idx_1}=double(any(sort_table_extension{x,idx_1}));
    temp_hands=str2num(sort_table_extension{x,idx_2});
    sort_table_extension{x,idx_2}=double(temp_hands(1) && temp_hands(2));    
end
% and format everthing to string
for x=1:numel(sort_table_extension)
    if ~ischar(sort_table_extension{x})
        sort_table_extension{x}=num2str(sort_table_extension{x});
    end
end
sort_table_extension(:,end+1)=sorting_table(:,sort_target_idx);



for x=1:size(sort_table_extension,2)
    col_idx=find(~cellfun(@isempty,strfind(sorting_table(1,:),sort_table_extension{1,x})));
    if isempty(col_idx)
    col_idx=size(sorting_table,2)+1;
    end
    sorting_table(:,col_idx)=sort_table_extension(:,x);
end
[unique_rows, ~ , idx_st]=uniqueRowsCA(sort_table_extension(2:end,:));
sort_table_conditions=[sort_table_extension(1,:); unique_rows];


% now, find number of cells and sites for each of the conditions... :)
sort_table_conditions{1,end+1}='Sessions';
sort_table_conditions{1,end+1}='Sites';
sort_table_conditions{1,end+1}='Neurons';
for r=2:size(sort_table_conditions,1)
    current_entries=[false; idx_st==r-1]; % the first row qith titles makes this a bit messy
    
    Sessions=unique(cellfun(@(x) num2str(x),sorting_table(current_entries,sort_session_idx),'UniformOutput',false));
    Sites=unique(sorting_table(current_entries,sort_site_idx));
    Cells=unique(sorting_table(current_entries,sort_neuron_idx));    
    Cells(~cellfun(@isempty,strfind(Cells,'_00')))=[]; % removing "no cell", only LFP site entries


sort_table_conditions{r,end-2}=numel(Sessions);    
sort_table_conditions{r,end-1}=numel(Sites);
sort_table_conditions{r,end}=numel(Cells);
    
end

xlswrite([dropboxpath filesep monkey(1:3) '_datasets_overview.xlsx'],sort_table_conditions,'summary');
end