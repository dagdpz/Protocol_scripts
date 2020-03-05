function a = DAG_clean_protocol(monkey,dates)
a=[];
monkey='Cornelius';
dates=[20130901 20160226];
drive=DAG_get_server_IP;
switch monkey
    case 'Linus'
        data_path= strcat(drive,':', filesep,'Data', filesep, monkey, '_microstim_with_parameters', filesep);
        prot_path=[drive,':', filesep,'microstim_behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
    case 'Curius'
        data_path= strcat(drive,':', filesep,'Data', filesep, monkey, '_microstim_with_parameters', filesep);
        prot_path=[drive,':', filesep,'microstim_behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
    case 'Cornelius'
        data_path= strcat(drive, 'Data', filesep, monkey, filesep);
        prot_path=[drive, 'Protocols', filesep, monkey, filesep, monkey, '_protocol.xls'];
        %data=[pwd, filesep, monkey, '_protocol.xls'];
end



[num,masterstring_orig,RAW_orig] = xlsread(prot_path,'mastertable');
cell_from_protocol=num(:,1:2);


session_folders_dir=dir(data_path);
session_folders={session_folders_dir([session_folders_dir.isdir]).name};
valid_indexes=cellfun(@(x) ~isempty(str2double(x)),session_folders);
session_folders=session_folders(valid_indexes);
valid_indexes=cellfun(@(x) str2double(x)>=dates(1) && str2double(x)<=dates(2),session_folders);
session_folders=session_folders(valid_indexes);
clear cell_from_data
idx=0;
for f=1:numel(session_folders)
    sub_dir=dir([data_path filesep session_folders{f} filesep '*.mat']);
    files_in_folder={sub_dir.name};
    for file=1:numel(files_in_folder)
        idx=idx+1;
        cell_from_data(idx,1)=str2num(session_folders{f});
        cell_from_data(idx,2)=str2num((files_in_folder{file}(end-5:end-4)));        
        idx_valid=ismember(cell_from_protocol(:,1),cell_from_data(idx,1)) & ismember(cell_from_protocol(:,2),cell_from_data(idx,2)) ;
        if sum(idx_valid)<1
            disp(['No protocol entry for date ' num2str(cell_from_data(idx,1)) ' run ' num2str(cell_from_data(idx,2))])
        elseif sum(idx_valid)>1
            disp(['More than one protocol entry for date ' num2str(cell_from_data(idx,1)) ' run ' num2str(cell_from_data(idx,2))])            
        end
        if idx_valid~=0
            
        end
    end
end
