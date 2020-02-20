function CalibrationInfo=DAG_get_reward_calibration_info

CalibrationInfo=struct();
current_user = getUserName; 
calibration_file = ['C:\Users' filesep current_user filesep 'Dropbox' filesep 'DAG' filesep 'DAG_toolbox' filesep 'Reward_system_scripts' filesep 'calibration.xlsx'];
[~,SHEETS] = xlsfinfo(calibration_file);

counter=1;
for sheetindex=1:numel(SHEETS)
    sheetname=SHEETS{sheetindex};
    seperator_idx=findstr(sheetname,'_');
    if numel(seperator_idx)<3
        continue
    end
    sheetsetup   =str2mat(sheetname(seperator_idx(1)+1:seperator_idx(2)-1));
    if isstr(sheetsetup) && strcmp(sheetsetup,'UMG')
        sheetsetup=-1;
    else        
        sheetsetup   =str2num(sheetsetup);
    end
    
    sheetdate =str2num(sheetname(1:seperator_idx(1)-1));
    sheetvoltage =sheetname(seperator_idx(2)+1:seperator_idx(end)-1);
    sheetvoltage(findstr(sheetvoltage,'_'))='.';
    sheetvoltage=str2num(sheetvoltage);
            
    CalibrationInfo(counter).ms_cal = xlsread(calibration_file, sheetname, 'B3:P3');
    CalibrationInfo(counter).ml_cal = xlsread(calibration_file, sheetname, 'B7:P7');
    
    CalibrationInfo(counter).date=sheetdate;
    %CalibrationInfo(counter).sheetname=sheetname;
    CalibrationInfo(counter).setup=sheetsetup;
    CalibrationInfo(counter).voltage=sheetvoltage;
    counter=counter+1;
end
