function ml=DAG_get_ml_from_CalibrationInfo(CalibrationInfo,valve_opening_time,Session,setup,voltage)
if setup==4
    setup=2;
end


CalibrationDate=[CalibrationInfo.date];
CalibrationSetup=[CalibrationInfo.setup];
CalibrationVoltage=[CalibrationInfo.voltage];
SetupVoltageIndex=CalibrationSetup==setup & CalibrationVoltage==voltage;
ValidCalibrationDates=CalibrationDate(SetupVoltageIndex);

idx=find(SetupVoltageIndex);
if all(ValidCalibrationDates>Session)    
    CalibrationIndex=idx(1);
elseif all(ValidCalibrationDates<=Session)
    CalibrationIndex=idx(end);    
else
    CalibrationIndex=idx(diff(ValidCalibrationDates>Session)==1);    
end

ms_cal=CalibrationInfo(CalibrationIndex).ms_cal;
ml_cal=CalibrationInfo(CalibrationIndex).ml_cal;
p = polyfit(ms_cal,ml_cal,1); % linear calibration
ms = valve_opening_time*1000;
ml = p(1)*ms + p(2);
end
