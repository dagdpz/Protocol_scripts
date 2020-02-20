function DAG_protocol_update(monkey,dates)
%DAG_protocol_update('Cornelius',[20170919 20170919]);


%%
switch monkey
    case {'Curius'}
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};
        
        reward_voltage=[11.3 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
            case {'Curius_phys'}
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};
        
        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
                    case {'Linus_phys'}
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};
        
        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
    case 'Cornelius'
%         mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
%             'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
%             'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
%             'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
%             'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
%             'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
%             'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
%             'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
%             'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
%             'microstim_fraction','microstim_start', 'microstim_state',  ...
%             'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
%             'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};
%      mastertable(1,:)={'Session','Run','N_Trials','hits','hits_per_total','completed','hits_per_completed','Type','Effector','Setup',...
%             'initiated','hits_per_initiated','Reward_time','Reward_ml','start_run','end_run','choice_fraction',...            
%             'ST_successrate','SED_successrate','SHD_successrate','DT_successrate','TED_successrate','THD_successrate','DED_successrate','DHD_successrate',...
%             'left_chosen_base','Sensor_L','Sensor_R','x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
%             'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
%             'fix_hold','fix_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
%             'cue_hold','cue_hold_var','mem_hold','mem_hold_var','ITI_success','ITI_fail'};
        
         mastertable(1,:)={'Session','Run','N_Trials','hits','hits_per_total','completed','hits_per_completed','Type','Effector','Setup',...
            'initiated','hits_per_initiated','Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','Sensor_L','Sensor_R','x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'cue_hold','cue_hold_var','mem_hold','mem_hold_var','ITI_success','ITI_fail'};

        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
     case 'Cornelius_phys'
         mastertable(1,:)={'Session','Run','N_Trials','hits','hits_per_total','completed','hits_per_completed','Type','Effector','Setup',...
            'initiated','hits_per_initiated','Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','Sensor_L','Sensor_R','x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'cue_hold','cue_hold_var','mem_hold','mem_hold_var','ITI_success','ITI_fail'};

        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
         
     case 'Cornelius_ina'
         mastertable(1,:)={'Session','Run','N_Trials','hits','hits_per_total','completed','hits_per_completed','Type','Effector','Setup',...
            'initiated','hits_per_initiated','Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','Sensor_L','Sensor_R','x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'cue_hold','cue_hold_var','mem_hold','mem_hold_var','ITI_success','ITI_fail'};

        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
    case 'Bacchus'
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};

        reward_voltage=[11.3 11.3 11.3 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
    case {'Linus','Linus_ina'}
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};

        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
        
    case {'Flaffus','Flaffus_ina', 'Flaffus_phys', 'Tesla_phys'}
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad','eye_fix_siz','eye_tar_siz','hnd_fix_siz','hnd_tar_siz',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};

        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
    case 'Tesla'
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};

        reward_voltage=[6.5 11.3 11.3 11.3 11.3]; % setup1 setup2 setup3 setup4 else
    case 'Magnus'
        mastertable(1,:)={'Session','Run','N_Trials','initiated','completed','hits','hits_per_total','hits_per_initiated','hits_per_completed','Type','Effector','Setup',...
            'Reward_time','Reward_ml','start_run','end_run','choice_fraction',...
            'left_chosen_base','left_chosen_stimulated', 'reward_modulation_fraction','high_reward_chosen_fraction',...
            'Sensor_L','Sensor_R','demanded_hand_fraction','reach_hand_fraction','R_hand_chosen_fraction',...
            'x_positions','y_positions','eye_fix_rad','eye_tar_rad','hnd_fix_rad','hnd_tar_rad',...
            'n_targets', 'cue_x_pos','cue_y_pos','convexities','convex_sides','targets_insp',...
            'fix_acq_eye','fix_acq_hnd','tar_acq_eye','tar_acq_hnd','tar_inv_acq_eye','tar_inv_acq_hnd',...
            'fix_hold','fix_hold_var','cue_hold','cue_hold_var','mem_hold','mem_hold_var','del_hold','del_hold_var','tar_inv_hold','tar_inv_hold_var','tar_hold','tar_hold_var',...
            'ITI_success','ITI_success_var','ITI_fail','ITI_fail_var','grace_time_eye','grace_time_hand',...
            'microstim_fraction','microstim_start', 'microstim_state',  ...
            'Eye_RT_R_base_mean', 'Eye_RT_R_base_std', 'Eye_RT_L_base_mean', 'Eye_RT_L_base_std','Hnd_RT_R_base_mean', 'Hnd_RT_R_base_std','Hnd_RT_L_base_mean', 'Hnd_RT_L_base_std',...
            'Eye_RT_R_stim_mean', 'Eye_RT_R_stim_std', 'Eye_RT_L_stim_mean', 'Eye_RT_L_stim_std','Hnd_RT_R_stim_mean', 'Hnd_RT_R_stim_std','Hnd_RT_L_stim_mean', 'Hnd_RT_L_stim_std'};
        
        reward_voltage=[11.6 11.3 12 11.3 11.3]; % setup1 setup2 setup3 setup4 else
end




%%
start_date=dates(1);
while true
    month=floor((start_date-floor(start_date/10000)*10000)/100);
    if month == 12
        end_date=(floor(start_date/10000)+1)*10000+100;
    else
        end_date=(floor(start_date/100)+1)*100;
    end
    if end_date>dates(2)
        DAG_protocol_from_trials(monkey,[start_date dates(2)], mastertable, reward_voltage);
        break
    end
    DAG_protocol_from_trials(monkey,[start_date end_date], mastertable, reward_voltage);
    start_date=end_date;
end