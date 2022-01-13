clear;close all; 
% run each datasets separately 
fileidx = [36]
bleachtime = 24; % in frame number
bleach_var = 0.1; % in frame number
for fileiter = 1:length(fileidx)
    fileidxeach = fileidx(fileiter);
getfilename = sprintf('dwell_%d_size_group_*_Q.mat',fileidxeach);
getfilename = dir(getfilename);
filenames = getfilename;
    %%
    Tracks_cell = cell(1,1);
    training_cell = cell(1,1);
    for fileiter = 1:length(filenames)
      file = filenames(fileiter);
      ld=importdata(file.name);
      a = ld.Tracks_pred;
      Tracks_cell{fileiter,1}= a;
      training_cell{fileiter,1}=ld.TrainingFinal;
    end
     Tracks_cell = vertcat(Tracks_cell{:});
      training_cell = vertcat(training_cell{:});
      %for iter = 1:length(Tracks_cell)
   
          Tracks_cell2 = Tracks_cell;
          training_cell2 = training_cell;
          if ~isempty(Tracks_cell2)
              
      
          Trackmate_Analysis_ML_FINAL(Tracks_cell2,training_cell2,fileidxeach,bleachtime,bleach_var)
          end
   
end

function Trackmate_Analysis_ML_FINAL(Tracks_cell2,training_cell2,iter,bleachtime,bleach_var)
% skip_filter_input = inputdlg('Would you like to skip filtering tests?','Existing Data?',[1 50],{'N'});
skip_filter_input{1} = 'N'
inten_use = 'Y';

if skip_filter_input{1} == 'N' | skip_filter_input{1} == 'n'
%     input_GMM_clustering = inputdlg({'Intensities components','Time Interval','Truncation Point'}, 'Options',...
%         [1 50; 1 50; 1 50], {'4', '1', '3'});
    input_GMM_clustering{1} = '4';
    input_GMM_clustering{2} = '1';
    input_GMM_clustering{3} = '3';
    time_int = 1;
    truncation_pt = 3;
   
    %%
   
    tracks_conc = Tracks_cell2
    num_tracks = length(tracks_conc(:,1));
    disp(num2str(num_tracks));
 
    
    On_time_final = (tracks_conc(:,14))*time_int;
    cat_training = training_cell2
%     inten_use = inputdlg({'Use Mean Intensity?'}, 'GMM Intensity',...
%         [1 50], {'Y'});
    if inten_use == 'Y'
        intensities_final =  cat_training(:,6);
    else
        intensities_final = cat_training(:,7);
    end
    
 
    intensities_final_bound = intensities_final;
     tracks_final_bound = tracks_conc;
     On_time_final_bound = On_time_final;
%     disp('classification tracks')
    disp(length(intensities_final_bound(:,1)))
    



    %% 
    intensities_models_tested = str2num(input_GMM_clustering{1});
    [BestModel_intensities, numComponents_intensities] = GMM_BIC ( intensities_final_bound,intensities_models_tested, true);


    idx_int = cluster(BestModel_intensities, intensities_final_bound);
    cluster_array_int = zeros(length( intensities_final_bound),numComponents_intensities);
    Int_clust = zeros(length( intensities_final_bound),numComponents_intensities);
    Int_values = cell(numComponents_intensities,1);

    for j = 1:numComponents_intensities
        cluster_array_int(:,j) = (idx_int==j);
        Int_clust(:,j) = cluster_array_int(:,j).* intensities_final_bound;
        Int_values{j} = nonzeros(Int_clust(:,j));
    end

    num_of_bins2 = ceil(sqrt(numel( intensities_final_bound))); 
    bin_width = (max( intensities_final_bound)-min( intensities_final_bound))/num_of_bins2;
    mean_intensities = zeros(numComponents_intensities,1);
    figure,
    for i=1:numComponents_intensities
        histogram(Int_values{i},'BinWidth',bin_width,'Normalization','count') %might want to incorporate bin width instead
        hold on
        mean_intensities(i) = mean(Int_values{i});
    end
    xlabel('Intensity (A.U)')
    ylabel('Counts')
    hold off
    single_intensities_ID = min(mean_intensities);
    if numComponents_intensities > 2
        unique_intensities = unique(mean_intensities);
        single_intensities_ID = unique_intensities(1);
    end

    Int_col = find(mean_intensities == single_intensities_ID);

    Single_molecules = Int_clust(:,Int_col);
    find_single_molecules = find(Single_molecules);
    %Quality_Tracks_seg_bound_single = Quality_Tracks_seg_bound_total(find_single_molecules,:);
    On_time_bound_single = On_time_final_bound (find_single_molecules,:);
    intensities_bound_single = intensities_final_bound (find_single_molecules, :);
    disp(mean(intensities_bound_single))
    disp(std(intensities_bound_single))
    tracks_bound_single = tracks_final_bound(find_single_molecules,:);
    if length (intensities_final_bound) < 100
        %Quality_Tracks_seg_bound_single = Quality_Tracks_seg_bound_total;
        On_time_bound_single = On_time_final_bound;
        intensities_bound_single = intensities_final_bound;
        tracks_bound_single =  tracks_final_bound;
    end

    %save( 'Quality_Tracks_seg_bound_single.mat', 'Quality_Tracks_seg_bound_single')
else 
    filenames_tracks = uigetfile('*TrackMate_tracks_bound_single.mat', 'Pick the segmented tracks .mat files','Multiselect', 'on');
    tracks_load = load(filenames_tracks);
    tracks_bound_single = tracks_load.tracks_bound_single;
   
     filename_On_time_final = uigetfile('*TrackMate_On_time_bound_single.mat', 'Pick On time bound single file');
        On_time_load = load(filename_On_time_final);
         On_time_bound_single = On_time_load.On_time_bound_single;
     filename_intensities_final = uigetfile('*TrackMate_intensities_bound_single.mat', 'Pick intensities bound single file');
         intensities_load = load(filename_intensities_final);
         intensities_bound_single = intensities_load.intensities_bound_single;
%      time_input = inputdlg({'Time Interval', 'Truncation Point'}, 'Time Settings', [1 50;1 50], {'1', '3'});
     time_int = 1
     truncation_pt = 3
end

%% 
figure,
histogram(On_time_bound_single,'BinMethod','sqrt','Normalization','pdf');
% input_test = inputdlg('Would you like to test for two exponentials?','Two Exponentials',[1 50],{'N'});
input_test{1} = 'N';
if input_test{1} == 'Y' | input_test{1}== 'y' 
    
        Trackmate_Analysis_Step2 (On_time_bound_single)
    
    
else
    histogram(On_time_bound_single,'BinMethod','sqrt','Normalization','pdf');
%     input_outlier = inputdlg('Eliminate Outliers?', 'Outlier Removals', [1 50], {'Y'});
input_outlier{1} = 'Y';
    if input_outlier{1} =='Y'| input_outlier{1}=='y'
    TF = isoutlier (On_time_bound_single, 'quartiles','ThresholdFactor',4.0);

    On_time_bound_single_filtered = On_time_bound_single(TF~=1,1);
  
    intensities_bound_single_filtered = intensities_bound_single (TF~=1,1);
    tracks_bound_single_filtered = tracks_bound_single(TF~=1,1);
    [est_filtered, ci_filtered, se_filtered] = Fitting_truncExponential (On_time_bound_single_filtered, time_int,truncation_pt);
    
%     input_err = inputdlg({'Do you want to calculate bound time?'}, 'Error Calculator', [1 50], {'Y'});
input_err{1} = 'Y';
    if input_err{1} == 'Y'| input_err {1} == 'y'
%         input_bleach = inputdlg({'Bleach Time', 'Variation in bleach' },'Errors', [1 50; 1 50], {'20', '0.10'});
        [Tbound_filt, Tbound_ci_filt, Tbound_err_filt] = Bound_time_estimator_no_bounds(On_time_bound_single_filtered, est_filtered, bleachtime,bleach_var, truncation_pt);
%         waitfor(msgbox({num2str(Tbound_filt),strcat(num2str(Tbound_ci_filt(1)), ':', num2str(Tbound_ci_filt(2))), num2str(Tbound_err_filt)}, 'Bound Time'));
        
%         input_save_filter = inputdlg('Save Results?', 'Save', [1 50], {'Y'});
    input_save_filter{1} = 'Y';
        if input_save_filter{1} == 'Y' | input_save_filter{1} == 'y'
        Results = struct('BoundTimeFiltered', Tbound_filt, 'BoundTimeCIFiltered', Tbound_ci_filt, 'BoundTimeSTDErrorFiltered', Tbound_err_filt,...
            'TrackDurationFiltered', est_filtered, 'TrackDurationCIFiltered',ci_filtered, 'TrackDurationSEFiltered', se_filtered);
        

%         save_dir_input = inputdlg('Pick a name for the folder','Save Folder', [1 50], {'Analysis Files'});
save_dir_input{1} = ['Analysis Files_' num2str(iter)];
        mkdir(save_dir_input{1}) 
        save(strcat(save_dir_input{1},'/','Results.mat'),'Results');
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single_filtered.mat'),'On_time_bound_single_filtered')
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single_filtered.mat'),'intensities_bound_single_filtered')
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single_filtered.mat'),'tracks_bound_single_filtered')
  
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single.mat'),'On_time_bound_single')
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single.mat'),'intensities_bound_single')
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single.mat'),'tracks_bound_single')
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_TOTAL.mat'),'tracks_conc')
       
        end
    else
%         input_save_filter = inputdlg('Save Results?', 'Save', [1 50], {'Y'});
    input_save_filter = 'Y';
        if input_save_filter{1} == 'Y' | input_save_filter{1} == 'y'
%         save_dir_input = inputdlg('Pick a name for the folder','Save Folder', [1 50], {'Analysis Files'});
save_dir_input = 'Analysis Files';
        mkdir(save_dir_input{1}) 
        Results = struct('TrackDurationFiltered', est_filtered, 'TrackDurationCIFiltered',ci_filtered, 'TrackDurationSEFiltered', se_filtered);
        save(strcat(save_dir_input{1},'/','Results_%d.mat'),'Results',iter);
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single_filtered_%d.mat'),'On_time_bound_single_filtered',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single_filtered_%d.mat'),'intensities_bound_single_filtered',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single_filtered_%d.mat'),'tracks_bound_single_filtered',iter)
  
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single_%d.mat'),'On_time_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single_%d.mat'),'intensities_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single_%d.mat'),'tracks_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_TOTAL_%d.mat'),'tracks_conc',iter)
        %save(strcat(save_dir_input{1},'/','PMtracker_PSFS_TOTAL.mat'),'psfs_final')
        end
    end

    
    else
        [est_original, ci_original, se_original] = Fitting_truncExponential (On_time_bound_single,time_int,truncation_pt);
%         input_err = inputdlg({'Do you want to calculate bound time?'}, 'Error Calculator', [1 50], {'Y'});
input_err = 'Y';
    if input_err{1} == 'Y'| input_err {1} == 'y'
%         input_bleach = inputdlg({'Bleach Time', 'Variation in bleach' },'Errors', [1 50; 1 50], {'20','0.10'});
        [Tbound_original, Tbound_ci_original, Tbound_err_original] = Bound_time_estimator_no_bounds(On_time_bound_single, est_original, 20,0.1, truncation_pt);

% waitfor(msgbox({num2str(Tbound_original),strcat(num2str(Tbound_ci_original(1)), ':', num2str(Tbound_ci_original(2))), num2str(Tbound_err_original)}, 'Bound Time'));
%         input_save = inputdlg('Save Results?', 'Save', [1 50], {'Y'});
    input_save = 'Y';
        if input_save{1} == 'Y' | input_save{1} == 'y'
%         save_dir_input = inputdlg('Pick a name for the folder','Save Folder', [1 50], {'Analysis Files'});
save_dir_input = 'Analysis Files';
        mkdir(save_dir_input{1}) 
        Results = struct('BoundTime', Tbound_original, 'BoundTimeCI', Tbound_ci_original, 'BoundTimeSTDError',  Tbound_err_original,...
            'TrackDuration', est_original, 'TrackDurationCI',ci_original, 'TrackDurationSE', se_original);
         Results = struct('TrackDuration', est_original, 'TrackDurationCI',ci_original, 'TrackDurationSE', se_original);
        save(strcat(save_dir_input{1},'/','Results_%d.mat'),'Results',iter);
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_final_%d.mat'),'On_time_final_bound',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_final_%d.mat'),'intensities_final_bound',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single_%d.mat'),'On_time_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single_%d.mat'),'intensities_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single_%d.mat'),'tracks_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_TOTAL_%d.mat'),'tracks_conc',iter)
        
        end
    else
%         input_save = inputdlg('Save Results?', 'Save', [1 50], {'Y'});
    input_save = 'Y';
        if input_save{1} == 'Y' | input_save{1} == 'y'
%         save_dir_input = inputdlg('Pick a name for the folder','Save Folder', [1 50], {'Analysis Files'});
save_dir_input = 'Analysis Files';
        mkdir(save_dir_input{1}) 
        Results = struct('TrackDuration', est_original, 'TrackDurationCI',ci_original, 'TrackDurationSE', se_original);
        save(strcat(save_dir_input{1},'/','Results_%d.mat'),'Results',iter);
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_final_%d.mat'),'On_time_final_bound',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_final_%d.mat'),'intensities_final_bound',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_On_time_bound_single_%d.mat'),'On_time_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_intensities_bound_single_%d.mat'),'intensities_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_bound_single_%d.mat'),'tracks_bound_single',iter)
        save(strcat(save_dir_input{1},'/','Trackmate_tracks_TOTAL_%d.mat'),'tracks_conc',iter)
        
        end
    end
    end   
end
end





