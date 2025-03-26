 Setup basic parameters

addpath(genpath('/home/users/gomezmi/conn'));
addpath(genpath('/home/users/gomezmi/spm12'));

NSUBJECTS = 167;
cwd = pwd; 
FUNCTIONAL_FILE = cellstr(conn_dir('sub*Deformed-To-Template40.nii'));  
STRUCTURAL_FILE = cellstr(conn_dir('sub-*_T2W-To-Template40.nii'));
MOTION_FILE = cellstr(conn_dir('sub*motion_1600_6.tsv')); 
TR = 0.392; 



% Load left and Right hemisphere masks (assumed to be binary masks) 
left_hemisphere_mask = spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/left_NoMasks.nii'));
Right_hemisphere_mask = spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/right_NoMasks.nii'));
LB_mask = spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Broca_left_week40_GM.nii'));
RB_mask =spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Broca_right_week40_GM.nii'));
LW_mask =spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Wernicke_left_week40_GM.nii'));
RW_mask =spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Wernicke_right_week40_GM.nii'));
LA_mask =spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Auditory_left_week40_GM.nii'));
RA_mask =spm_read_vols(spm_vol('/home/shared/pi-bokdea/dHCP/language/week40/Auditory_right_week40_GM.nii'));

% Flatten the masks to get linear indices of the voxels in each hemisphere 
left_voxel_indices = find(left_hemisphere_mask); 
Right_voxel_indices = find(Right_hemisphere_mask);
LB_voxel_indices = find(LB_mask);
RB_voxel_indices = find(RB_mask);
LW_voxel_indices = find(LW_mask);
RW_voxel_indices = find(RW_mask);
LA_voxel_indices = find(LA_mask);
RA_voxel_indices = find(RA_mask);

% Get the number of left and Right voxels 
n_left_voxels = length(left_voxel_indices);
n_Right_voxels = length(Right_voxel_indices);
n_LB_voxels = length(LB_voxel_indices);
n_RB_voxels = length(RB_voxel_indices);
n_LW_voxels = length(LW_voxel_indices);
n_RW_voxels = length(RW_voxel_indices);
n_LA_voxels = length(LA_voxel_indices);
n_RA_voxels = length(RA_voxel_indices);

subj_data = readlines('/home/shared/pi-bokdea/dHCP/language/week40/fconn/gracemia_denoised.txt');
num_subj=length(subj_data)
%%% The output_list file contains the name of the output (use ID number code) i.e. 00065757
output_table = readlines('/home/shared/pi-bokdea/dHCP/language/week40/fconn/gracemia_output.txt');

for x=1:num_subj;
subj= subj_data{x} ;
output= output_table{x} ;


%%%%start for loop here    
    % Load the NIfTI output file from the 3dTproject command, first 
 nifti_file_path=fullfile(cwd, [num2str(subj)]);

    subject_bold_data = spm_read_vols(spm_vol(nifti_file_path));  
    

    % Reshape the subject data to make it easier to work with 
    [nx, ny, nz, nt] = size(subject_bold_data);
    subject_bold_data = reshape(subject_bold_data, [], nt);  % Flatten the spatial dimensions

    % Extract time series for all left and Right hemisphere voxels 
    left_voxel_ts = subject_bold_data(left_voxel_indices, :);  % Size: [n_left_voxels, nt]
    Right_voxel_ts = subject_bold_data(Right_voxel_indices, :);  % Size: [n_Right_voxels, nt]
    LB_voxel_ts = subject_bold_data(LB_voxel_indices, :);
    RB_voxel_ts = subject_bold_data(RB_voxel_indices, :);
    LW_voxel_ts = subject_bold_data(LW_voxel_indices, :);
    RW_voxel_ts = subject_bold_data(RW_voxel_indices, :);
    LA_voxel_ts = subject_bold_data(LA_voxel_indices, :);
    RA_voxel_ts = subject_bold_data(RA_voxel_indices, :);

    % Normalize the time series to have zero mean and unit variance (z-scoring)
    left_voxel_ts = zscore(left_voxel_ts, 0, 2);  % Normalize along the time dimension
    Right_voxel_ts = zscore(Right_voxel_ts, 0, 2);
    LB_voxel_ts = zscore(LB_voxel_ts, 0, 2);
    RB_voxel_ts = zscore(RB_voxel_ts, 0, 2);
    LW_voxel_ts = zscore(LW_voxel_ts, 0, 2);
    RW_voxel_ts = zscore(RW_voxel_ts, 0, 2);
    LA_voxel_ts = zscore(LA_voxel_ts, 0, 2);
    RA_voxel_ts = zscore(RA_voxel_ts, 0, 2);
    % Initialize the output correlation map
    avg_correlation_map = zeros(size(left_voxel_indices));
    avg_LB_correlation_map = zeros(size(LB_voxel_indices));
    avg_RB_correlation_map = zeros(size(RB_voxel_indices));
    avg_LW_correlation_map = zeros(size(LW_voxel_indices));
    avg_RW_correlation_map = zeros(size(RW_voxel_indices));
    avg_LA_correlation_map = zeros(size(LA_voxel_indices));
    avg_RA_correlation_map = zeros(size(RA_voxel_indices));

    % Process in chunks LBtoRB
        chunk_size = 1000; 
    for k = 1:chunk_size:n_LB_voxels
        chunk_indices = k:min(k+chunk_size-1, n_LB_voxels);

        % Calculate correlation matrix for the current chunk
        correlation_matrix_chunk = LB_voxel_ts(chunk_indices, :) * RB_voxel_ts' / (nt - 1);

        % Calculate the mean correlation for each voxel in the current chunk
        avg_LB_correlation_map(chunk_indices) = mean(correlation_matrix_chunk, 2);
    end

    % Initialize the output volume with zeros
    full_correlation_map = zeros(nx*ny*nz, 1);
    full_correlation_map(LB_voxel_indices) = avg_LB_correlation_map;

    % Reshape the correlation map back to 3D
    full_correlation_map = reshape(full_correlation_map, [nx, ny, nz]);

    % Save the average correlation map for this subject
output_vol = spm_vol(nifti_file_path);  % Copy the header information from the input file

% Check if output_vol is a structure array and select the first element
if numel(output_vol) > 1
    output_vol = output_vol(1);
end

% change filename below to what you would like - suggest change ~ID to ID number
output_vol.fname = fullfile(cwd, ['avg_correlation_subject_' num2str(output) '_LB_to_RB.nii']);
spm_write_vol(output_vol, full_correlation_map);

   % Process in chunks LBtoLW
        chunk_size = 1000; 
    for k = 1:chunk_size:n_LB_voxels
        chunk_indices = k:min(k+chunk_size-1, n_LB_voxels);

        % Calculate correlation matrix for the current chunk
        correlation_matrix_chunk = LB_voxel_ts(chunk_indices, :) * LW_voxel_ts' / (nt - 1);

        % Calculate the mean correlation for each voxel in the current chunk
        avg_LB_correlation_map(chunk_indices) = mean(correlation_matrix_chunk, 2);
    end

    % Initialize the output volume with zeros
    full_correlation_map = zeros(nx*ny*nz, 1);
    full_correlation_map(LB_voxel_indices) = avg_LB_correlation_map;

    % Reshape the correlation map back to 3D
    full_correlation_map = reshape(full_correlation_map, [nx, ny, nz]);

    % Save the average correlation map for this subject
output_vol = spm_vol(nifti_file_path);  % Copy the header information from the input file

% Check if output_vol is a structure array and select the first element
if numel(output_vol) > 1
    output_vol = output_vol(1);
end

% change filename below to what you would like - suggest change ~ID to ID number
output_vol.fname = fullfile(cwd, ['avg_correlation_subject_' num2str(output) '_LB_to_LW.nii']);
spm_write_vol(output_vol, full_correlation_map);

   % Process in chunks LBtoRW
        chunk_size = 1000; 
    for k = 1:chunk_size:n_LB_voxels
        chunk_indices = k:min(k+chunk_size-1, n_LB_voxels);

        % Calculate correlation matrix for the current chunk
        correlation_matrix_chunk = LB_voxel_ts(chunk_indices, :) * RW_voxel_ts' / (nt - 1);

        % Calculate the mean correlation for each voxel in the current chunk
        avg_LB_correlation_map(chunk_indices) = mean(correlation_matrix_chunk, 2);
    end

    % Initialize the output volume with zeros
    full_correlation_map = zeros(nx*ny*nz, 1);
    full_correlation_map(LB_voxel_indices) = avg_LB_correlation_map;

    % Reshape the correlation map back to 3D
    full_correlation_map = reshape(full_correlation_map, [nx, ny, nz]);

    % Save the average correlation map for this subject
output_vol = spm_vol(nifti_file_path);  % Copy the header information from the input file

% Check if output_vol is a structure array and select the first element
if numel(output_vol) > 1
    output_vol = output_vol(1);
end

% change filename below to what you would like - suggest change ~ID to ID number
output_vol.fname = fullfile(cwd, ['avg_correlation_subject_' num2str(output) '_LB_to_RW.nii']);
spm_write_vol(output_vol, full_correlation_map);

   % Process in chunks LBtoLA
        chunk_size = 1000; 
    for k = 1:chunk_size:n_LB_voxels
        chunk_indices = k:min(k+chunk_size-1, n_LB_voxels);

        % Calculate correlation matrix for the current chunk
        correlation_matrix_chunk = LB_voxel_ts(chunk_indices, :) * LA_voxel_ts' / (nt - 1);

        % Calculate the mean correlation for each voxel in the current chunk
        avg_LB_correlation_map(chunk_indices) = mean(correlation_matrix_chunk, 2);
    end

    % Initialize the output volume with zeros
    full_correlation_map = zeros(nx*ny*nz, 1);
    full_correlation_map(LB_voxel_indices) = avg_LB_correlation_map;

    % Reshape the correlation map back to 3D
    full_correlation_map = reshape(full_correlation_map, [nx, ny, nz]);

    % Save the average correlation map for this subject
output_vol = spm_vol(nifti_file_path);  % Copy the header information from the input file

% Check if output_vol is a structure array and select the first element
if numel(output_vol) > 1
    output_vol = output_vol(1);
end

% change filename below to what you would like - suggest change ~ID to ID number
output_vol.fname = fullfile(cwd, ['avg_correlation_subject_' num2str(output) '_LB_to_LA.nii']);
spm_write_vol(output_vol, full_correlation_map);
  % Process in chunks LBtoRA
        chunk_size = 1000; 
    for k = 1:chunk_size:n_LB_voxels
        chunk_indices = k:min(k+chunk_size-1, n_LB_voxels);

        % Calculate correlation matrix for the current chunk
        correlation_matrix_chunk = LB_voxel_ts(chunk_indices, :) * RA_voxel_ts' / (nt - 1);

        % Calculate the mean correlation for each voxel in the current chunk
        avg_LB_correlation_map(chunk_indices) = mean(correlation_matrix_chunk, 2);
    end

    % Initialize the output volume with zeros
    full_correlation_map = zeros(nx*ny*nz, 1);
    full_correlation_map(LB_voxel_indices) = avg_LB_correlation_map;

    % Reshape the correlation map back to 3D
    full_correlation_map = reshape(full_correlation_map, [nx, ny, nz]);

    % Save the average correlation map for this subject
output_vol = spm_vol(nifti_file_path);  % Copy the header information from the input file

% Check if output_vol is a structure array and select the first element
if numel(output_vol) > 1
    output_vol = output_vol(1);
end

% change filename below to what you would like - suggest change ~ID to ID number
output_vol.fname = fullfile(cwd, ['avg_correlation_subject_' num2str(output) '_LB_to_RA.nii']);
spm_write_vol(output_vol, full_correlation_map);
end 


