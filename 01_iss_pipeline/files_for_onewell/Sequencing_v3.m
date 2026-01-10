% in situ sequencing pipeline
% Ideally keep one file for one experiment
% Five major functions: 1)Tiling, 2)Decoding, 3)Thresholding, 4)Plotting 
% All variables ending with _YN are yes or no (1 or 0) questions.
%
% Sequencing v3
% Xiaoyan, 2017


clear, clc; close all; drawnow;

x_num='027'
directory=strcat('Z:\03_Team\In_situ_sequencing\[5]program_run\220310c_4\in_situ\particle0\Results\tile27');
%directory=strcat('Z:\03_Team\In_situ_sequencing\[5]program_run\220310c_4\',x_num);
cd(directory);


% Choose functions to run
run_Tiling_YN = 0; 
run_Decode_YN = 01;
run_Threshold_YN =01;
run_Plotting_Global_YN =01;
%================================================
% set parameters
%----------------------------
% Tiling_Sequencing
    t.folder_image = strcat(directory,'\image'); % preferably full path name
    t.filename_base_prefix = 'base';  % keep single quote marks
        t.in_subfolder_YN = 0;
    t.filename_channel_prefix = strcat('xy',x_num,'c');
    t.filename_suffix = '.tif';
    t.base_start = 1;     t.base_end = 4;       
    t.channel_start = 1;  t.channel_end = 6;
    t.tile_size = 510;
    t.channel_order = { 'G','C','A','Nuclei','T','General_stain'};
    t.CSV_filename_prefix = 'Tiled';
%----------------------------
% Decoding_Sequencing
    d.input_file = 'blobs.csv';
    % don't change this unless your input file has some weird form
    d.General_Alignment_ACGT_column_number = [4,5,6,7,8,9];    % use 0 if any of them is MISSING in the file
    d.XYPosition_ParentCell_column_number = [10,11,12];
    
    d.num_hybs = 4;
    d.taglist = 'ID_NatMeth.csv';   % old .m taglist or .csv file with columns: code, name, symbol(optional), no header
    d.csv_file_contain_tile_position = 'Tiled.csv';
    d.output_directory = 'Decoding';   
    % options
    d.check_parent_cell_YN = 0;       
    d.check_alignment_YN = 0;
        alignment_min_threshold = 1.8;
    d.abnormal_sequencing_YN = 0;
        d.sequencing_order = '12340';  % keep the quote marks, same length as

%----------------------------3
% Threshold_Sequencing
    q.quality_threshold = 0.3;        
    q.general_stain_threshold = 0;
%----------------------------
% Plotting_global_Sequencing
    p.background_image = strcat(directory, '\image\base1finalimage.tif'); 
    p.scale = 1;		% image scale
    p.I_want_to_plot_on_white_backgound = 1; %background image가 없으면 1로 해라. 있으면 0으로 하고 위에 경로설정 제대로 하면 됨.
    % options
    p.exclude_NNNN_YN = 0;
    p.plot_reads_beforeQT_YN = 0;
    p.plot_ref_general_stain = 0; 
%================================================
if run_Tiling_YN || run_Decode_YN || run_Threshold_YN || run_Plotting_Global_YN
else
    error('Choose at least one function.');
end

if run_Tiling_YN
    seqtiling(t);
end
if run_Decode_YN
    decoding(d);
end
if run_Threshold_YN
    qthreshold(d.output_directory, q);
end 
if run_Plotting_Global_YN
    seqplotting(d.output_directory, d.taglist, q, p);
end

%% Save figures

figHandles = findall(0,'Type','figure'); 
if not(isfolder(strcat(directory,'\figure')))
    mkdir(strcat(directory,'\figure'))
end 
for i = 1 : numel(figHandles)
export_fig(strcat(directory,'\figure\',num2str(i)),'-png',figHandles(i))
end
clear;   
    