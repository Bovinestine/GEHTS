clear all; 
clc;

%% Tiled File

directory = './'

%cd(strcat(directory,'\in_situ1'));
%lob1=readtable('blobs.csv')

%cd(strcat(directory,'\in_situ2'))
%Blob2=readtable('blobs.csv');
%for i=1:height(Blob2)
%Blob2.ImageNumber(i)=Blob2.ImageNumber(i)+600;
%end 

%Blob= vertcat(Blob1,Blob2);
%cd(strcat(directory,'\in_situ'))
%writetable(Blob, 'blobs.csv')

%% copy needed files

cd(strcat(directory,'\in_situ'))
dir_num=dir;
folder_names='parti'
base_max=4;



%% Result files distribution _ 약물처리일때
im_nums=0;

for l=3:numel(dir_num)
    if dir_num(l).name(1:5)==folder_names;
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        cd(particle_dir);
        T=readtable('Tiled.csv','Delimiter',',');
        im_nums=im_nums+height(T)
        cd(dir_num(l).folder);
        X= readtable('blobs.csv');
        MAT=table();
        for m=1:length(X.ImageNumber)
            if (X.ImageNumber(m)<=im_nums) && (X.ImageNumber(m)>(im_nums-height(T)))
                X.ImageNumber(m)=X.ImageNumber(m)-(im_nums-height(T));
                X.Metadata_position(m)=X.Metadata_position(m)-(im_nums-height(T))/base_max;
                MAT=[MAT;X(m,:)];
            end
        end
        cd(particle_dir);
        writetable(MAT,'blobs.csv')
    end
end
strcat("result file distributed")
%copyfile(strcat(directory,'\in_situ\blobs.csv'),strcat(directory,'\in_situ\particle0\blobs.csv'))



%% Tile distribution
cd(strcat(directory,'\in_situ'))

for l=  3:numel(dir_num)
    if dir_num(l).name(1:5)==folder_names;
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        cd(particle_dir);
        if not(isfolder('Results'))
             mkdir('Results');
        end
        if not(isfile('blobs.csv'))
            continue
        end
        B=readtable('blobs.csv','Delimiter',','); %% _ <- delete
        T=readtable('Tiled.csv','Delimiter',',');
        pos_max=T.Metadata_position(length(T.Metadata_position));
        BLOB=table();
        for j=1:pos_max
            f_name=strcat(particle_dir,'\Results\tile',num2str(j));
                if not(isfolder(f_name))
                    mkdir(f_name);
                end 
            for i=1:length(B.ImageNumber)
            if(B.Metadata_position(i)==j)
                B.ImageNumber(i)=B.ImageNumber(i)-(B.Metadata_position(i)-1)*base_max;
                B.Metadata_position(i)=1;
                BLOB=[BLOB;B(i,:)];
            end
            end
            writetable(BLOB,strcat(f_name,'\blobs.csv'))
            BLOB=table();
            TILE=table();
            for m=1:length(T.Metadata_position)
                if T.Metadata_position(m) == j
                     T.Metadata_position(m)=1;
                     T.Tile_xPos(m)=0;
                     T.Tile_yPos(m)=0;
                    TILE=[TILE;T(m,:)];
                end
            end
            writetable(TILE,strcat(f_name,'\Tiled.csv'))
        end
    end 
end
strcat("Tiles distributed")

%% copy needed files
for l=3:numel(dir_num)
    if dir_num(l).name(1:5)==folder_names
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        cd(particle_dir)
        if not(isfile('ID_NatMeth_4.csv'))
            copyfile(strcat(directory,'\files\ID_NatMeth_4.csv'),'ID_NatMeth_4.csv')
        end
        cd(strcat(particle_dir,'\Results'))
        tile_dir= dir;

        for i=3:numel(tile_dir)
            cd(strcat(particle_dir,'\Results\',tile_dir(i).name))
        if not(isfile('histread.m'))
            copyfile(strcat(directory,'\files\histread.m'),'histread.m')
        end 
        if not(isfile('ID_NatMeth_1.csv'))
            copyfile(strcat(directory,'\files\ID_NatMeth_4.csv'),'ID_NatMeth_4.csv')
        end
        end
                                                                                        
    end 
end 
strcat("Files copied")

    

%% background image  %no run


%% Run Sequencing.m files
 set(groot, 'defaultFigureVisible', 'off');
for l=3:numel(dir_num)
    if dir_num(l).name(1:5)==folder_names
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        cd(strcat(particle_dir,'\Results'))
        tile_dir= dir;
        for i=3:numel(tile_dir)
        if tile_dir(i).name(1:4)=='tile'
        cd(strcat(particle_dir,'\Results\',tile_dir(i).name))
    % Decoding_Sequencing
     A=readtable('blobs.csv');
     if isempty(A)  
        continue
     end
     d.input_file = 'blobs.csv';

     % don't change this unless your input file has some weird form
     d.General_Alignment_ACGT_column_number = [4,5,6,7,8,9];    % use 0 if any of them is MISSING in the file
     d.XYPosition_ParentCell_column_number = [10,11,12];
    
     d.num_hybs = base_max;
     d.taglist = 'ID_NatMeth_4.csv';   % old .m taglist or .csv file with columns: code, name, symbol(optional), no header
     d.csv_file_contain_tile_position = 'Tiled.csv';
    d.output_directory = 'Decoding';   
    % options
    d.check_parent_cell_YN = 0;       
    d.check_alignment_YN = 0;
        alignment_min_threshold = 1.8;
    d.abnormal_sequencing_YN = 0;
        d.sequencing_order = '12340';  % keep the quote marks, same length as
    decoding(d);

    strcat("decoding fin")
    %----------------------------
    % Threshold_Sequencing
        q.quality_threshold = 0.3;        
        q.general_stain_threshold = 0;
        qthreshold(d.output_directory, q);   
        strcat("threshold set")
            end
        end 
    end
end 



%% Gene Count

cd(directory)
dCsv=dir('*.csv');
T=readtable(dCsv.name);
T=T{:,:};
A=T;
B=fliplr(T);
P=[];
for i=1:size(T)[1];
        if rem(i,2) ==1
            P=[P ; B(i,:)];
        else
            P=[P ; A(i,:)];
        end
end

cd(strcat(directory,'\in_situ'))
for l=3 :numel(dir_num)
    if dir_num(l).name(1:5)==folder_names
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        [~,txt,~]=xlsread(strcat(particle_dir,'\ID_NatMeth_4.csv'));
        cd(strcat(particle_dir,'\Results'))
        tile_dir= dir;
        Gene_counts=table();
        drug={'*drug no',};
        particle_no=regexp(dir_num(l).name,'\d*','match');
        p_num=str2num(cell2mat(particle_no));
        [row,col]= find(P==p_num);
        p_pos=table();
        p_pos.row=row;
        p_pos.col=col;
        p_pos=sortrows(p_pos,'row');
        gene_names=unique(string(txt(:,2)));
        gene_names=[gene_names;'NNNN';'Homomer']
        Gene_counts.GeneName=gene_names;
            
    for i=3:numel(tile_dir)
        if not(isfolder(strcat(particle_dir,'\Results\',tile_dir(i).name)))
            COUNTS=table();
            COUNTS.Count=zeros(height(COUNTS),1)
      
        else
        
        cd(strcat(particle_dir,'\Results\',tile_dir(i).name))
        A=[];
        if isfile('blobs.csv')
        A=readtable('blobs.csv');
        end
            if isempty(A)
                if tile_dir(i).name(1:4)=='tile'
                strcat(tile_dir(i).name);
                COUNTS=table();
                COUNTS.GeneName=gene_names;
                COUNTS.Count=zeros(height(COUNTS),1);
                end
            else 
                if tile_dir(i).name(1:4)=='tile'
                    strcat(tile_dir(i).name)
                cd(strcat(particle_dir,'\Results\',tile_dir(i).name,'\Decoding'))
                COUNTS=readtable('QT_0.3_0_gene_n_count.csv','Delimiter',',');
                    for b=1:length(COUNTS.GeneName)
                        if contains(COUNTS.GeneName(b),gene_names)
                        else
                            COUNTS(b,:)=[];
                        end 
                    end
                end 

                for a=1:length(gene_names)
                    if contains(cellstr(gene_names(a)),COUNTS.GeneName)
                    else
                        COUNTS.GeneName(length(COUNTS.GeneName)+1)=cellstr(gene_names(a));
                        gene_names(a);
                    end
                end 
            end
        end
        
        tile_num=regexp(tile_dir(i).name,'\d*','match');
        t_num=str2num(cell2mat(tile_num));
        try 
            w_num=10*(p_pos.row(t_num)-1)+p_pos.col(t_num); %%20;
        
            xy_po= sprintf('%03d',w_num);
            COUNTS=sortrows(COUNTS,'GeneName');
            Gene_counts.GeneName=COUNTS.GeneName;
            xypo=strcat('no',xy_po);
            Gene_counts.(xypo) = COUNTS.Count;
            drug={drug{:,:},p_num};
        catch
 
       end
        
        if i==numel(tile_dir)
            strcat("a")
        Gene_counts=[Gene_counts;drug]
        writetable(Gene_counts,strcat(particle_dir,'\gene_counts.csv'))
        end
    end
    end
end


%% Merge files 
G1=table();

G1.GeneName=COUNTS.GeneName;
for l=11:11
    if dir_num(l).name(1:5)==folder_names
        particle_dir=strcat(dir_num(l).folder,'\',dir_num(l).name);
        cd(particle_dir)
        G=readtable('gene_counts.csv');
        G1=outerjoin(G1,G,'MergeKeys',true)
    end
end 

G2=table2cell(G1);

G3=G2';
g_names=G1.Properties.VariableNames
g_names1=g_names';
G4=[g_names1, G3];


G5=sortrows(G4,1)
cd(dir_num(l).folder)
writetable( cell2table(G5), 'well_gene_counts38.csv', 'writevariablenames', false, 'quotestrings', true)


