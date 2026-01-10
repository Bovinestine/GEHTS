clear all;
clc;
cd('Z:\03_Team\In_situ_sequencing\[5]program_run\221017_c4')

%% Read particle map %%
dCsv=dir('*.csv');
T=readtable(dCsv.name);
T=T{:,:};
A=T;
B=fliplr(T);
P=[];
for i=1:size(T)[1];
        if rem(i,2) ==0
            P=[P ; B(i,:)];
        else
            P=[P ; A(i,:)];
        end
end



%% Read raw images %%
directory = cd;
cd(strcat(directory,'\raw_image'));
dList=dir('*.tif');


n=0;
base_max= 4;
ch_max= 6;
tile_max=100;

%% Create channel no.5
for base=1:base_max
    for xy_no= 1:tile_max   
    for ch =5 :53
        base=num2str(base);
        xy=sprintf('%03d',xy_no);
        ch=num2str(ch);
    im_black = zeros(512,512);
    imwrite(uint16(im_black),strcat('base',base,'xy',xy,'c',ch,'.tif')) % name match
    end 
    end
end

%% Create Cy7 channel

dList=dir('*.tif');
image_max=zeros(512,'uint16');
for n= 1:length(dList)
    nums=regexp(dList(n).name,'\d*','match');
    base_num=cell2mat(nums(1,1));
    xypo=cell2mat(nums(1,2));
    chpo=cell2mat(nums(1,3));
    image_data=imread(dList(n).name);
    leng=size(image_data);

        for i=1:leng(1)
            for j=1:leng(2)
                image_max(i,j)=max(image_max(i,j),image_data(i,j));
            end
        end

    if chpo=='5'
        imwrite(uint16(image_max),strcat('base',base_num,'xy',xypo,'c6.tif'))
        image_max=zeros(512,'uint16');
    elseif chpo=='6'
        image_max=zeros(512,'uint16');

    end
end
dList=dir('*.tif');


%% Create folder to be aligned %% 
newF=strcat(directory,'\in_situ'); 
if not(isfolder(newF))
    mkdir(newF);
end
cd(newF);
%% Distribute images %% 

for i=1:size(P)[1];
    for j=1:size(P)[2];
        particle_name=strcat('particle',num2str(P(i,j)));
        if not(isfolder(particle_name))
            mkdir(particle_name);
            mkdir(strcat(particle_name,'\','image'));
            
            for m=1:base_max;
            for k=1:ch_max;
                f_name = strcat(particle_name,'\image\','Tiled_base',num2str(m),'_c', num2str(k));
                mkdir(f_name);
            end
            
            end
        end
        n=n+1; 
        if n>tile_max;
            n=1;
        end 
        xy_po=sprintf('%03d',n);
        
        for m=1:base_max
            for l=1:numel(dList)/base_max;
                nums=regexp(dList(numel(dList)/base_max*(m-1)+l).name,'\d*','match');
                base_num=cell2mat(nums(1,1));
                xypo=cell2mat(nums(1,2));
                ch_num=cell2mat(nums(1,3));
                if xypo==xy_po
                    d=dir([strcat(particle_name,'\image\','Tiled_base',base_num,'_c',ch_num) '/*tif']);
                    tile_num=numel(d)+1;
                    f_name=strcat(dList(numel(dList)/base_max*(m-1)+l).folder,'\',dList(numel(dList)/base_max*(m-1)+l).name);
                    dir_name=strcat(particle_name,'\image\','Tiled_base',base_num,'_c',ch_num,'\','tile',num2str(tile_num),'.tif');
                    copyfile(f_name,dir_name);
                end
            end
        end
        
    end
end
strcat("Image distributed")

%% fill in the tile images to squared number %% 
dir_num=dir;
imsize=512;

for l= 3:numel(dir_num);
    if dir_num(l).name(1:5) == 'parti'
    for c=1:ch_max;
        for b=1:base_max;
            dir_name=strcat(dir_num(l).name,'\image\','Tiled_base',num2str(b),'_c',num2str(c));
            d=dir([dir_name '/*tif']);
            for t= length(d):tile_max
                if sqrt(t)==floor(sqrt(t))
                    break
                end 
            end 
            for tile_num=numel(d)+1:t
                f_name=strcat('tile',num2str(tile_num),'.tif');
                im_black = zeros(imsize,imsize);
                imwrite(uint16(im_black),strcat(dir_name,'\',f_name))
            end 
        end
    end
    end
end
strcat("tiling finished")

%% create Tiled.csv files %% 

for l=3:numel(dir_num);
    if dir_num(l).name(1:5) == 'parti'
    particle_dir=strcat(newF,'\',dir_num(l).name);
    cd(particle_dir)
    M={'Metadata_position' 'Tile_xPos' 'Tile_yPos' 'Hyb_step' 'Image_PathName_G' 'Image_FileName_G' 'Image_PathName_C' 'Image_FileName_C' 'Image_PathName_A' 'Image_FileName_A' 'Image_PathName_Nuclei' 'Image_FileName_Nuclei' 'Image_PathName_T' 'Image_FileName_T' 'Image_PathName_Spec_blob' 'Image_FileName_Spec_blob' 'Image_PathName_General_blob' 'Image_FileName_General_blob'};
    dir_name=strcat(particle_dir,'\image\Tiled_base1_c1');
    d=dir([dir_name '/*tif']);
    meta=1;
    xpo=-512;
    ypo=0;
    y=0;
    for j=1:numel(d);
        h=1;
        xpo=xpo+imsize;
        if xpo>=sqrt(numel(d))*imsize
            xpo=0;
            y=y+1;
        end 
        ypo=y*imsize;
        for i=1:base_max;
          tile=strcat('tile',num2str(meta),'.tif');
          M=[M;{meta,xpo,ypo,h,strcat(particle_dir,'\image\Tiled_base',num2str(h),'_c1\'),tile,strcat(particle_dir,'\image\Tiled_base',num2str(h),'_c2\'),tile,strcat(particle_dir,'\image\Tiled_base',num2str(h),'_c3\'),tile,strcat(particle_dir,'\image\Tiled_base1_c4\'),tile,strcat(particle_dir,'\image\Tiled_base',num2str(h),'_c5\'),tile,strcat(particle_dir,'\image\Tiled_base',num2str(h),'_c6\'),tile ,strcat(particle_dir,'\image\Tiled_base1_c6\'),tile}];
          h=h+1;
        end 

        meta=meta+1;
        
    end
    M2=cell2table(M(2:end,:),'VariableNames',M(1,:));
    writetable(M2,'Tiled.csv')
    cd(newF);
    end
end
strcat("csv file created")

%% Create merged Tiled.csv files %%

X=[];
for l=3:numel(dir_num);
    if dir_num(l).name(1:5) == 'parti'
    particle_dir=strcat(newF,'\',dir_num(l).name);
    cd(particle_dir)
    T=readtable('Tiled.csv','Delimiter',',');
    X=[X;T];
    meta=1;
    for x=1:height(X)/base_max;
        for i=1:base_max;
            X.Metadata_position(base_max*(x-1)+i)=meta;
        end
        meta=meta+1;
    end
    end
end
cd(newF);
writetable(X,'Tiled.csv')
strcat("merged_csv file created")

