clear all;
clc;
%% read all files
directory='Z:\03_Team\In_situ_sequencing\[5]program_run\220310c_4';
cd(strcat(directory,'\raw_image'));

dList=dir('*.tif');

well_num='027';
dest= strcat(directory,'\',well_num);
if(not(isdir(dest)))
    mkdir(dest);
    mkdir(strcat(dest,'\image'))
end

xy_max=27;
base_max=4;

%% move selected files
for n= 1:length(dList)
    nums=regexp(dList(n).name,'\d*','match');
    base_num=cell2mat(nums(1,1));
    xypo=cell2mat(nums(1,2));
    chpo=cell2mat(nums(1,3));
    image_data=imread(dList(n).name);
    leng=size(image_data);
    if xypo==well_num
        copyfile(strcat(dList(n).folder,'\',dList(n).name),strcat(dest,'\image\base',base_num,'xy',xypo,'c',chpo,'.tif'));
    end
    
end 

strcat('move_file_finished')

%% read all files
cd(strcat(dest,'\image'));


%% Create channel no.5
for base=1:base_max
    for xy_no= xy_max:xy_max   
    for ch =5 :5
        base=num2str(base)
        xy=sprintf('%03d',xy_no)
        ch=num2str(ch)
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
        imwrite(uint16(10*image_max),strcat('base',base_num,'finalimage.tif'))
        image_max=zeros(512,'uint16');
    end
end

strcat('move_file_finished')
