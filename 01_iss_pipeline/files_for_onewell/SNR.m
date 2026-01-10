clear all; 
clc;

%% SNR

directory = 'Z:\03_Team\In_situ_sequencing\[5]program_run\220713t_1\in_situ'; % Directory where blobs.csv file exists
cd(directory);

X=readtable('blobs.csv');
base_num=1;

for i=1:base_num
    intensity=max(max(X.Intensity_MaxIntensity_aligned_A,X.Intensity_MaxIntensity_aligned_C),X.Intensity_MaxIntensity_aligned_G);
    norm_dist=string(strcat('=NORM.DIST(A',num2str(2),',AVERAGE(A:A),STDEV.P(A:A),FALSE)'));
    average=mean(intensity(X.ObjectNumber(length(X.ImageNumber))*(i-1)+1:X.ObjectNumber(length(X.ImageNumber))*i));
    int=intensity(X.ObjectNumber(length(X.ImageNumber))*(i-1)+1);
    for n=1:X.ObjectNumber(length(X.ImageNumber))
        m=X.ObjectNumber(length(X.ImageNumber))*(i-1)+n;
        if X.ImageNumber(m)==i
            n_d = strcat('=NORM.DIST(A',num2str(n+2),',AVERAGE(A:A),STDEV.P(A:A),FALSE)');
            norm_dist=[norm_dist;n_d];
            average=[average;mean(intensity(X.ObjectNumber(length(X.ImageNumber))*(i-1)+1:X.ObjectNumber(length(X.ImageNumber))*i))];
            int=[int;intensity(m)];
        end
    end

    val=table(int,norm_dist,average);

    writetable(val,strcat('base',num2str(i),'SNR.csv'))
    end


%% Intensity

G=readtable(strcat(directory,'\decoding\QT_0.3_0_gene_n_count.csv'));
statistics=table();
for i=1:base_num
    int=intensity(X.ObjectNumber(length(X.ImageNumber))*(i-1)+1:X.ObjectNumber(length(X.ImageNumber))*i);
    Quantile_1st = quantile(int,0.25);
    Quantile_2nd = quantile(int,0.5);
    Quantile_3rd = quantile(int,0.75);
    Average = mean(int);
    Blob_num= length(int);
    stats = table(Blob_num,Average,Quantile_1st, Quantile_2nd, Quantile_3rd);
    for j=1:length(G.GeneName)
    stats.(strcat(char(G.GeneName(j)))) = G.Count(j);
    end 
    statistics=[statistics;stats]
end
writetable(statistics,'Stats.csv')

