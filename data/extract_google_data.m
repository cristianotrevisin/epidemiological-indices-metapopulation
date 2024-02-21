clc
close all
clear all

% Read from file
t1 = readtable('2020_IT_Region_Mobility_Report.csv');
t2 = readtable("2021_IT_Region_Mobility_Report.csv");
t3 = readtable("2022_IT_Region_Mobility_Report.csv");

% Select Veneto and columns of interest
t1 = t1(ismember(t1.sub_region_1, 'Veneto'),:);
t2 = t2(ismember(t2.sub_region_1, 'Veneto'),:);
t3 = t3(ismember(t3.sub_region_1, 'Veneto'),:);
t1 = t1(:,[6 9 10:15]); t2 = t2(:,[6 9 10:15]); t3 = t3(:,[6 9 10:15]);

% Make province tables
t1_VR = t1(ismember(t1.iso_3166_2_code, 'IT-VR'),:);
t2_VR = t2(ismember(t2.iso_3166_2_code, 'IT-VR'),:);
t3_VR = t3(ismember(t3.iso_3166_2_code, 'IT-VR'),:);
t_VR = [t1_VR;t2_VR;t3_VR]; clear t1_VR t2_VR t3_VR; t_VR(:,1)=[];

t1_VI = t1(ismember(t1.iso_3166_2_code, 'IT-VI'),:);
t2_VI = t2(ismember(t2.iso_3166_2_code, 'IT-VI'),:);
t3_VI = t3(ismember(t3.iso_3166_2_code, 'IT-VI'),:);
t_VI = [t1_VI;t2_VI;t3_VI]; clear t1_VI t2_VI t3_VI; t_VI(:,1)=[];

t1_BL = t1(ismember(t1.iso_3166_2_code, 'IT-BL'),:);
t2_BL = t2(ismember(t2.iso_3166_2_code, 'IT-BL'),:);
t3_BL = t3(ismember(t3.iso_3166_2_code, 'IT-BL'),:);
t_BL = [t1_BL;t2_BL;t3_BL]; clear t1_BL t2_BL t3_BL; t_BL(:,1)=[];

t1_TV = t1(ismember(t1.iso_3166_2_code, 'IT-TV'),:);
t2_TV = t2(ismember(t2.iso_3166_2_code, 'IT-TV'),:);
t3_TV = t3(ismember(t3.iso_3166_2_code, 'IT-TV'),:);
t_TV = [t1_TV;t2_TV;t3_TV]; clear t1_TV t2_TV t3_TV; t_TV(:,1)=[];

t1_VE = t1(ismember(t1.iso_3166_2_code, 'IT-VE'),:);
t2_VE = t2(ismember(t2.iso_3166_2_code, 'IT-VE'),:);
t3_VE = t3(ismember(t3.iso_3166_2_code, 'IT-VE'),:);
t_VE = [t1_VE;t2_VE;t3_VE]; clear t1_VE t2_VE t3_VE; t_VE(:,1)=[];

t1_PD = t1(ismember(t1.iso_3166_2_code, 'IT-PD'),:);
t2_PD = t2(ismember(t2.iso_3166_2_code, 'IT-PD'),:);
t3_PD = t3(ismember(t3.iso_3166_2_code, 'IT-PD'),:);
t_PD = [t1_PD;t2_PD;t3_PD]; clear t1_PD t2_PD t3_PD; t_PD(:,1)=[];

t1_RO = t1(ismember(t1.iso_3166_2_code, 'IT-RO'),:);
t2_RO = t2(ismember(t2.iso_3166_2_code, 'IT-RO'),:);
t3_RO = t3(ismember(t3.iso_3166_2_code, 'IT-RO'),:);
t_RO = [t1_RO;t2_RO;t3_RO]; clear t1_RO t2_RO t3_RO; t_RO(:,1)=[];

% We only consider Workplace Mobility (6th column of remaining tables)
Time_GMD = datetime(2020,02,15):datetime(2022,10,15);
GMD = zeros(7,height(t_VR));

GMD(1,:) = table2array(t_VR(:,6))';
GMD(2,:) = table2array(t_VI(:,6))';
GMD(3,:) = table2array(t_BL(:,6))';
GMD(4,:) = table2array(t_TV(:,6))';
GMD(5,:) = table2array(t_VE(:,6))';
GMD(6,:) = table2array(t_PD(:,6))';
GMD(7,:) = table2array(t_RO(:,6))';

GMD = 1+GMD/100; 

save('google-data.mat','Time_GMD',"GMD")