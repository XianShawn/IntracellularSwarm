% This main function is used to calculate recorded young's modulus using
% linear curve fitting

%% clear
clear all;
close all
clearvars -global

%% initialize
global select_extend1_withdraw2     %
global re_select_roi_N0_L1_C2       % -1 auto, 2 manul 
global show_figure_on1_off0
show_figure_on1_off0=1

re_select_roi_N0_L1_C2=2            % set everything manul
InMain_select_extend1_withdraw2=2   

select_extend1_withdraw2=InMain_select_extend1_withdraw2

% set parameters for both tip and sample
para.indent_data_length=1024;  % length of indent length
para.R=2000                    % tip radius of the probe, in nanometer
para.v_sample=0.5;            % posion ratio of the sample
para.L_sample=0.000003;       % the thickness of the sample

para.v_tip=0.5;               % consider the conflict of both tip and sample
para.E_tip=135;%Gpa  
para.u = 0.00003;

%% load data

%T24
%pa = 'C:\Users\amnl\Desktop\Xian\data_insitu nuclear mechanics data\XianWang\T24\TxtNov9\TxtNov9\';
 pa = 'C:\Users\amnl\Desktop\cancer strain stress\Changhong\tencurves\';
% pa = 'C:\Users\amnl\Desktop\Xian\data_insitu nuclear mechanics data\XianWang\T24\Txtnov27_updated\Txtnov27_updated\';

%RT4
%pa = 'C:\Users\amnl\Desktop\Xian\data_insitu nuclear mechanics data\XianWang\RT4\TxtOct27\';
%pa = 'C:\Users\amnl\Desktop\Xian\data_insitu nuclear mechanics data\XianWang\RT4\TxtNov24\'
%pa = 'C:\Users\amnl\Desktop\Xian\data_insitu nuclear mechanics data\XianWang\RT4\TxtOct31_updated\';

filename='*.txt'
[filename,pa]=uigetfile([pa filename])
%filename = 'indent_RT4_Oct27.Spot7.V45.235M.txt';
% load txt in column

%define a global pathway
global pfn   
FN=dir([pa filename])
pfn=[pa FN.name]


% z_piezo_NM:height in nm, prc_readout: force in pN 
%[z_piezo_NM,prc_readout,paras]=read_indentation_file_brucker(pfn);
[z_piezo_NM,prc_readout,z_tip_NM,paras]=read_indentation_file_brucker2(pfn);

% show data

show_indentation_data(z_piezo_NM,prc_readout);


%% level_indentation_data

% [z_piezo_NM_c,prc_readout_adjusted_c]=level_indentation_data(z_piezo_NM,prc_readout);
% 
% Displacement=z_piezo_NM_c{select_extend1_withdraw2};      % convert and use extend data
% Force=prc_readout_adjusted_c{select_extend1_withdraw2};   % convert and use extend data

[z_piezo_NM_c,prc_readout_adjusted_c,z_tip_NM_c]=level_indentation_data2(z_piezo_NM,prc_readout,z_tip_NM);

Displacement=z_piezo_NM_c{select_extend1_withdraw2}-z_tip_NM_c{select_extend1_withdraw2};      % convert and use extend data
Force=prc_readout_adjusted_c{select_extend1_withdraw2};   % convert and use extend data


%% select select indentation roi, this ROI will be used to calculate Young's modulus

    [sDisplacement,sForce,ind]=manual_select_line_roi(Displacement,Force,'select indentation roi',para.indent_data_length,'brucker');

    [Esample,EL,EH,cfL,gofR2]=fit_youngs_modulus_linear(sDisplacement,sForce,para,0,1)

    
%% calculate viscoelasticity 

    [E_sample2,EL,EH,Y1,cfL,gofR2] = fit_viscoelasticity_nonlinear(sDisplacement, sForce, para, 0,1)

%% save data
% define a struct named youngs and save all the variables related to youngs
    youngs.Esample = Esample;
    youngs.EL = EL;
    youngs.EH = EH;
    youngs.cfL = cfL;
    youngs.gofR2 = gofR2;
    youngs.sDisplacement = sDisplacement;
    youngs.sForce = sForce;
%   save('Youngs30_1',youngs);
    
%     disp('saving to T24');
%     cd('C:\Users\amnl\Desktop\Xian\matlab_code\analyze_data_indent\matlab_data_T24');
    
     disp('saving to RT4');
     cd('C:\Users\amnl\Desktop\Xian\matlab_code\analyze_data_indent\matlab_data_RT4');
    
    index_filename1 = findstr(filename,'Spot');
    if(isempty(index_filename1))
        index_filename1 = findstr(filename,'spot');
    end
    
    index_filename2 = findstr(filename,'V');    
    if(isempty(index_filename2))
        index_filename2 = findstr(filename,'v');
    end
    
    % index_filename2 = index_filename2(2)% this is for November
    filename(index_filename1:(index_filename2+2))
    save([filename(index_filename1:(index_filename2+2)) '.mat'], 'youngs');
    
    cd('C:\Users\amnl\Desktop\Xian\matlab_code\analyze_data_indent')
%% load all the Esample at the same speed
d = dir('*V30.mat');          % only looking for .mat-Files           
Number_mat = length(d);    % number of .mat-Files         
for i=1:Number_mat                    
    load(['Spot' num2str(i) '.V30.mat'])
    if(gofR2<0.95)
        disp(i);
    end
    E30(i) = Esample;
end

%%  load Esample for one spot

i = 8; % difine the spot to load

% load the file separately and save Esmaple to matrix E
load(['spot' num2str(i) '_v30.mat']);
E_single(1) = youngs.Esample;
load(['spot' num2str(i) '_v45.mat']);
E_single(2) = youngs.Esample;
load(['spot' num2str(i) '_v60.mat']);
E_single(3) = youngs.Esample;

figure(1)
bar(E_single, 0.5);
title(i);

