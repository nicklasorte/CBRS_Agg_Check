clear;
clc;
close all;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
top_start_clock=clock;
folder1='C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\3.5GHz ITM Reliability Webster';
cd(folder1)
addpath(folder1)
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Basic_Functions')  %%%%%%Another Github repo
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\General_Terrestrial_Pathloss')  %%%%%%%%This is another Github repo
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\General_Movelist')  %%%%%%%%This is another Github repo
addpath('C:\Local Matlab Data\Local MAT Data')  %%%%%%%One Drive Error and need these files
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sims on Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drive

%%%%%%%Loading  DPA Geography
load('cell_expand_all_dpa.mat','cell_expand_all_dpa')
load('cell_portal_dpa_data.mat','cell_portal_dpa_data')
load('downsampled_east10km.mat','downsampled_east10km')
load('downsampled_west10km.mat','downsampled_west10km')
load('us_cont.mat','us_cont')
load('cell_e_dpa_data.mat','cell_e_dpa_data')%%%%%%%%Pull in the Neighborhoods
%%%%%%%%%%%%%%%%%%%cell_e_dpa_data
    % %1)Name,
    % % 2)Lat
    % %%3)Lon,
    % % 4) Radar Threshold,
    % % 5) Radar Height,
    % % 6)Radar Beamwidth,
    % % 7) Min Azimuth
    % % 8) Max Azimuth
    %%%9) CatB Above 6m dist
    %%%10) CatB Below 6m dist
    %%%11) CatA inside Above 6m dist
    %%%12) CatA inside Below 6m dist
    %%%13) Low Freq
    %%%14) High Freq
    %%%15) Cell Geometry
    %%%16) CatA Outdoor Above 6m dist
    %%%17) CatA Outdoor Below 6m dist


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
load('Cascade_new_full_census_2010.mat','new_full_census_2010')%%%%%%%Geo Id, Center Lat, Center Lon,  NLCD (1-4), Population
load('cell_census_centroid.mat','cell_census_centroid')
cell_bs_data=cell_census_centroid;
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev=111; %%%%%%East and West
freq_separation=0; %%%%%%%Assuming co-channel
array_bs_eirp=horzcat(47,47,47); %%%%%EIRP [dBm/100MHz] for Rural, Suburan, Urban 
network_loading_reduction=8
fdr_dB=0  %%%%%%%%%%%Loss from FDR
pol_mismatch=0 %%%%%%%%%Polarization mismatch dB
mitigation_dB=0%vertcat(0,30); %%%%% Beam Muting or PRB Blanking (or any other mitigation mechanism):  30 dB reduction %%%%%%%%%%%%Consider have this be an array, 3dB step size, to get a more granular insight into how each 3dB mitigation reduces the coordination zone.
tf_opt=0; %%%%This is for the optimized move list, (not WinnForum)
margin=1%2; %%%2dB %dB margin for aggregate interference
reliability=[1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,91,92,93,94,95,96,97,98,99]'; %%%A custom ITM range to interpolate from
move_list_reliability=50;
agg_check_reliability=reliability;
building_loss=15;
Tpol=1; %%%polarization for ITM
FreqMHz=3550;
confidence=50;
mc_percentile=100
move_list_mc_percentile=mc_percentile
mc_size=1
move_list_mc_size=mc_size
agg_check_mc_percentile=95;
agg_check_mc_size=1000;
deployment_percentage=100;
%%%%sim_radius_km=200; %%%%%%%%Change Sim_radius to the [neighborhood distance x 1.5] and make it unique to each DPA::::        binary_dist_array=[2,4,8,16,32,64,128,256,512,1024,2048];
tf_census=1; %%%%%%%If 0, then us randomized real.
base_station_height=30%NaN(1,1); %%%%If NaN, then keep the normal heights
tf_cbsd_mask=1;
tf_clutter=0;%1;  %%%%%%This if for P2108.
num_pts=8%%%%8;  %%%%%%%Number of DPA Sample Points along the front edge, Typically 8
number_rand_pts=NaN(1);
min_ant_loss=40; %%%%%%%%Main to side gain: 40dB
dpa_idx=1:1:find(contains(cell_expand_all_dpa(:,1),'BremertonEverett'))
dpa_data_idx=1:1:find(contains(cell_e_dpa_data(:,1),'Bremerton'))
sim_folder1='Z:\Matlab2025 Sims\3.5GHz CBRS Reliability Sims'%%%'C:\Local Matlab Data\3.5GHz Reliability Sims'%%'Z:\Matlab2025 Sims\3.5GHz CBRS Reliability Sims'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
horzcat(cell_expand_all_dpa(dpa_idx,1),cell_e_dpa_data(dpa_data_idx,1))
% % % 
% % % 'check'
% % % pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tf_census==0
    %%%Just select Census Lat/Lon or the Randomized Real
    tic;
    disp_progress(app,'Loading Randomized Real . . .')
    %%load('cell_err_data.mat','cell_err_data') %%%%%%%%Placeholder of the 5G deployment
    load('cell_err_data_single_sector.mat','cell_err_data_single_sector')
    cell_err_data=cell_err_data_single_sector;
    toc; %%%%%%%15 Seconds
    cell_bs_data=cell_err_data;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
load('aas_zero_elevation_data.mat','aas_zero_elevation_data')
toc;
%%%%1) Azimuth -180~~180
%%%2) Rural
%%%3) Suburban
%%%4) Urban
%%%%AAS Reduction in Gain to Max Gain (0dB is 0dB reduction, which equates to the make antenna gain of 25dB)
%%%%Need to normalize to zero after the "downtilt reductions" are calculated
%%%%To simplify the data, this is gain at the horizon. 50th Percentile

aas_zero_elevation_data(1,:)
%%%%%%%%%%Set all gains to 0dB, since there will be no azimuths
aas_zero_elevation_data(:,[2:4])=0;
aas_zero_elevation_data(1,:)
bs_down_tilt_reduction=abs(max(aas_zero_elevation_data(:,[2:4]))) %%%%%%%%Downtilt dB Value for Rural/Suburban/Urban
norm_aas_zero_elevation_data=horzcat(aas_zero_elevation_data(:,1),aas_zero_elevation_data(:,[2:4])+bs_down_tilt_reduction);
max(norm_aas_zero_elevation_data(:,[2:4])) %%%%%This should be [0 0 0]
array_bs_eirp_reductions=(array_bs_eirp-bs_down_tilt_reduction-network_loading_reduction-fdr_dB-pol_mismatch)-mitigation_dB %%%%%Rural, Suburban, Urban cols:(1-3), No Mitigations/Mitigations rows:(1-2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%1)Name, 
% %%%%%%2)Lat/Lon, 
% %%%%%%3) Radar Threshold, 
% %%%%%%4) Radar Height, 
% %%%%%%5)Radar Beamwidth, 
% %%%%%%6)tf_ship (for protection points)
% %%%%%%7)min azimuth
% %%%%%%8)max azimuth
% %%%%%%9)CatB Neighborhood Distance 2.0
cell_dpa_geo=vertcat(cell_expand_all_dpa([dpa_idx],[1,3])); %%%%%Selecting the DPAs to calculate the Neighborhood: Name/Geo Points
cell_threshold=cell(1,1);
cell_dpa_geo(:,3)=cell_e_dpa_data(dpa_idx,4);
cell_dpa_geo(:,4)=cell_e_dpa_data(dpa_idx,5);
cell_dpa_geo(:,5)=cell_e_dpa_data(dpa_idx,6);
cell_tf_ship1=cell(1,1);
cell_tf_ship1{1}=1;
cell_tf_ship0=cell(1,1);
cell_tf_ship0{1}=0;
for i=1:1:length(dpa_idx)
    if contains(cell_dpa_geo(i,1),'East') || contains(cell_dpa_geo(i,1),'West')
        cell_dpa_geo(i,6)=cell_tf_ship1; 
    else
         cell_dpa_geo(i,6)=cell_tf_ship0; 
    end
end
cell_dpa_geo(:,7)=cell_e_dpa_data(dpa_idx,7);
cell_dpa_geo(:,8)=cell_e_dpa_data(dpa_idx,8);
cell_dpa_geo(:,9)=cell_e_dpa_data(dpa_idx,9);
cell_dpa_geo'



%cell_e_dpa_data(dpa_idx,:)'
%%%%%%%%%%%%%%%%%%%cell_e_dpa_data
    % %1)Name,
    % % 2)Lat
    % %%3)Lon,
    % % 4) Radar Threshold,
    % % 5) Radar Height,
    % % 6)Radar Beamwidth,
    % % 7) Min Azimuth
    % % 8) Max Azimuth
    %%%9) CatB Above 6m dist
    %%%10) CatB Below 6m dist
    %%%11) CatA inside Above 6m dist
    %%%12) CatA inside Below 6m dist
    %%%13) Low Freq
    %%%14) High Freq
    %%%15) Cell Geometry
    %%%16) CatA Outdoor Above 6m dist
    %%%17) CatA Outdoor Below 6m dist


%%%%%%Create a Rev Folder
cd(sim_folder1);
pause(0.1)
tempfolder=strcat('Rev',num2str(rev));
[status,msg,msgID]=mkdir(tempfolder);
rev_folder=fullfile(sim_folder1,tempfolder);
cd(rev_folder)
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Saving the simulation files in a folder for the option to run from a server
tic;
save('cell_dpa_geo.mat','cell_dpa_geo')
save('reliability.mat','reliability')
save('move_list_reliability.mat','move_list_reliability')
save('confidence.mat','confidence')
save('FreqMHz.mat','FreqMHz')
save('Tpol.mat','Tpol')
save('building_loss.mat','building_loss')
save('mc_percentile.mat','mc_percentile')
save('move_list_mc_percentile.mat','move_list_mc_percentile')
save('move_list_mc_size.mat','move_list_mc_size')
save('agg_check_mc_percentile.mat','agg_check_mc_percentile')
save('agg_check_mc_size.mat','agg_check_mc_size')
save('mc_size.mat','mc_size')
save('margin.mat','margin')
save('deployment_percentage.mat','deployment_percentage')
%save('sim_radius_km.mat','sim_radius_km')
save('array_bs_eirp_reductions.mat','array_bs_eirp_reductions') %%%%%Rural, Suburban, Urban cols:(1-3), No Mitigations/Mitigations rows:(1-2)
save('norm_aas_zero_elevation_data.mat','norm_aas_zero_elevation_data')
save('agg_check_reliability.mat','agg_check_reliability')
save('tf_clutter.mat','tf_clutter')
save('mitigation_dB.mat','mitigation_dB')
save('tf_opt.mat','tf_opt')
toc;




%%%%%%%%%%%%%
[num_loc,~]=size(cell_dpa_geo);
location_table=table([1:1:num_loc]',cell_dpa_geo(:,1))
array_bs_latlon=cell2mat(cell_bs_data(:,[5,6]));
size(array_bs_latlon)

for base_idx=1:1:num_loc
    close all;
    strcat(num2str(base_idx/num_loc*100),'%')
    
    temp_cell_geo_data=cell_dpa_geo(base_idx,:)
    data_label1=erase(temp_cell_geo_data{1}," ");  %%%Remove the Spaces
    
    %%%%%%%%%Step 1: Make a Folder for this single DPA
    cd(rev_folder);
    pause(0.1)
    tempfolder2=strcat(data_label1);
    [status,msg,msgID]=mkdir(tempfolder2);
    sim_folder=fullfile(rev_folder,tempfolder2);
    cd(sim_folder)
    pause(0.1)
    
    
    base_polygon=temp_cell_geo_data{2};  %%%%%%DPA or Base
    save(strcat(data_label1,'_base_polygon.mat'),'base_polygon')
    
    tf_ship=temp_cell_geo_data{6};
    if tf_ship==1
        %%%%Find the inner edge
        uni_base_polygon=unique(base_polygon,'stable','rows');
        inner_edge=vertcat(downsampled_east10km,downsampled_west10km);
        [inner_line,inner_corner1,inner_corner2]=find_dpa_line_overlap(inner_edge,uni_base_polygon);
        base_protection_pts=curvspace(inner_line,num_pts);
    else
        [num_ppts,~]=size(base_polygon);
        if num_ppts==1
            %Do nothing
            base_protection_pts=base_polygon;
        else
            temp_pts=curvspace(base_polygon,num_pts+1);
            base_protection_pts=temp_pts([1:num_pts],:);
            figure;
            hold on;
            plot(base_polygon(:,2),base_polygon(:,1),'-k')
            plot(base_protection_pts(:,2),base_protection_pts(:,1),'or','LineWidth',2)
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%Create pp_pts
    [num_ppts,~]=size(base_polygon);
    if num_ppts==1
        %Do nothing
    else
        if ~isnan(number_rand_pts)==1
            %%%%%Uniform Random Points
            %%%%%%DPA Bounds
            x_max=max(base_polygon(:,2));
            x_min=min(base_polygon(:,2));
            y_max=max(base_polygon(:,1));
            y_min=min(base_polygon(:,1));
            
            rng(rev);%For Repeatability
            %%%Preallocate
            marker1=1;
            rand_pts_uni=NaN(number_rand_pts,2);
            while (marker1<=number_rand_pts)
                %%%Generate Random Points inside DPA
                x_rand=rand(1);
                y_rand=rand(1);
                
                x_pt=x_rand*(x_max-x_min)+x_min;
                y_pt=y_rand*(y_max-y_min)+y_min;
                
                %%%%Check to see if it falls inside the DPA
                tf1=inpolygon(x_pt,y_pt,base_polygon(:,2),base_polygon(:,1));
                
                if tf1==1
                    rand_pts_uni(marker1,:)=horzcat(y_pt,x_pt);
                    marker1=marker1+1;
                    
                end
            end
            
            close all;
            figure;
            hold on;
            scatter(rand_pts_uni(:,2),rand_pts_uni(:,1),10,'or')
            scatter(base_protection_pts(:,2),base_protection_pts(:,1),10,'db')
            plot(base_polygon(:,2),base_polygon(:,1),'-k','LineWidth',2)
            grid on;
            base_protection_pts=vertcat(base_protection_pts,rand_pts_uni);
            size(base_protection_pts)
            pause(0.1)
        end
    end
    radar_height=temp_cell_geo_data{4};
    base_protection_pts(:,3)=radar_height;
    save(strcat(data_label1,'_base_protection_pts.mat'),'base_protection_pts') %%%%%Save the Protection Points
    

    radar_threshold=temp_cell_geo_data{3};
    save(strcat(data_label1,'_radar_threshold.mat'),'radar_threshold')
    radar_beamwidth=temp_cell_geo_data{5};
    save(strcat(data_label1,'_radar_beamwidth.mat'),'radar_beamwidth')
    save(strcat(data_label1,'_min_ant_loss.mat'),'min_ant_loss')
    min_azimuth=temp_cell_geo_data{7};
    save(strcat(data_label1,'_min_azimuth.mat'),'min_azimuth')
    max_azimuth=temp_cell_geo_data{8};
    save(strcat(data_label1,'_max_azimuth.mat'),'max_azimuth')
    neighborhood_dist=temp_cell_geo_data{9};
    save(strcat(data_label1,'_neighborhood_dist.mat'),'neighborhood_dist')


    figure;
    hold on;
    plot(base_polygon(:,2),base_polygon(:,1),'-r')
    plot(base_protection_pts(:,2),base_protection_pts(:,1),'ok')
    grid on;
    size(base_protection_pts)
    plot_google_map('maptype','terrain','APIKey','AIzaSyCgnWnM3NMYbWe7N4svoOXE7B2jwIv28F8') %%%Google's API key made by nick.matlab.error@gmail.com
    filename1=strcat('Operational_Area_',data_label1,'.png');
    pause(0.1)
    saveas(gcf,char(filename1))


    %%%%%%%%Sim Bound
    if any(isnan(base_polygon))
        base_polygon=base_polygon(~isnan(base_polygon(:,1)),:);

        figure;
        plot(base_polygon(:,2),base_polygon(:,1),'-r')
    end
    sim_radius_km=neighborhood_dist*1.1;
    [sim_bound]=calc_sim_bound(app,base_polygon,sim_radius_km,data_label1);

    %%%%%%%Filter Base Stations that are within sim_bound
    tic;
    bs_inside_idx=find(inpolygon(array_bs_latlon(:,2),array_bs_latlon(:,1),sim_bound(:,2),sim_bound(:,1))); %Check to see if the points are in the polygon
    toc;
    size(bs_inside_idx)
    temp_sim_cell_bs_data=cell_bs_data(bs_inside_idx,:);


    %%%%%%%%%%%%Downsample deployment
    [num_inside,~]=size(bs_inside_idx)
    sample_num=ceil(num_inside*deployment_percentage/100)
    rng(rev+base_idx); %%%%%%%For Repeatibility
    rand_sample_idx=datasample(1:num_inside,sample_num,'Replace',false);
    size(temp_sim_cell_bs_data)
    temp_sim_cell_bs_data=temp_sim_cell_bs_data(rand_sample_idx,:);
    size(temp_sim_cell_bs_data)
    temp_lat_lon=cell2mat(temp_sim_cell_bs_data(:,[5,6]));


    figure;
    hold on;
    plot(temp_lat_lon(:,2),temp_lat_lon(:,1),'ob')
    plot(sim_bound(:,2),sim_bound(:,1),'-r','LineWidth',3)
    plot(base_protection_pts(:,2),base_protection_pts(:,1),'sr','Linewidth',4)
    grid on;
    plot_google_map('maptype','terrain','APIKey','AIzaSyCgnWnM3NMYbWe7N4svoOXE7B2jwIv28F8') %%%Google's API key made by nick.matlab.error@gmail.com
    filename1=strcat('Sim_Area_Deployment_',data_label1,'.png');
    pause(0.1)
    saveas(gcf,char(filename1))

    %%%%%%%%%%Add an index for R/S/U (NLCD)
    rural_idx=find(contains(temp_sim_cell_bs_data(:,11),'R'));
    sub_idx=find(contains(temp_sim_cell_bs_data(:,11),'S'));
    urban_idx=find(contains(temp_sim_cell_bs_data(:,11),'U'));
    [num_bs,num_col]=size(temp_sim_cell_bs_data);
    array_ncld_idx=NaN(num_bs,1);
    array_ncld_idx(rural_idx)=1;
    array_ncld_idx(sub_idx)=2;
    array_ncld_idx(urban_idx)=3;
    cell_ncld=num2cell(array_ncld_idx);

   
    clutter_dB=zeros(num_bs,1);  %%%%%%Clutter is done somewhere else.
    array_eirp_bs=NaN(num_bs,2); %%%%%1)No Mitigations, 2)Mitigations --> 14 and 15 of cell
    for i=1:1:num_bs
        temp_nlcd_idx=array_ncld_idx(i);
        array_eirp_bs(i,:)=array_bs_eirp_reductions(:,temp_nlcd_idx)-clutter_dB(i);
    end
    cell_eirp1=num2cell(array_eirp_bs(:,1));
    cell_eirp2=cell_eirp1;
    sim_cell_bs_data=horzcat(temp_sim_cell_bs_data,cell_ncld,cell_eirp1,cell_eirp2);
    size(sim_cell_bs_data)

    sim_cell_bs_data(1,:)

    %%%1) LaydownID
    %%%2) FCCLicenseID
    %%%3) SiteID
    %%%4) SectorID
    %%%5) SiteLatitude_decDeg
    %%%6) SiteLongitude_decDeg
    %%%7) SE_BearingAngle_deg
    %%%8) SE_AntennaAzBeamwidth_deg
    %%%9) SE_DownTilt_deg  %%%%%%%%%%%%%%%%%(Check for Blank)
    %%%10) SE_AntennaHeight_m
    %%%11) SE_Morphology
    %%%12) SE_CatAB
    %%%%%%%%%%13) NLCD idx
    %%%%%%%%%14) EIRP (no mitigations)
    %%%%%%%%%15) EIRP (mitigations)

    tic;
    save(strcat(data_label1,'_sim_cell_bs_data.mat'),'sim_cell_bs_data')
    toc; %%%%%%%%%3 seconds


    %%%%%%%%%%%%%%%Also include the array of the list_catb (order) that we
    %%%%%%%%%%%%%%%usually use for the other sims. (As this will be used
    %%%%%%%%%%%%%%%for the path loss and move list.)

    sim_cell_bs_data(1,:)
    [num_tx,~]=size(sim_cell_bs_data)

    sim_array_list_bs=horzcat(cell2mat(sim_cell_bs_data(:,[5,6,10,14])),NaN(num_tx,1),array_ncld_idx,cell2mat(sim_cell_bs_data(:,[7,15])));
    [num_bs_sectors,~]=size(sim_array_list_bs);
    sim_array_list_bs(:,5)=1:1:num_bs_sectors;
    % % %      %%%%array_list_bs  %%%%%%%1) Lat, 2)Lon, 3)BS height, 4)BS EIRP Adjusted 5) Nick Unique ID for each sector, 6)NLCD: R==1/S==2/U==3, 7) Azimuth 8)BS EIRP Mitigation
    %%%%%%%%If there is no mitigation EIRPs, make all of these NaNs (column 8)

    %%%%%%%%%%%Put the rest of the Link Budget Parameters in this list

    %%%%%%%9) EIRP dBm:         array_bs_eirp
    sim_array_list_bs(rural_idx,9)=array_bs_eirp(1);
    sim_array_list_bs(sub_idx,9)=array_bs_eirp(2);
    sim_array_list_bs(urban_idx,9)=array_bs_eirp(3);

    %%%%%%%10) AAS (Vertical) dB Reduction: (Downtilt)   %%%%%%%%Downtilt dB Value for Rural/Suburban/Urban
    sim_array_list_bs(rural_idx,10)=bs_down_tilt_reduction(1);
    sim_array_list_bs(sub_idx,10)=bs_down_tilt_reduction(2);
    sim_array_list_bs(urban_idx,10)=bs_down_tilt_reduction(3);

    %%%%%%%%%11)Clutter
    sim_array_list_bs(:,11)=clutter_dB;

    %%%%%%%%%12)Network Loading and TDD (dB)
    sim_array_list_bs(:,12)=network_loading_reduction;

    %%%%%%%%%13)FDR (dB)
    sim_array_list_bs(:,13)=fdr_dB;

    %%%%%%%%%14)Polarization (dB)
    sim_array_list_bs(:,14)=pol_mismatch;


    %%%%%%%%%15)Mitigation Reduction (dB) (Not doing mitigations this way)
    sim_array_list_bs(:,15)=0;

    sim_array_list_bs(1,:)
    size(sim_array_list_bs)


    if ~isnan(base_station_height)
        sim_array_list_bs(:,3)=base_station_height;
        'Change all BS height'
        unique(sim_array_list_bs(:,3))
        %pause;
    end
    unique(sim_array_list_bs(:,3))

    tic;
    save(strcat(data_label1,'_sim_array_list_bs.mat'),'sim_array_list_bs')
    toc; %%%%%%%%%3 seconds
        % % %      %%%%array_list_bs  %%%%%%%1) Lat, 2)Lon, 3)BS height, 4)BS EIRP Adjusted 5) Nick Unique ID for each sector, 6)NLCD: R==1/S==2/U==3, 7) Azimuth 8)BS EIRP Mitigation


        'Check for nans in power'
        unique(sim_array_list_bs(:,4))
        any(isnan(sim_array_list_bs(:,4)))
        unique(sim_array_list_bs(:,3))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parallel_flag=1%0%1 %%%%%0 --> serial, 1 --> parallel (if you have parallel toolbox)
[workers,parallel_flag]=check_parallel_toolbox(app,parallel_flag)
tf_server_status=0

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Where the magic happens
wrapper_move_list_agg_check_rev2(app,parallel_flag,rev_folder,tf_server_status,workers)

end_clock=clock;
total_clock=end_clock-top_start_clock;
total_seconds=total_clock(6)+total_clock(5)*60+total_clock(4)*3600+total_clock(3)*86400;
total_mins=total_seconds/60;
total_hours=total_mins/60;
if total_hours>1
    strcat('Total Hours:',num2str(total_hours))
elseif total_mins>1
    strcat('Total Minutes:',num2str(total_mins))
else
    strcat('Total Seconds:',num2str(total_seconds))
end
cd(folder1)
'Done'


    








