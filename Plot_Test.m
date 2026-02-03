%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing

CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [2004 7 1 10 0]; %Time as Year, Month, Day, Hour, Minute
R12 = 64.8;

Source_Site = Chesapeake_Site; % Define the source site for ray tracing

%Generate the Ionospheric Grid
[iono_pf_grid, iono_pf_grid_5, collision_freq, iono_te_grid, start_height,...
    height_inc, range_inc, irreg] = Ionospheric_Grid_2D(Source_Site,...
    Auburn_Site,UT);


%Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Source_Site,Auburn_Site);

%Initialize the Raytrace Model
bearing = azim; %Set the bearing for the transmission
irregs_flag = 1; %Pass the Irregularity Information
nhops = 10; %Only interested in the First Hop
tol = 1e-7; %Set Tolerance (What does this do?)

elevs_init = 1:90;
freq_init = 5; %Set Frequency
freqs_init = freq_init.*ones(size(elevs_init));

%Convert Plasma Frequency Grid to Electron Number Density
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


[ray_data_init, ray_path_data_init, ray_state_vec_init] = ...
          raytrace_2d(origin_lat, origin_long, elevs_init, bearing, freqs_init, nhops, ...
               tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
               collision_freq, start_height, height_inc, range_inc, irreg);

% plot the rays and ionosphere
figure(1)
UT_str = [num2str(UT(3)) '/' num2str(UT(2)) '/' num2str(UT(1)) '  ' ...
          num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') 'UT'];
freq_str = [num2str(freqs_init) 'MHz'];
R12_str = num2str(R12);
lat_str = num2str(origin_lat);
lon_str = num2str(origin_long);
bearing_str = num2str(azim);
fig_str = [UT_str '   ' freq_str '   R12 = ' R12_str '   lat = ' lat_str ...
           ', lon = ' lon_str ', bearing = ' bearing_str];
set(gcf, 'name', fig_str)
start_range = 0;
end_range = 3000;
end_range_idx = fix((end_range-start_range) ./ range_inc) + 1;
start_ht = start_height;
start_ht_idx = 1;
end_ht = 400;
end_ht_idx = fix(end_ht ./ height_inc) + 1;
iono_pf_subgrid = iono_pf_grid(start_ht_idx:end_ht_idx, 1:end_range_idx);
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, ray_path_data_init, 'color', [1, 1, 0.99], ...
    'linewidth', 2);

set(gcf,'units','normal')
pos = get(gcf,'position');
pos(2) = 0.55;
set(gcf,'position', pos)

% uncomment the following to print figure to hi-res ecapsulated postscript
% and PNG files
set(gcf, 'paperorientation', 'portrait')
set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 61 18])
set(gcf, 'papertype', 'a4') 
% print -depsc2 -loose -opengl test.ps 
% print -dpng test.png



