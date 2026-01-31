%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing

CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [2004 7 1 10 0]; %Time as Year, Month, Day, Hour, Minute

%Generate the Ionospheric Grid
[iono_en_grid, iono_en_grid_5, collision_freq, iono_te_grid, start_height,...
    height_inc, range_inc, irreg] = Ionospheric_Grid_2D(CorpusChristi_Site,...
    Auburn_Site,UT);


%Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(CorpusChristi_Site,Auburn_Site);

%Initialize the Raytrace Model
bearing = azim; %Set the bearing for the transmission
irregs_flag = 1; %Pass the Irregularity Information
nhops = 1; %Only interested in the First Hop
tol = 1e-6; %Set Tolerance (What does this do?)

freq = 1; %Set Dummy Frequency
elevs = 45; %Set Dummy Elevation

[ray_data, ray_path_data, ray_state_vec] = ...
          raytrace_2d(origin_lat, origin_long, elevs, bearing, freq, nhops, ...
               tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
               collision_freq, start_height, height_inc, range_inc, irreg);


%Seeking Loop