function [iono_en_grid, iono_en_grid_5, collison_freq, iono_te_grid, start_height, height_inc, range_inc, irreg] = Ionospheric_Grid_2D(Source_Location,Destination_Location,UT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    Source_Location
    Destination_Location
    UT
end

arguments (Output)
    iono_en_grid
    iono_en_grid_5
    collison_freq
    iono_te_grid
    start_height
    height_inc
    range_inc
    irreg
end

kp = 0;
doppler_flag = 1;
R12 = 84.1;
num_range = 501;
num_heights = 500;
max_range = 10000; % maximum group range for data grids (Km) (SEE IMPORTANT NOTE ABOVE)
start_height = 0; % desired start height of iono_pf_grid & bfield grids (km) (SEE IMPORTANT NOTE ABOVE)
range_inc = max_range ./ (num_range - 1); % desired range step of iono_pf_grid, bfield, iono_parms, dec, dip, and irreg_strngth arrays (km) (SEE IMPORTANT NOTE ABOVE)
height_inc = 3; % desired height step of iono_pf_grid & bfield grids (km) (SEE IMPORTANT NOTE ABOVE)


[azim, Target_Distance, origin_lat, origin_lon] =...
    Bearing_Calculator(Source_Location,Destination_Location);       

[iono_en_grid, iono_en_grid_5, collison_freq, irreg, iono_te_grid] = ...
          gen_iono_grid_2d(origin_lat, origin_lon, R12, UT, azim, max_range, ...
                           num_range, range_inc, start_height, height_inc, ...
                           num_heights, kp, doppler_flag);


end