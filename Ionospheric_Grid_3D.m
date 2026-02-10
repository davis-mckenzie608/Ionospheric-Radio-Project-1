function [iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms, geomag_grid_parms] = Ionospheric_Grid_3D(UT)
%UNTITLED3 Summary of this function goes here
%   generates the params. for the area of interest at the given UT
arguments (Input)
    UT
end

arguments (Output)
    iono_pf_grid
    iono_pf_grid_5
    collision_freq
    Bx
    By
    Bz
    iono_grid_parms
    geomag_grid_parms
end
R12 = 64.8;
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 201;           
lat_start = 26.0;
lat_inc = 0.2;
num_lat = 60.0;
lon_start= -99.0;
lon_inc = 1.0;
num_lon = 24.0;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 0.2;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc); 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

doppler_flag = 1; % Set the Doppler flag for grid generation

[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
                     geomag_grid_parms, doppler_flag);
end