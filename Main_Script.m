%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing

CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [1969 7 20 18 0]; %Time as Year, Month, Day, Hour, Minute

Source_Site = CorpusChristi_Site; % Define the source site for ray tracing

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
nhops = 1; %Only interested in the First Hop
tol = 1e-7; %Set Tolerance (What does this do?)

elevs_init = 20;
freq_init = 5; %Set Frequency
freqs_init = freq_init.*ones(size(elevs_init));

%Convert Plasma Frequency Grid to Electron Number Density
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


[ray_data_init, ray_path_data_init, ray_state_vec_init] = ...
          raytrace_2d(origin_lat, origin_long, elevs_init, bearing, freqs_init, nhops, ...
               tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
               collision_freq, start_height, height_inc, range_inc, irreg);

%Loop Using Already Initialized Ionospheric Model

elevs = 1:0.5:90;
frequencies = 1:0.1:15;  % 1 to 20 MHz

Closest_Apogee = zeros(size(frequencies));
Closest_Distance = zeros(size(frequencies));
Closest_Elevation = zeros(size(frequencies));

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    
    [ray_data, ray_path_data, ray_state_vec] = ...
        raytrace_2d(origin_lat, origin_long, elevs, bearing, freqs, nhops, ...
            tol, irregs_flag);
    
    % Extract arrays
    ground_ranges = [ray_data.ground_range];
    apogees = [ray_data.apogee];
    
    % Find valid indices (where apogee is a real number)
    valid_idx = ~isnan(apogees) & ~isinf(apogees);
    
    if any(valid_idx)
        % Filter to only valid data
        valid_ranges = ground_ranges(valid_idx);
        valid_apogees = apogees(valid_idx);
        valid_elevations = [ray_data(valid_idx).initial_elev];
        
        % Find closest distance among valid points
        [~, idx] = min(abs(valid_ranges - Target_Distance));
        
        % Store results
        Closest_Apogee(i) = valid_apogees(idx);
        Closest_Distance(i) = valid_ranges(idx);
        Closest_Elevation(i) = valid_elevations(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee(i) = NaN;
        Closest_Distance(i) = NaN;
        Closest_Elevation(i) = NaN;
        warning('No valid rays found for %.1f MHz', freq);
    end
end
