%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing
CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [2004 7 1 10 0]; %Time as Year, Month, Day, Hour, Minute

Site_Selector = 2; % Set to 1 for Chesapeake, Set to 2 For Corpus Christi

if Site_Selector ==1
    Source_Site = Chesapeake_Site;
end
if Site_Selector ==2
    Source_Site = CorpusChristi_Site;
end

%Generate the Ionospheric Grid
[iono_pf_grid, iono_pf_grid_5, collision_freq, iono_te_grid, start_height,...
    height_inc, range_inc, irreg] = Ionospheric_Grid_2D(Source_Site,...
    Auburn_Site,UT);

%Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Source_Site,Auburn_Site);

%Initialize the Raytrace Model
bearing = azim;
irregs_flag = 1;
nhops = 1;
tol = [1e-7 0.025 25];

elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));

%Convert Plasma Frequency Grid to Electron Number Density
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data_init, ray_path_data_init, ray_state_vec_init] = ...
    raytrace_2d(origin_lat, origin_long, elevs_init, bearing, freqs_init, nhops, ...
        tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
        collision_freq, start_height, height_inc, range_inc, irreg);

%Loop Using Already Initialized Ionospheric Model
elevs = 1:0.05:90;
frequencies = 1:0.01:15;  % 1 to 15 MHz

Closest_Apogee = zeros(size(frequencies));
Closest_Distance = zeros(size(frequencies));
closest_group_range = zeros(size(frequencies));
closest_phase_path = zeros(size(frequencies));
closest_geometric_path_length = zeros(size(frequencies));

max_distance_error = 30;  % km

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
        raytrace_2d(origin_lat, origin_long, elevs, bearing, freqs, nhops, ...
            tol, irregs_flag);
    
    % Extract arrays
    ray_labels = [ray_data.ray_label];
    ground_ranges = [ray_data.ground_range];
    apogees = [ray_data.apogee];
    group_ranges = [ray_data.group_range];
    phase_paths = [ray_data.phase_path];
    geometric_path_lengths = [ray_data.geometric_path_length];

    
    % Filter Indices
    valid_idx = (ray_labels == 1) & (abs(ground_ranges - Target_Distance) <= max_distance_error);
    
    if any(valid_idx)
        % Filter to only valid data
        valid_ranges = ground_ranges(valid_idx);
        valid_apogees = apogees(valid_idx);
        valid_group_ranges = group_ranges(valid_idx);
        valid_phase_paths = phase_paths(valid_idx);
        valid_geometric_path_lengths = geometric_path_lengths(valid_idx);
        
        % Find closest distance among valid points
        [~, idx] = min(abs(valid_ranges - Target_Distance));
        
        % Store results
        Closest_Apogee(i) = valid_apogees(idx);
        Closest_Distance(i) = valid_ranges(idx);
        closest_group_range(i) = valid_group_ranges(idx);
        closest_phase_path(i) = valid_phase_paths(idx);
        closest_geometric_path_length(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee(i) = NaN;
        Closest_Distance(i) = NaN;
        closest_group_range(i) = NaN;
        closest_phase_path(i) = NaN;
        closest_geometric_path_length(i) = NaN;
    end
end

% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length
Calc_Height_Group = sqrt((closest_group_range.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase = sqrt((closest_phase_path.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo = sqrt((closest_geometric_path_length.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord
Calc_Height_Group = Calc_Height_Group - s_mat;
Calc_Height_Phase = Calc_Height_Phase - s_mat;
Calc_Height_Geo = Calc_Height_Geo - s_mat;

%% Plot
figure;
hold on;

scatter(frequencies, Closest_Apogee, 'DisplayName', 'Ray Trace Apogee');
scatter(frequencies, real(Calc_Height_Group), 'DisplayName', 'Calculated Height (Group Path)');
scatter(frequencies, real(Calc_Height_Geo), 'DisplayName', 'Calculated Height (Geometric Path)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UT', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time
if Site_Selector == 1
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nChesapeake to Auburn at %s', date_str);
else
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nCorpus Christi to Auburn at %s', date_str);
end
title(title_str);