%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing

tic 

CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [2004 7 1 10 0]; %Time as Year, Month, Day, Hour, Minute

Site_Selector = 1; % Set to 1 for Chesapeake, Set to 2 For Corpus Christi
Simulation_Selector = 2; %Set to 1 for 2D Simulation, 2 For 3D Simulation, 3 for Chesapeake Anisotropy, 4 for Corpus Christi Anisotropy

if Site_Selector ==1
    Source_Site = Chesapeake_Site;
end
if Site_Selector ==2
    Source_Site = CorpusChristi_Site;
end

if Simulation_Selector ==1;
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
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee = zeros(size(frequencies));
Closest_Distance = zeros(size(frequencies));
closest_group_range = zeros(size(frequencies));
closest_phase_path = zeros(size(frequencies));
closest_geometric_path_length = zeros(size(frequencies));
closest_lats = zeros(size(frequencies));
closest_lons = zeros(size(frequencies));

max_distance_error = 5;  % km

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
    lats = [ray_data.lat];
    lons = [ray_data.lon];

    
    % Filter Indices
    valid_idx = (ray_labels == 1) & (abs(ground_ranges - Target_Distance) <= max_distance_error);
    
    if any(valid_idx)
        % Filter to only valid data
        valid_ranges = ground_ranges(valid_idx);
        valid_apogees = apogees(valid_idx);
        valid_group_ranges = group_ranges(valid_idx);
        valid_phase_paths = phase_paths(valid_idx);
        valid_geometric_path_lengths = geometric_path_lengths(valid_idx);
        valid_lats = lats(valid_idx);
        valid_lons = lons(valid_idx);
        
        % Find closest distance among valid points
        [~, idx] = min(abs(valid_ranges - Target_Distance));
        
        % Store results
        Closest_Apogee(i) = valid_apogees(idx);
        Closest_Distance(i) = valid_ranges(idx);
        closest_group_range(i) = valid_group_ranges(idx);
        closest_phase_path(i) = valid_phase_paths(idx);
        closest_geometric_path_length(i) = valid_geometric_path_lengths(idx);
        closest_lats(i) = valid_lats(idx);
        closest_lons(i) = valid_lons(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee(i) = NaN;
        Closest_Distance(i) = NaN;
        closest_group_range(i) = NaN;
        closest_phase_path(i) = NaN;
        closest_geometric_path_length(i) = NaN;
        closest_lats(i) = NaN;
        closest_lons(i) = NaN;
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
fig = figure;
hold on;

p{1} = scatter(frequencies, Closest_Apogee, 'DisplayName', 'Ray Trace Apogee');
p{2} = scatter(frequencies, real(Calc_Height_Group), 'DisplayName', 'Calculated Height (Group Path)');
p{3} = scatter(frequencies, real(Calc_Height_Geo), 'DisplayName', 'Calculated Height (Geometric Path)');
%p{4} = geoscatter(closest_lats,closest_lons,[],"m","d");

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time
if Site_Selector == 1
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nChesapeake to Auburn at %s', date_str);
else
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nCorpus Christi to Auburn at %s', date_str);
end
title(title_str);
ax = gca;
lgd = legend;
IEEE_Format_Plot_v2('figure', fig, 'whr', 1, 'scale', 2, 'curve', p, 'axis', ax, 'legend', lgd);

 
geo1 = figure;

geoscatter(closest_lats,closest_lons,[],"blue","x");


hold on;
geoscatter(32.604722, 360-85.486667,[],"red","x");
hold off;
if Site_Selector == 1 %Chesapeake Title
    title(sprintf('Landing Point of Traced Rays from\nChesapeake, VA to Auburn, AL at %s', date_str));
else
    title(sprintf('Landing Point of Traced Rays from\nCorpus Christi, TX to Auburn, AL at %s', date_str));
end
legend('Ray Landing Points','Broun Hall','Location','best');
gx = gca; 
gx.LongitudeLabel.FontSize = 18; 
gx.LatitudeLabel.FontSize = 18;
gx.Title.FontSize = 18;
gx.FontName = 'Times New Roman';
gx.MapCenter = [32.604722, 360-85.486667];
gx.Basemap = 'topographic';
gx.GridLineStyle = '--';
gx.GridAlpha = 0.4;
gx.Legend.FontSize = 12;



end

if Simulation_Selector == 2;
% Constants
origin_ht = 0.0;

% Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Source_Site,Auburn_Site);
% Generate the Ionospheric Parameters
[iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms,...
    geomag_grid_parms] = Ionospheric_Grid_3D(UT);
% Initialize Ionospheric Model
elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));
ray_bears = ones(size(elevs_init)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs_init, ray_bears,...
                        freqs_init, OX_mode, nhops, tol, iono_en_grid, ...
     	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
                        Bx, By, Bz, geomag_grid_parms);
%% O-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_O = zeros(size(frequencies));
Closest_Distance_O = zeros(size(frequencies));
closest_group_range_O = zeros(size(frequencies));
closest_phase_path_O = zeros(size(frequencies));
closest_geometric_path_length_O = zeros(size(frequencies));
closest_lats_O = zeros(size(frequencies));
closest_lons_O = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 3; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
    % Extract arrays
    ray_labels = [ray_data.ray_label];
    ground_ranges = [ray_data.ground_range];
    apogees = [ray_data.apogee];
    group_ranges = [ray_data.group_range];
    phase_paths = [ray_data.phase_path];
    geometric_path_lengths = [ray_data.geometric_path_length];
    lats_O = [ray_data.lat];
    lons_O = [ray_data.lon];

    
    % Filter Indices
    valid_idx = (ray_labels == 1) & (abs(ground_ranges - Target_Distance) <= max_distance_error);
    
    if any(valid_idx)
        % Filter to only valid data
        valid_ranges = ground_ranges(valid_idx);
        valid_apogees = apogees(valid_idx);
        valid_group_ranges = group_ranges(valid_idx);
        valid_phase_paths = phase_paths(valid_idx);
        valid_geometric_path_lengths = geometric_path_lengths(valid_idx);
        valid_lats_O = lats_O(valid_idx);
        valid_lons_O = lons_O(valid_idx);
        
        % Find closest distance among valid points
        [~, idx] = min(abs(valid_ranges - Target_Distance));
        
        % Store results
        Closest_Apogee_O(i) = valid_apogees(idx);
        Closest_Distance_O(i) = valid_ranges(idx);
        closest_group_range_O(i) = valid_group_ranges(idx);
        closest_phase_path_O(i) = valid_phase_paths(idx);
        closest_geometric_path_length_O(i) = valid_geometric_path_lengths(idx);
        closest_lats_O(i) = valid_lats_O(idx);
        closest_lons_O(i) = valid_lons_O(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_O(i) = NaN;
        Closest_Distance_O(i) = NaN;
        closest_group_range_O(i) = NaN;
        closest_phase_path_O(i) = NaN;
        closest_geometric_path_length_O(i) = NaN;
        closest_lats_O(i) = NaN;
        closest_lons_O(i) = NaN;
    end
end
%% X-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_X = zeros(size(frequencies));
Closest_Distance_X = zeros(size(frequencies));
closest_group_range_X = zeros(size(frequencies));
closest_phase_path_X = zeros(size(frequencies));
closest_geometric_path_length_X = zeros(size(frequencies));
closest_lats_X = zeros(size(frequencies));
closest_lons_X = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = -1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 3; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
    % Extract arrays
    ray_labels = [ray_data.ray_label];
    ground_ranges = [ray_data.ground_range];
    apogees = [ray_data.apogee];
    group_ranges = [ray_data.group_range];
    phase_paths = [ray_data.phase_path];
    geometric_path_lengths = [ray_data.geometric_path_length];
    lats_X = [ray_data.lat];
    lons_X = [ray_data.lon];

    
    % Filter Indices
    valid_idx = (ray_labels == 1) & (abs(ground_ranges - Target_Distance) <= max_distance_error);
    
    if any(valid_idx)
        % Filter to only valid data
        valid_ranges = ground_ranges(valid_idx);
        valid_apogees = apogees(valid_idx);
        valid_group_ranges = group_ranges(valid_idx);
        valid_phase_paths = phase_paths(valid_idx);
        valid_geometric_path_lengths = geometric_path_lengths(valid_idx);
        valid_lats_X = lats_X(valid_idx);
        valid_lons_X = lons_X(valid_idx);
        
        % Find closest distance among valid points
        [~, idx] = min(abs(valid_ranges - Target_Distance));
        
        % Store results
        Closest_Apogee_X(i) = valid_apogees(idx);
        Closest_Distance_X(i) = valid_ranges(idx);
        closest_group_range_X(i) = valid_group_ranges(idx);
        closest_phase_path_X(i) = valid_phase_paths(idx);
        closest_geometric_path_length_X(i) = valid_geometric_path_lengths(idx);
        closest_lats_X(i) = valid_lats_X(idx);
        closest_lons_X(i) = valid_lons_X(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_X(i) = NaN;
        Closest_Distance_X(i) = NaN;
        closest_group_range_X(i) = NaN;
        closest_phase_path_X(i) = NaN;
        closest_geometric_path_length_X(i) = NaN;
        closest_lats_X(i) = NaN;
        closest_lons_X(i) = NaN;
    end
end
% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
% O-Mode
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range_O));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range_O));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length - O-Mode
Calc_Height_Group_O = sqrt((closest_group_range_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_O = sqrt((closest_phase_path_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_O = sqrt((closest_geometric_path_length_O.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - O-Mode
Calc_Height_Group_O = Calc_Height_Group_O - s_mat;
Calc_Height_Phase_O = Calc_Height_Phase_O - s_mat;
Calc_Height_Geo_O = Calc_Height_Geo_O - s_mat;

% Calculate the Triangular Height for Each Path Length - X-Mode
Calc_Height_Group_X = sqrt((closest_group_range_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_X = sqrt((closest_phase_path_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_X = sqrt((closest_geometric_path_length_X.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - X-Mode
Calc_Height_Group_X = Calc_Height_Group_X - s_mat;
Calc_Height_Phase_X = Calc_Height_Phase_X - s_mat;
Calc_Height_Geo_X = Calc_Height_Geo_X - s_mat;

% Plot
fig2 = figure;
hold on;

p2{1} = scatter(frequencies, Closest_Apogee_O, 'DisplayName', 'Ray Trace Apogee (O-Mode)');
p2{2} = scatter(frequencies, Closest_Apogee_X, 'DisplayName', 'Ray Trace Apogee (X-Mode)');
p2{3} = scatter(frequencies, real(Calc_Height_Geo_O), 'DisplayName', 'Calculated Height (Geometric Path, O-Mode)');
p2{4} = scatter(frequencies, real(Calc_Height_Geo_X), 'DisplayName', 'Calculated Height (Geometric Path, X-Mode)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);


% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time
if Site_Selector == 1
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nChesapeake to Auburn at %s', date_str);
else
    title_str = sprintf('Ionospheric Reflection Height vs Frequency\nCorpus Christi to Auburn at %s', date_str);
end
title(title_str);
ax2 = gca;
lgd2 = legend;
IEEE_Format_Plot_v2('figure', fig2, 'whr', 1, 'scale', 2, 'curve', p2, 'axis', ax2, 'legend', lgd2);

ThreeD_lats = [closest_lats_X closest_lats_O];
ThreeD_lons = [closest_lons_X closest_lons_O];

geo2 = figure;
geoscatter(ThreeD_lats,ThreeD_lons,[],"blue","x");


hold on;
geoscatter(32.604722, -85.486667,[],"red","o","filled");
hold off;
if Site_Selector == 1 %Chesapeake Title
    title(sprintf('Landing Point of Traced Rays from\nChesapeake, VA to Auburn, AL at %s', date_str));
else
    title(sprintf('Landing Point of Traced Rays from\nCorpus Christi, TX to Auburn, AL at %s', date_str));
end
legend('Ray Landing Points','Broun Hall','Location','best');
gx = gca; 
gx.LongitudeLabel.FontSize = 18; 
gx.LatitudeLabel.FontSize = 18;
gx.Title.FontSize = 18;
gx.FontName = 'Times New Roman';
gx.MapCenter = [32.604722, -85.486667];
gx.Basemap = 'topographic';
gx.GridLineStyle = '--';
gx.GridAlpha = 0.4;
gx.Legend.FontSize = 12;
gx.ZoomLevel = 13;
end    
if Simulation_Selector==3
% Analyze Anisotropy for Chesapeake, VA and Auburn, AL
% Constants
origin_ht = 0.0;

%% First Do Chesapeake to Auburn
% Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Chesapeake_Site,Auburn_Site);
% Generate the Ionospheric Parameters
[iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms,...
    geomag_grid_parms] = Ionospheric_Grid_3D(UT);
% Initialize Ionospheric Model
elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));
ray_bears = ones(size(elevs_init)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs_init, ray_bears,...
                        freqs_init, OX_mode, nhops, tol, iono_en_grid, ...
     	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
                        Bx, By, Bz, geomag_grid_parms);
%% O-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_O = zeros(size(frequencies));
Closest_Distance_O = zeros(size(frequencies));
closest_group_range_O = zeros(size(frequencies));
closest_phase_path_O = zeros(size(frequencies));
closest_geometric_path_length_O = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_O(i) = valid_apogees(idx);
        Closest_Distance_O(i) = valid_ranges(idx);
        closest_group_range_O(i) = valid_group_ranges(idx);
        closest_phase_path_O(i) = valid_phase_paths(idx);
        closest_geometric_path_length_O(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_O(i) = NaN;
        Closest_Distance_O(i) = NaN;
        closest_group_range_O(i) = NaN;
        closest_phase_path_O(i) = NaN;
        closest_geometric_path_length_O(i) = NaN;
    end
end
%% X-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_X = zeros(size(frequencies));
Closest_Distance_X = zeros(size(frequencies));
closest_group_range_X = zeros(size(frequencies));
closest_phase_path_X = zeros(size(frequencies));
closest_geometric_path_length_X = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = -1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_X(i) = valid_apogees(idx);
        Closest_Distance_X(i) = valid_ranges(idx);
        closest_group_range_X(i) = valid_group_ranges(idx);
        closest_phase_path_X(i) = valid_phase_paths(idx);
        closest_geometric_path_length_X(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_X(i) = NaN;
        Closest_Distance_X(i) = NaN;
        closest_group_range_X(i) = NaN;
        closest_phase_path_X(i) = NaN;
        closest_geometric_path_length_X(i) = NaN;
    end
end
% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
% O-Mode
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range_O));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range_O));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length - O-Mode
Calc_Height_Group_O_True = sqrt((closest_group_range_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_O_True = sqrt((closest_phase_path_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_O_True = sqrt((closest_geometric_path_length_O.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - O-Mode
Calc_Height_Group_O_True = Calc_Height_Group_O_True - s_mat;
Calc_Height_Phase_O_True = Calc_Height_Phase_O_True - s_mat;
Calc_Height_Geo_O_True = Calc_Height_Geo_O_True - s_mat;

% Calculate the Triangular Height for Each Path Length - X-Mode
Calc_Height_Group_X_True = sqrt((closest_group_range_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_X_True = sqrt((closest_phase_path_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_X_True = sqrt((closest_geometric_path_length_X.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - X-Mode
Calc_Height_Group_X_True = Calc_Height_Group_X_True - s_mat;
Calc_Height_Phase_X_True = Calc_Height_Phase_X_True - s_mat;
Calc_Height_Geo_X_True = Calc_Height_Geo_X_True - s_mat;

%% Then Do Auburn To Chesapeake
% Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Auburn_Site,Chesapeake_Site);
% Generate the Ionospheric Parameters
[iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms,...
    geomag_grid_parms] = Ionospheric_Grid_3D(UT);
% Initialize Ionospheric Model
elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));
ray_bears = ones(size(elevs_init)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs_init, ray_bears,...
                        freqs_init, OX_mode, nhops, tol, iono_en_grid, ...
     	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
                        Bx, By, Bz, geomag_grid_parms);
%% O-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_O = zeros(size(frequencies));
Closest_Distance_O = zeros(size(frequencies));
closest_group_range_O = zeros(size(frequencies));
closest_phase_path_O = zeros(size(frequencies));
closest_geometric_path_length_O = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_O(i) = valid_apogees(idx);
        Closest_Distance_O(i) = valid_ranges(idx);
        closest_group_range_O(i) = valid_group_ranges(idx);
        closest_phase_path_O(i) = valid_phase_paths(idx);
        closest_geometric_path_length_O(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_O(i) = NaN;
        Closest_Distance_O(i) = NaN;
        closest_group_range_O(i) = NaN;
        closest_phase_path_O(i) = NaN;
        closest_geometric_path_length_O(i) = NaN;
    end
end
%% X-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_X = zeros(size(frequencies));
Closest_Distance_X = zeros(size(frequencies));
closest_group_range_X = zeros(size(frequencies));
closest_phase_path_X = zeros(size(frequencies));
closest_geometric_path_length_X = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = -1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_X(i) = valid_apogees(idx);
        Closest_Distance_X(i) = valid_ranges(idx);
        closest_group_range_X(i) = valid_group_ranges(idx);
        closest_phase_path_X(i) = valid_phase_paths(idx);
        closest_geometric_path_length_X(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_X(i) = NaN;
        Closest_Distance_X(i) = NaN;
        closest_group_range_X(i) = NaN;
        closest_phase_path_X(i) = NaN;
        closest_geometric_path_length_X(i) = NaN;
    end
end
% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
% O-Mode
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range_O));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range_O));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length - O-Mode
Calc_Height_Group_O_Reverse = sqrt((closest_group_range_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_O_Reverse = sqrt((closest_phase_path_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_O_Reverse = sqrt((closest_geometric_path_length_O.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - O-Mode
Calc_Height_Group_O_Reverse = Calc_Height_Group_O_Reverse - s_mat;
Calc_Height_Phase_O_Reverse = Calc_Height_Phase_O_Reverse - s_mat;
Calc_Height_Geo_O_Reverse = Calc_Height_Geo_O_Reverse - s_mat;

% Calculate the Triangular Height for Each Path Length - X-Mode
Calc_Height_Group_X_Reverse = sqrt((closest_group_range_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_X_Reverse = sqrt((closest_phase_path_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_X_Reverse = sqrt((closest_geometric_path_length_X.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - X-Mode
Calc_Height_Group_X_Reverse = Calc_Height_Group_X_Reverse - s_mat;
Calc_Height_Phase_X_Reverse = Calc_Height_Phase_X_Reverse - s_mat;
Calc_Height_Geo_X_Reverse = Calc_Height_Geo_X_Reverse - s_mat;


% Plot
fig3 = figure;
hold on;

p3{1} = scatter(frequencies, Calc_Height_Geo_X_True, 'DisplayName', 'Calculated X-Mode Reflection Height (Chesapeake -> Auburn)');
p3{2} = scatter(frequencies, Calc_Height_Geo_X_Reverse, 'DisplayName', 'Calculated X-Mode Reflection Height (Auburn -> Chesapeake)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time

    title_str = sprintf('Anisotropy Analysis\nAuburn <-> Chesapeake at %s', date_str);

title(title_str);
ax3 = gca;
lgd3 = legend;
IEEE_Format_Plot_v2('figure', fig3, 'whr', 1, 'scale', 2, 'curve', p3, 'axis', ax3, 'legend', lgd3);

fig4 = figure;
hold on;

p4{1} = scatter(frequencies, Calc_Height_Geo_O_True, 'DisplayName', 'Calculated O-Mode Reflection Height (Chesapeake -> Auburn)');
p4{2} = scatter(frequencies, Calc_Height_Geo_O_Reverse, 'DisplayName', 'Calculated O-Mode Reflection Height (Auburn -> Chesapeake)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time

    title_str = sprintf('Anisotropy Analysis\nAuburn <-> Chesapeake at %s', date_str);

title(title_str);
ax4 = gca;
lgd4 = legend;
IEEE_Format_Plot_v2('figure', fig4, 'whr', 1, 'scale', 2, 'curve', p4, 'axis', ax4, 'legend', lgd4);
end

if Simulation_Selector==4
% Analyze Anisotropy for Corpus Christi, TX and Auburn, AL
% Constants
origin_ht = 0.0;

%% First Do Corpus Christi to Auburn
% Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(CorpusChristi_Site,Auburn_Site);
% Generate the Ionospheric Parameters
[iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms,...
    geomag_grid_parms] = Ionospheric_Grid_3D(UT);
% Initialize Ionospheric Model
elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));
ray_bears = ones(size(elevs_init)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs_init, ray_bears,...
                        freqs_init, OX_mode, nhops, tol, iono_en_grid, ...
     	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
                        Bx, By, Bz, geomag_grid_parms);
%% O-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_O = zeros(size(frequencies));
Closest_Distance_O = zeros(size(frequencies));
closest_group_range_O = zeros(size(frequencies));
closest_phase_path_O = zeros(size(frequencies));
closest_geometric_path_length_O = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_O(i) = valid_apogees(idx);
        Closest_Distance_O(i) = valid_ranges(idx);
        closest_group_range_O(i) = valid_group_ranges(idx);
        closest_phase_path_O(i) = valid_phase_paths(idx);
        closest_geometric_path_length_O(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_O(i) = NaN;
        Closest_Distance_O(i) = NaN;
        closest_group_range_O(i) = NaN;
        closest_phase_path_O(i) = NaN;
        closest_geometric_path_length_O(i) = NaN;
    end
end
%% X-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_X = zeros(size(frequencies));
Closest_Distance_X = zeros(size(frequencies));
closest_group_range_X = zeros(size(frequencies));
closest_phase_path_X = zeros(size(frequencies));
closest_geometric_path_length_X = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = -1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_X(i) = valid_apogees(idx);
        Closest_Distance_X(i) = valid_ranges(idx);
        closest_group_range_X(i) = valid_group_ranges(idx);
        closest_phase_path_X(i) = valid_phase_paths(idx);
        closest_geometric_path_length_X(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_X(i) = NaN;
        Closest_Distance_X(i) = NaN;
        closest_group_range_X(i) = NaN;
        closest_phase_path_X(i) = NaN;
        closest_geometric_path_length_X(i) = NaN;
    end
end
% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
% O-Mode
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range_O));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range_O));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length - O-Mode
Calc_Height_Group_O_True = sqrt((closest_group_range_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_O_True = sqrt((closest_phase_path_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_O_True = sqrt((closest_geometric_path_length_O.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - O-Mode
Calc_Height_Group_O_True = Calc_Height_Group_O_True - s_mat;
Calc_Height_Phase_O_True = Calc_Height_Phase_O_True - s_mat;
Calc_Height_Geo_O_True = Calc_Height_Geo_O_True - s_mat;

% Calculate the Triangular Height for Each Path Length - X-Mode
Calc_Height_Group_X_True = sqrt((closest_group_range_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_X_True = sqrt((closest_phase_path_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_X_True = sqrt((closest_geometric_path_length_X.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - X-Mode
Calc_Height_Group_X_True = Calc_Height_Group_X_True - s_mat;
Calc_Height_Phase_X_True = Calc_Height_Phase_X_True - s_mat;
Calc_Height_Geo_X_True = Calc_Height_Geo_X_True - s_mat;

%% Then Do Auburn To Corpus Christi
% Calculate Targeting Information
[azim, Target_Distance, origin_lat, origin_long] =...
    Bearing_Calculator(Auburn_Site,CorpusChristi_Site);
% Generate the Ionospheric Parameters
[iono_pf_grid,iono_pf_grid_5, collision_freq, Bx, By, Bz, iono_grid_parms,...
    geomag_grid_parms] = Ionospheric_Grid_3D(UT);
% Initialize Ionospheric Model
elevs_init = 20;
freq_init = 5;
freqs_init = freq_init.*ones(size(elevs_init));
ray_bears = ones(size(elevs_init)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

[ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs_init, ray_bears,...
                        freqs_init, OX_mode, nhops, tol, iono_en_grid, ...
     	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
                        Bx, By, Bz, geomag_grid_parms);
%% O-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_O = zeros(size(frequencies));
Closest_Distance_O = zeros(size(frequencies));
closest_group_range_O = zeros(size(frequencies));
closest_phase_path_O = zeros(size(frequencies));
closest_geometric_path_length_O = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = 1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_O(i) = valid_apogees(idx);
        Closest_Distance_O(i) = valid_ranges(idx);
        closest_group_range_O(i) = valid_group_ranges(idx);
        closest_phase_path_O(i) = valid_phase_paths(idx);
        closest_geometric_path_length_O(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_O(i) = NaN;
        Closest_Distance_O(i) = NaN;
        closest_group_range_O(i) = NaN;
        closest_phase_path_O(i) = NaN;
        closest_geometric_path_length_O(i) = NaN;
    end
end
%% X-Mode Virtual Heights
% Loop Using Already Initialized Ionospheric Model

elevs = 1:0.05:90;
frequencies = 1:0.1:20;  % 1 to 15 MHz

Closest_Apogee_X = zeros(size(frequencies));
Closest_Distance_X = zeros(size(frequencies));
closest_group_range_X = zeros(size(frequencies));
closest_phase_path_X = zeros(size(frequencies));
closest_geometric_path_length_X = zeros(size(frequencies));

ray_bears = ones(size(elevs)); % initial bearing of rays
ray_bears = ray_bears.*azim;
OX_mode = -1;
nhops = 1;
tol = [1e-7 0.01 25];  
max_distance_error = 100; % Max distance ray i allowed to miss target by (in km)

for i = 1:length(frequencies)
    freq = frequencies(i);
    freqs = freq .* ones(size(elevs));
    [ray_data, ray_path_data, ray_state_vec] = ...
            raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears,...
                        freqs, OX_mode, nhops, tol);
    
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
        Closest_Apogee_X(i) = valid_apogees(idx);
        Closest_Distance_X(i) = valid_ranges(idx);
        closest_group_range_X(i) = valid_group_ranges(idx);
        closest_phase_path_X(i) = valid_phase_paths(idx);
        closest_geometric_path_length_X(i) = valid_geometric_path_lengths(idx);
    else
        % No valid rays for this frequency
        Closest_Apogee_X(i) = NaN;
        Closest_Distance_X(i) = NaN;
        closest_group_range_X(i) = NaN;
        closest_phase_path_X(i) = NaN;
        closest_geometric_path_length_X(i) = NaN;
    end
end
% Calculate The Chord Length from Ground Distance (Straight-Line Distance)
% O-Mode
R = earthRadius('Kilometer');
Chord_Length = 2.*R.*sin(Target_Distance./(2.*R));
Chord_Length_Mat = ones(size(closest_group_range_O));
Chord_Length_Mat = Chord_Length.*Chord_Length_Mat;

% Calculate the Curvature Correction
s = R - sqrt(R.^2 - (Chord_Length.^2 ./ 4));
s_mat = ones(size(closest_group_range_O));
s_mat = s.*s_mat;

% Calculate the Triangular Height for Each Path Length - O-Mode
Calc_Height_Group_O_Reverse = sqrt((closest_group_range_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_O_Reverse = sqrt((closest_phase_path_O.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_O_Reverse = sqrt((closest_geometric_path_length_O.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - O-Mode
Calc_Height_Group_O_Reverse = Calc_Height_Group_O_Reverse - s_mat;
Calc_Height_Phase_O_Reverse = Calc_Height_Phase_O_Reverse - s_mat;
Calc_Height_Geo_O_Reverse = Calc_Height_Geo_O_Reverse - s_mat;

% Calculate the Triangular Height for Each Path Length - X-Mode
Calc_Height_Group_X_Reverse = sqrt((closest_group_range_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Phase_X_Reverse = sqrt((closest_phase_path_X.^2 - Chord_Length_Mat.^2)./4);
Calc_Height_Geo_X_Reverse = sqrt((closest_geometric_path_length_X.^2 - Chord_Length_Mat.^2)./4);

% Correct for the "Thickness" of Earth at the Center of the Chord - X-Mode
Calc_Height_Group_X_Reverse = Calc_Height_Group_X_Reverse - s_mat;
Calc_Height_Phase_X_Reverse = Calc_Height_Phase_X_Reverse - s_mat;
Calc_Height_Geo_X_Reverse = Calc_Height_Geo_X_Reverse - s_mat;

% Plot
fig5 = figure;
hold on;

p5{1} = scatter(frequencies, Calc_Height_Geo_X_True, 'DisplayName', 'Calculated X-Mode Reflection Height (Corpus Christi -> Auburn)');
p5{2} = scatter(frequencies, Calc_Height_Geo_X_Reverse, 'DisplayName', 'Calculated X-Mode Reflection Height (Auburn -> Corpus Christi)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time

    title_str = sprintf('Anisotropy Analysis\nAuburn <-> Corpus Chrsti at %s', date_str);

title(title_str);
ax5 = gca;
lgd5 = legend;
IEEE_Format_Plot_v2('figure', fig5, 'whr', 1, 'scale', 2, 'curve', p5, 'axis', ax5, 'legend', lgd5);

fig6 = figure;
hold on;

p6{1} = scatter(frequencies, Calc_Height_Geo_O_True, 'DisplayName', 'Calculated O-Mode Reflection Height (Corpus Christi -> Auburn)');
p6{2} = scatter(frequencies, Calc_Height_Geo_O_Reverse, 'DisplayName', 'Calculated O-Mode Reflection Height (Auburn -> Corpus Christi)');

hold off;

xlabel('Frequency (MHz)');
ylabel('Height (km)');
legend('Location', 'best');
grid on;
xlim([0 max(frequencies)]);
ylim([0 400]);

% Format date and time string
date_str = sprintf('%04d-%02d-%02d %02d:%02d UTC', UT(1), UT(2), UT(3), UT(4), UT(5));

% Set title with site info and date/time

    title_str = sprintf('Anisotropy Analysis\nAuburn <-> Corpus Christi at %s', date_str);

title(title_str);
ax6 = gca;
lgd6 = legend;
IEEE_Format_Plot_v2('figure', fig6, 'whr', 1, 'scale', 2, 'curve', p6, 'axis', ax6, 'legend', lgd6);
end

toc