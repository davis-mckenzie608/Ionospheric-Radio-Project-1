%ELEC 6970 - Ionospheric Radio: Project 1
%JDM 0206 - Spring 2026
%Ray Tracing

CorpusChristi_Site = [27 47 45 1 97 24 14 -1];
Chesapeake_Site = [36 46 06 1 76 17 13 -1];
Auburn_Site = [32 36 17 1 85 29 12 -1];
UT = [2004; 7; 1; 10; 0]; %Time as Year, Month, Day, Hour, Minute

[iono_en_grid, iono_en_grid_5, collison_freq, irreg, iono_te_grid] =...
    Ionospheric_Grid_2D(CorpusChristi_Site,Auburn_Site,UT);