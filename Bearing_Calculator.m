function [TargetHeading,TargetDistance] = Bearing_Calculator(LatDegSrc,...
    LatMinSrc,LatSecSrc,LatPolSrc,LongDegSrc,LongMinSrc,LongSecSrc,...
    LongPolSrc,LatDegDest,LatMinDest,LatSecDest,LatPolDest,LongDegDest,...
    LongMinDest,LongSecDest,LongPolDest)
%BEARING_CALCULATOR Calculates True North Heading from Source to
%Destination 
%   Uses Degree, Minute, Second (DMS) valued latitude and longitude for any two
%   points on the globe to calculate a heading at the source pointing to
%   the destination in relation to true north. Function will have
%   unexpected/suboptimal behavior in polar regions or around the
%   international date line. Distance accuracy will decrease as points move
%   further from the equator.
arguments (Input)
    LatDegSrc %The Degree value of the source latitude
    LatMinSrc %The Arcminute value of the source latitude
    LatSecSrc %The Arcsecond value of the source latitude
    LatPolSrc %The Polarity of source latitude (1=North, -1=South)
    LongDegSrc %The Degree value of the source longitude
    LongMinSrc %The Arcminute value of the source longitude
    LongSecSrc %The Arcsecond value of the source longitude
    LongPolSrc %The Polarity of source longitude (1=East, -1=West)
    LatDegDest %The Degree value of the destination latitude
    LatMinDest %The Arcminute value of the destination latitude
    LatSecDest %The Arcsecond value of the destination latitude
    LatPolDest %The Polarity of destination latitude (1=North, -1=South)
    LongDegDest %The Degree value of the destination longitude
    LongMinDest %The Arcminute value of the destination longitude
    LongSecDest %The Arcsecond value of the destination longitude
    LongPolDest %The Polarity of destination longitude (1=East, -1=West)
end

arguments (Output)
    TargetHeading
    TargetDistance
end

%Converts DMS to signed Decimal Coordinates for Source
LatSrcDEC = (LatDegSrc + (LatMinSrc + (LatSecSrc/60))/60);
LongSrcDEC = (LongDegSrc + (LongMinSrc + (LongSecSrc/60))/60);
LatSrcDEC = LatSrcDEC*LatPolSrc;
LongSrcDEC = LongSrcDEC*LongPolSrc;

%Converts DMS to signed Decimal Coordinates for Destination
LatDestDEC = (LatDegDest + (LatMinDest + (LatSecDest/60))/60);
LongDestDEC = (LongDegDest + (LongMinDest + (LongSecDest/60))/60);
LatDestDEC = LatDestDEC*LatPolDest;
LongDestDEC = LongDestDEC*LongPolDest;

lat1 = deg2rad(LatSrcDEC);
lon1 = deg2rad(LongSrcDEC);
lat2 = deg2rad(LatDestDEC);
lon2 = deg2rad(LongDestDEC);
dLat = lat2 - lat1;
dLon = lon2 - lon1;

% Calculate bearing using forward azimuth formula
x = sin(dLon) * cos(lat2);
y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
TargetHeading = rad2deg(atan2(x, y));
TargetHeading = mod(TargetHeading, 360);

% Calculate distance using Haversine formula
R = 6371;  % Earth's mean radius in kilometers
a = sin(dLat/2)^2 + cos(lat1) * cos(lat2) * sin(dLon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
TargetDistance = R * c;
end