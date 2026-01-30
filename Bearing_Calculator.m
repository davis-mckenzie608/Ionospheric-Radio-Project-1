function [TargetHeading,TargetDistance] = Bearing_Calculator(LatDegSrc,LatMinSrc,...
    LatSecSrc,LatPolSrc,LongDegSrc,LongMinSrc,LongSecSrc,LongPolSrc,...
    LatDegDest,LatMinDest,LatSecDest,LatPolDest,LongDegDest,LongMinDest,...
    LongSecDest,LongPolDest)
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
LongSrcDEC = (LongDegSrc + (LongMinSrc + (LongSecSrc/60)/60));
LatSrcDEC = LatSrcDEC*LatPolSrc;
LongSrcDEC = LongSrcDEC*LongPolSrc;
Src = [LatSrcDEC,LongSrcDEC];

%Converts DMS to signed Decimal Coordinates for Destination
LatDestDEC = (LatDegDest + (LatMinDest + (LatSecDest/60))/60);
LongDestDEC = (LongDegDest + (LongMinDest + (LongSecDest/60)/60));
LatDestDEC = LatDestDEC*LatPolDest;
LongDestDEC = LongDestDEC*LongPolDest;
Dest = [LatDestDEC,LongDestDEC];

ElevDecl = Dest(1) - Src(1);
Azm = Dest(2) - Src(2);

TargetHeading = atan(ElevDecl/Azm);
TargetDistanceDeg = sqrt(ElevDecl^2 + Azm^2);
TargetDistance = TargetDistanceDeg * 111; % Convert degrees to km (approx.)
end