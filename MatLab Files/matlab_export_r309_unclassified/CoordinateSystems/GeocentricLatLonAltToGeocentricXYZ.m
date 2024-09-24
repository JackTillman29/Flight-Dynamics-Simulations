function [r] = GeocentricLatLonAltToGeocentricXYZ(lambda,phigc,range)
    r = [ ...
        range .* cos(phigc) .* cos(lambda);
        range .* cos(phigc) .* sin(lambda);
        range .* sin(phigc) ...
        ];
end