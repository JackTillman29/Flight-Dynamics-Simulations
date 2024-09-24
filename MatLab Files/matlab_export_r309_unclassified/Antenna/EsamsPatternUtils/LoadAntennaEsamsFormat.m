function ant = LoadAntennaEsamsFormat(antfile)

antArray = AntennaPatternReader(antfile);
%convert esams table to matlab arrays usable for interp2
[antAzVec,antElVec,antPat] = esamsToMatlabTable(antArray);
[antAz,antEl] = meshgrid(antAzVec,antElVec);
antPat = antPat.';

ant.AzVec = antAzVec;
ant.ElVec = antElVec;
ant.Az    = antAz;
ant.El    = antEl;
ant.Pat   = antPat;

end