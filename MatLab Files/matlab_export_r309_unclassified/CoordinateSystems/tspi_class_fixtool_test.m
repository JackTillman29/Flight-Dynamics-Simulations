close all; clear all; clc

% add directory where "tspi_class.m" is located to the MATLAB search path
addpath('H:\MATLAB\CoordinateSystems\')

tp = tspi_class('unclass_sample_tspi.csv');

% update any search strings to match column headers of the specific file
tp.searchString_wgs84_alt = 'GEOHT84';
tp.searchString_roll      = 'ROLL';
tp.searchString_pitch     = 'PITCH';
tp.searchString_yaw       = 'YAW';

tp.fix_TSPI('platform',1,'run',1);
