close all;
clear all;
clc;

f = 1000:1000:5000;
P = 1./f


% determine dimensions
dt = min(P) / 360
dur = max(P)

t = single(0:dt:dur);

A = zeros(360,360,360,360,length(t));

