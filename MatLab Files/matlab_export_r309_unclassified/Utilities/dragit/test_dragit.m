% test "dragit"

close all; clear all; clc

Fs = 10;
dt = 1/Fs;
t = [0:dt:10];

x = randn(1,length(t));

figure
plot(t,x,'.-','MarkerSize',10)

% simply call "dragit" after plotting something
%  (dragit will look for the last line added to the current axes object)
dragit







