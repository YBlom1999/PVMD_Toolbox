function [DeviceTable, DeviceInfo] = read_device(DeviceName)
% Reads DeviceTable from xlsx file. All columns are read as string.
% Note: conversion from string to categorical (Layer, Unit) or double
% (Thickness) is NOT done here, but in DevicePlot2D and GENPRO_shell.

% load device structure from xlsx file
DeviceFilePath = [DeviceName,'.xlsx'];

%---options for readtable function---
opts = spreadsheetImportOptions("NumVariables", 3);
opts.VariableNames = ["Layer", "Thickness", "Unit"];
opts.VariableTypes = ["string", "string", "string"];
%---

DeviceTable = readtable(DeviceFilePath,opts,'Sheet','Layers'); % READ FILE

%---options for readtable function---
opts = spreadsheetImportOptions("NumVariables", 2);
opts.VariableTypes = ["string", "string"];
DeviceInfo = readtable(DeviceFilePath,opts,'Sheet','Info'); % READ FILE

DeviceTable(1,:) = [];      % ignore header (line 1)

end
 % TODO: test somewhere if it is a propper table (3 columns) and holds 
 % valid data.