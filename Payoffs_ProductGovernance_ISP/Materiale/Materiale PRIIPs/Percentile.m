function [PrcValue,Position] = Percentile(Data,Perc,RefCol)

% Manual percentile calculation
    NumSim = size(Data,1);
    Pointer = (1:1:NumSim)';
    Data = [Pointer,Data];
    DataSort = sortrows(Data,RefCol+1);
    Perc = floor(Perc/100*NumSim);
    Position = DataSort(Perc,1);
    PrcValue = DataSort(Perc,2:end);