% DESCRIPTION: 
%   This script converts an XLSX file exported from NISElements (Nikon) to 
%   a CSV file that can be imported by following scripts.
%    * Script_CalcConfASpringK.m
%    * CalcFitAnomalousMSD.m
%
% INPUT:
%   fileName : Input XLSX file name without extension. You should save the 
%              XLSX file in the same folder as this script.
%   pitch : Pixel size of image [um/pix].
%   interval : The frame interval time of the original movie [s].
%
% OUTPUT:
%   A trajectory CSV file in the following format.
%   ----------------------------------
%      TrjNo    1    1    2    2  ...
%      x_or_y   x    y    x    y
%      0.00     0.0  0.0  1.3  2.3
%      0.03     0.2  0.2  0.3  3.1
%      0.06     0.4  0.3  1.2  1.2
%    -----------------------------------
%
% CHANGELOG:
%   1.0 (2021-09-13) Released by Yuma Ito <yitou@bio.titech.ac.jp>
%   1.1 (2021-11-09) Updated comments by Yuma Ito <yitou@bio.titech.ac.jp>


% Parameters 
pitch = 0.0338; % [um/pix]
interval = 0.551; % [s]
fileName = 'eu2-1733_561';

% Start script
filePath = [pwd filesep fileName '.xlsx'];
t = readtable(filePath);
trjNoList = unique(t.ID);

% Read xlsx table
x = {}; y = {}; timeList = {};
for i = 1:length(trjNoList)
    trjNo = trjNoList(i);
    tTrj = t(t.ID == trjNo,:);
    xPix = tTrj.PositionX_px_;
    yPix = tTrj.PositionY_px_;
    stepNo = 0:length(xPix)-1;
    
    x{i} = xPix * pitch;
    y{i} = yPix * pitch;
    timeList{i} = stepNo * interval;
end

% Convert to cell 
[maxLen,posMax] = max(cellfun(@length,x));
Cell = cat(2,{'trjNo','x_or_y'},num2cell(timeList{posMax}))';
for i = 1:length(trjNoList)
    xCell = cat(1,num2cell(trjNoList(i)),{'x'});
    yCell = cat(1,num2cell(trjNoList(i)),{'y'});
    xCell = cat(1,xCell,num2cell(x{i}));
    yCell = cat(1,yCell,num2cell(y{i}));
    if size(xCell,1) < maxLen + 2
        padLen = maxLen + 2 - size(xCell,1);
        xCell = cat(1,xCell,num2cell(nan(padLen,1)));
        yCell = cat(1,yCell,num2cell(nan(padLen,1)));
    end

    Cell = cat(2,Cell,xCell);
    Cell = cat(2,Cell,yCell);
end

% Expor csv
filePath = [pwd filesep fileName '.csv'];
writecell(Cell,filePath)
disp('CSV file exported.')