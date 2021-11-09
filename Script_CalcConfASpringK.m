% DESCRIPTION: 
%   This script calculates Area of Confinement and Effective Spring
%   Coefficient from the x- and y- coordinates of the trajectories.
%
% INPUT:
%   fileName : Input CSV file name without extension. You should save the 
%              the trajectory CSV file in the same folder as this script.
%              The CSV file can be converted from NISElements (Nikon)
%              format using Script_ConvertTrajectoryData.m. 
%
% OUTPUT:
%   A result CSV file in the following format inserting the footer
%   '_ConfA_SpringK'.
%   ----------------------------------
%   TrjNo  Rlong Rshort ConfA SpringK         
%   1      0.0   0.0    1.3   2.3
%   2      0.2   0.2    0.3   3.1
%   3      0.4   0.3    1.2   1.2
%   -----------------------------------
%   TrjNo: Trajectory number that is independent of the input CSV file.
%   Rlong: Semimajor axis of an ellipse containing 95% of Trajectory 
%          Points (in um)
%   Rshort: Semiminor axis of an ellipse containing 95% of Trajectory 
%           Points (in um)
%   ConfA: Area of an ellipse containing 95% of Trajectory Points 
%          (in um^2)
%   SpringK: Effective Spring Coefficient (in kBT/um^2)
%
% REFERENCE: 
%   1] Germier, T. et al. Biophys J 113, 1383-1394 (2017).
%      http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-c
%      ovariance-matrix
%   2] Shukron, O., Seeber, A., Amitai, A. & Holcman, D. Trends Genet 35, 
%      685-705 (2019).
%
% CHANGELOG:
%   1.0 (2021-09-13) Released by Yuma Ito <yitou@bio.titech.ac.jp>
%   1.1 (2021-11-09) Updated comments by Yuma Ito <yitou@bio.titech.ac.jp>


% Start script
fileName = 'eu2-1733_frm1-30';
CalcConfASpringK(fileName);

% ----- Functions -----
function CalcConfASpringK(fileName)
    % Execute all functions
    [timeList,trjList] = GetData(fileName);
    CheckTime(timeList)
    [rlList, rsList, confAList] = CalcConfA(trjList);
    kList = CalcSpringK(timeList, trjList);
    SaveConfASpringK(rlList,rsList,confAList,kList,fileName);
end

function [timeList,trjList] = GetData(fileName)
    % Read CSV file and create cell
    filePath = [pwd filesep fileName '.csv']
    fid = fopen(filePath);
    rows = textscan(fid,'%s');
    fclose(fid);
    rowCell = cellfun(@(x) split(x,','), rows,'UniformOutput',false);
    rowCell = rowCell{1};
    
    % Get time row
    timeCell = rowCell(3:end,1)';
    timeList = cellfun(@str2double,timeCell);
    
    % Get xy-coordinates
    totalTrj = (size(rowCell,2)-1)/2;
    trjList = {};
    for trjNo = 1:totalTrj
        xCell = rowCell(3:end,2*trjNo-1+1);
        yCell = rowCell(3:end,2*trjNo  +1);
        toValue = ~cellfun(@isempty,xCell);
        xCell = xCell(toValue);
        yCell = yCell(toValue);
        x = cellfun(@str2num,xCell);
        y = cellfun(@str2num,yCell);
        
        % Remove NaN
        x = rmmissing(x);
        y = rmmissing(y);
        
        % Set first position as (0,0)
        x = x - x(1);
        y = y - y(1); 
        trjList{trjNo} = [x y];
    end
end

function CheckTime(timeList)
    % Check time interval integrity
    if length(unique(round(diff(diff(timeList))*10^10))) > 1
        disp('Time interval is invalid. Spring K can not be calculated.')
    end
end

function [rlList, rsList, confAList] = CalcConfA(trjList)
    chisqr = 2.4477;
    rlList = {}; rsList = {}; confAList = {};
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        X = [trj(:,1) trj(:,2)];
        [~, eval] = eig(cov(X));
        [evecMax_col,~] = find(eval == max(max(eval)));
        evalMax = max(max(eval));
        if(evecMax_col == 1)
            evalMin = max(eval(:,2));
        else
            evalMin = max(eval(:,1));
        end
        a = chisqr * sqrt(evalMax);
        b = chisqr * sqrt(evalMin);
        area = pi * a * b;
        
        rlList{trjNo} = a;
        rsList{trjNo} = b;
        confAList{trjNo} = area;
    end
    rlList = cell2mat(rlList)';
    rsList = cell2mat(rsList)';
    confAList = cell2mat(confAList)';
end

function kList = CalcSpringK(timeList, trjList)
    %   Ref: Shukron et al Trends Genet 2019 Box 2
 
    kList = {};
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        xTrj = trj(:,1); yTrj = trj(:,2);
        dx = diff(xTrj); dy = diff(yTrj);
        xc = mean(xTrj); yc = mean(yTrj);
        
        % x = Distance from centroid (1-dimentional)
        x = [(xTrj(1:end-1)-xc); (yTrj(1:end-1)-yc)];
        
        % y = velocity = displacement/interavl (1-dimentional)
        dt = mean(diff(timeList));
        y = [dx; dy]/dt;
        
        X = [ones(length(x),1) x];
        b = X\y;
        
        % k = -slope/D, (1step displacement^2 = 4*D*t)
        disp = sqrt(dx.^2 + dy.^2);
        D = mean(disp.^2./(dt * 4));
        kList{trjNo} = -b(2)/D;
    end
    kList = cell2mat(kList)';
end

function SaveConfASpringK(rlList,rsList,confAList,kList,fileName)
    % Export as Table
    filePath = [pwd filesep fileName '_ConfA_SpringK.csv'];
    
    trjNo = 1:length(rlList);
    dataList = [trjNo' rlList rsList confAList kList];

    t = array2table(dataList);
    Header = {'TrjNo','Rlong','Rshort','ConfA','SpringK'};
    t.Properties.VariableNames = Header;
    writetable(t,filePath)
    disp('ConfA SpringK exported')
end
