% DESCRIPTION: 
%   This script calculates Ensemble-averaged Mean Square Displacement(MSD) 
%   and fit with Anomalous Diffusion model from the x- and y- coordinates 
%   of the trajectories.
%
% REQUIREMENTS:
%   Curve Fitting Toolbox
%
% INPUT:
%   fileName : Input CSV file name without extension. You should save the 
%              the trajectory CSV file in the same folder as this script.
%              The CSV file can be converted from NIS Elements (Nikon)
%              format using Script_ConvertTrajectoryData.m. 
%
% OUTPUT:
%   This script export following three result csv files.
%   [fileName]_MSD.csv : 
%       A result CSV file in the following format.
%       ----------------------------------
%       Time   Ave  SD   SE   1    2    3         
%       0.00   0.0  0.0  1.3  2.3  1.3  2.3
%       0.03   0.2  0.2  0.3  3.1  0.3  3.1
%       0.06   0.4  0.3  1.2  1.2  0.3  3.1
%       -----------------------------------
%   [fileName]_FitParam.csv : 
%       Table of fitted parameter.
%   [fileName]_FitModel.csv : 
%       Fitted model curve with the time interval of 1/precision.
%
% CHANGELOG:
%   1.0 (2021-09-13) Released by Yuma Ito <yitou@bio.titech.ac.jp>
%   1.1 (2021-11-09) Updated comments by Yuma Ito <yitou@bio.titech.ac.jp>


% Start script
fileName = 'eu2-1733_frm1-30';
CalcFitAnnomMSD(fileName);


% ----- Functions -----
function CalcFitAnnomMSD(fileName)
    % Execute all functions
    [timeList,trjList] = GetData(fileName);
    CheckTime(timeList)
    msdList = CalcMSD(trjList);
    [msdAve, msdSD, msdSEM, msdN] = AverageMSD(msdList);
    SaveMSD(timeList,msdList,msdAve,msdSD,msdSEM,msdN,fileName);
    
    [result, goodness] = FitMSD(timeList,msdAve);
    SaveFitResult(result,fileName);
    
    precision = 10;
    SaveFitModel(timeList,result,precision,fileName)
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

function msdList = CalcMSD(trjList)
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        x = trj(:,1)'; y = trj(:,2)';
        msd = [];
        for i = 1:length(trj)
            msd = [msd mean((x(i:end)-x(1:end-i+1)).^2 + (y(i:end)-y(1:end-i+1)).^2)];
        end
        msdList{trjNo} = msd;
    end
end

function [msdAve, msdSD, msdSEM, msdN] = AverageMSD(msdList)
    %@Summarize statistical values
    lenList = cellfun(@length,msdList);
    maxLen = max(lenList);
    msdAve = zeros(1,maxLen); msdN = zeros(1,maxLen); 
    msdSD = zeros(1,maxLen); msdSEM = zeros(1,maxLen);
    for i = 1:maxLen
        toSel = lenList >= i;
        msdSel = msdList(toSel);
        msd = cellfun(@(x) x(i),msdSel);
        msdAve(i) = mean(msd);
        msdN(i) = length(msd);
        msdSD(i) = std(msd);
        msdSEM(i) = std(msd)/sqrt(length(msd));
    end
    
end

function SaveMSD(timeList,msdList,msdAve,msdSD,msdSEM,msdN,fileName)
    filePath = [pwd filesep fileName '_MSD.csv'];
    
    % Summarize statistical values
    dataList = [timeList' msdAve' msdSD' msdSEM' msdN'];
    
    % Add individual MSD
    totalRow = size(dataList,1);
    nameList = {};
    for trjNo = 1:length(msdList)
        msd = msdList{trjNo};
        if length(msd) < totalRow
            msd = [msd nan(1,totalRow-length(msd))];
        end
        dataList = [dataList msd'];
        nameList{trjNo} = ['trj' num2str(trjNo)];
    end
    
    % Save as table
    t = array2table(dataList);
    Header = cat(2,{'TimeInteval','EnsembleAveragedMSD','SD','SEM','Count'},nameList);
    t.Properties.VariableNames = Header;
    writetable(t,filePath)
    disp('MSD exported')
end

function [Result,Goodness] = FitMSD(timeList,msdAve)
    % Fit MSD with Anomalous diffusion model using curve fitting toolbox
    
    x = timeList';
    y = msdAve';
    
    model = '4*D*x^a';
    independent = 'x';
    coefficients = {'D','a'};
    
    Type = fittype(model,'dependent','y','independent',independent,'coefficient',coefficients); 
    Options = fitoptions(Type);

    dt = x(2)-x(1); 
    Options.StartPoint(1) = y(2)/(4*dt);
    Options.Lower(1) = 0;
    Options.Upper(1) = 100;

    Options.StartPoint(2) = 1;
    Options.Lower(2) = 0;
    Options.Upper(2) = 10;
        
    [Result,Goodness] = fit(x,y,Type,Options);
end

function SaveFitResult(result,fileName)
    filePath = [pwd filesep fileName '_MSD_FitParam.csv'];
    RowName = {'FitValue','ConfidenceInterval95%High','ConfidenceInterval95%Low'}';
    t = [cell2table(RowName) array2table([coeffvalues(result);confint(result)])];
    t.Properties.VariableNames = {'RowName','DiffusionCoefficient','AnomalousExponent'};
    writetable(t,filePath)
    disp('MSD FitResult exported')
end

function SaveFitModel(timeList,result,precision,fileName)
    filePath = [pwd filesep fileName '_MSD_FitModel.csv'];
    d = timeList(2)-timeList(1);
    x = timeList(1):d/precision:timeList(end);
    y = result(x);
    t = array2table([x' y]);
    t.Properties.VariableNames = {'TimeInterval','FitModel'};
    writetable(t,filePath)
    disp('MSD FitModel exported')
end
