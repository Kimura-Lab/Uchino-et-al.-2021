% 軌跡のxy座標のリストから平均MSDを算出してAnomalous diffusion model で fittingする。
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Curve Fitting Toolbox が必要。[アプリ] > [さらにアプリを取得] からインストール。
% 2) Input形式 このスクリプト本体[Script_CalcFitAnomalousMSD.m]と同じフォルダ内に、
%    [FileName].csvという以下の形式のcsvファイルを保存。
%    ファイル名は任意で、以下の fileName パラメータの値に入力する。
%
% ----- Input csv のファイル形式 -----
%   TrjNo    1    1    2    2  ...
%   x_or_y   x    y    x    y
%   0.00     0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2
% -----------------------------------
% ※ 1列目の3行目からは、ステップ時間を入力してください。
% ※ 2列目以降は、軌跡ごとに、その時間におけるx,y座標を入力してください。
% ※ 短い軌跡のx,y座標は空欄でもかまわないが、平均MSDのfittingは最長で計算される。
% 　　Fittingする長さを短くしたい場合は、Fittingしたくない時間以降の行を削除したcsvを用意してください。
% ※ 単位は時間 s, 距離 um にそろえてからcsvにしてください。軌跡の開始を 0 s とします。
% ※ 各軌跡の始点は (0,0) でなくて良い。プログラムが(0,0)に平行移動する。
 
% ----- Output のファイル形式
%   FileName_MSD.csv
%   Time   Ave  SD   SE   1    2    3         
%   0.00   0.0  0.0  1.3  2.3  1.3  2.3
%   0.03   0.2  0.2  0.3  3.1  0.3  3.1
%   0.06   0.4  0.3  1.2  1.2  0.3  3.1

%   FileName_FitParam.csv
%   
%
%
%   FileName_FitModel.csv
%   時間間隔を 1/precision にした Fitモデルの曲線。

% <<ファイル名を入力して Cntrol + Enter を押して実行>>
fileName = 'eu2-1733_frm1-20'
CalcFitAnnomMSD(fileName);


% ----- 以下は関数 -----
function CalcFitAnnomMSD(fileName)
    % 全ての処理を実行する関数
    
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
    
    % Csvファイルを読込み、Cell配列にする
    filePath = [pwd filesep fileName '.csv']
    fid = fopen(filePath);
    rows = textscan(fid,'%s');
    fclose(fid);
    rowCell = cellfun(@(x) split(x,','), rows,'UniformOutput',false);
    rowCell = rowCell{1};
    
    % 時間行を取り出す
    timeCell = rowCell(3:end,1)';
    timeList = cellfun(@str2double,timeCell);
    
    % 各軌跡のxy座標を取り出す
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
        
        %欠損値除去
        x = rmmissing(x);
        y = rmmissing(y);
        
        % 軌跡の最初の座標を(0,0)にする
        x = x - x(1);
        y = y - y(1); 
        trjList{trjNo} = [x y];
    end
    
end

function CheckTime(timeList)
    % 時間が等間隔か確認する
    %  浮動小数点誤差が10^-15程度あるので、差の差(=0)を10^10倍して四捨五入し、
    %  0になるかどうかで判断する
    if length(unique(round(diff(diff(timeList))*10^10))) > 1
        disp('時間間隔が均一ではありません。正しいMSDが計算できません。')
    end
end

function msdList = CalcMSD(trjList)
    % 超賢いMSD計算法
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
    %　統計量の計算
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
    
    % 統計量をまとめる
    dataList = [timeList' msdAve' msdSD' msdSEM' msdN'];
    
    % 個々のMSDを追加する。
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
    
    % Tableにして出力する
    t = array2table(dataList);
    Header = cat(2,{'TimeInteval','EnsembleAveragedMSD','SD','SEM','Count'},nameList);
    t.Properties.VariableNames = Header;
    writetable(t,filePath)
    disp('MSD exported')
end

function [Result,Goodness] = FitMSD(timeList,msdAve)
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
