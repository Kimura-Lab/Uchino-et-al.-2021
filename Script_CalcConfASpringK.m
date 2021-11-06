% 軌跡のxy座標のリストから Area of Confinement と Spring Constant を計算する。
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Input形式 
%    このスクリプト本体[Calc_Fit_Anomalous_MSD.m]と同じフォルダ内に、
%    FileName.csvという以下の形式のcsvファイルを保存。
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
% ※ 短い軌跡のx,y座標は空欄でもかまわない。
% ※ 単位は時間 s, 距離 um にそろえてからcsvにしてください。軌跡の開始を 0 s とします。
% ※ 各軌跡の始点は (0,0) でなくて良い。プログラムが(0,0)に平行移動する。
 
% ----- Output のファイル形式
%   FileName_ConfA_SpringK.csv
%   TrjNo  Rlong Rshort ConfA SpringK         
%   1      0.0   0.0    1.3   2.3
%   2      0.2   0.2    0.3   3.1
%   3      0.4   0.3    1.2   1.2
%   
%   TrjNo: 軌跡番号の連番。Input csvファイルの番号には対応せず、単に順番
%   Rlong: Trajectory Point の 95%が含まれる楕円の長半径 (単位は um)
%   Rshort: 楕円の短半径 (単位は um)
%   ConfA: Trajectory Point の 95%が含まれる楕円の面積 (単位は um^2)
%   SpringK: ばね定数。(単位は kbT/um^2)

% Reference
% Confinment area :
%   Germier, T. et al. Biophys J 113, 1383?1394 (2017).
%   http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix
%
% Spring coefficient : 
%   Shukron, O., Seeber, A., Amitai, A. & Holcman, D.
%   Trends Genet 35, 685?705 (2019).

% <<ファイル名を入力して Cntrol + Enter を押して実行>>
fileName = 'eu2-1733_frm1-30'
CalcConfASpringK(fileName);


% ----- 以下は関数 -----
function CalcConfASpringK(fileName)
    % 全ての処理を実行する関数
    [timeList,trjList] = GetData(fileName);
    CheckTime(timeList)
    [rlList, rsList, confAList] = CalcConfA(trjList);
    kList = CalcSpringK(timeList, trjList);
    SaveConfASpringK(rlList,rsList,confAList,kList,fileName);
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
        disp('時間間隔が均一ではありません。正しいSpringKが計算できません。')
    end
end

function [rlList, rsList, confAList] = CalcConfA(trjList)
    chisqr = 2.4477;
    rlList = {}; rsList = {}; confAList = {};
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        X = [trj(:,1) trj(:,2)];
        [evec, eval] = eig(cov(X));
        [evecMax_col,~] = find(eval == max(max(eval)));
        evecMax = evec(:, evecMax_col);
        evalMax = max(max(eval));
        if(evecMax_col == 1)
            evalMin = max(eval(:,2));
            evecMin = evec(:,2);
        else
            evalMin = max(eval(:,1));
            evecMin = evec(1,:);
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
    %   Shukron et al Trends Genet 2019 の Box 2 参照
    % - D * k * x + b となるので 傾き a = -D*k よって k = -a/D
    
    kList = {};
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        xTrj = trj(:,1); yTrj = trj(:,2);
        dx = diff(xTrj); dy = diff(yTrj);
        xc = mean(xTrj); yc = mean(yTrj);
        
        % x = Distance from centroid (1-dimentional)
        % Displaceと対応するので最後の１点を除く
        x = [(xTrj(1:end-1)-xc); (yTrj(1:end-1)-yc)];
        
        % y = velocity = displacement/interavl (1-dimentional)
        dt = mean(diff(timeList));
        y = [dx; dy]/dt;
        
        % 線形回帰
        X = [ones(length(x),1) x];
        b = X\y;
        
        % k = -傾き/D, (1step displacement^2 = 4*D*t)
        disp = sqrt(dx.^2 + dy.^2);
        D = mean(disp.^2./(dt * 4));
        kList{trjNo} = -b(2)/D;
    end
    kList = cell2mat(kList)';
end

function SaveConfASpringK(rlList,rsList,confAList,kList,fileName)
    filePath = [pwd filesep fileName '_ConfA_SpringK.csv'];
    
    trjNo = 1:length(rlList);
    dataList = [trjNo' rlList rsList confAList kList];
    
    % Tableにして出力する
    t = array2table(dataList);
    Header = {'TrjNo','Rlong','Rshort','ConfA','SpringK'};
    t.Properties.VariableNames = Header;
    writetable(t,filePath)
    disp('ConfA SpringK exported')
end
