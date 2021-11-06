% 軌跡解析データ XLSX ファイルから 解析用のcsv形式に変更する。
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Input形式 このスクリプト本体[ConvertTrajectoryData.m]と同じフォルダ内に、
%    [FileName].xlsx という以下の形式のcsvファイルを保存。
%    ファイル名は任意で、以下の fileName パラメータの値に入力する。
%
% 2) 以下のパラメーターを入力する。
%     pitch [um/pix]
%     interval [s]
%
% ----- Output のファイル形式
%   TrjNo    1    1    2    2  ...
%   x_or_y   x    y    x    y
%   0.00     0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2
% -----------------------------------
%
% ※ 複数の xlsx ファイルからまとめて解析したい場合は、
%    出力した１つのcsvファイルに単に列をコピーして追加していけばよい。
%    TrjNoはダブってよい。
% ----- 手作業でコピーした Output のファイル形式     
%            |> 1ファイル目のデータ         |> 2ファイル目のデータ
%   TrjNo    1    1    2    2    3    3    1    1 
%   x_or_y   x    y    x    y    x    y    x    y
%   0.00     0.0  0.0  1.3  2.3  0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1  0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2  0.4  0.3  1.2  1.2
% -----------------------------------
%
% ※ 解析する長さを短くしたい場合は、
%　  解析したくない時間以降の行を削除したcsvを用意してください。
%    例：eu2-1733.csv > 31frame目以降を削除して、eu2-1733_frm1-30.csv


% <<ピクセルサイズ、時間間隔、ファイル名を入力して Cntrol + Enter を押して実行>>
pitch = 0.0338; % [um/pix]
interval = 0.551; % [s]

fileName = 'eu2-1733_488';

filePath = [pwd filesep fileName '.xlsx'];
t = readtable(filePath);
trjNoList = unique(t.ID);

% テーブル読込
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

% Cell配列に成形
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

% csv データ書き出し
filePath = [pwd filesep fileName '.csv'];
writecell(Cell,filePath)
disp('CSVファイル書き出し完了')