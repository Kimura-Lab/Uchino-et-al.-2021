% �O�Չ�̓f�[�^ XLSX �t�@�C������ ��͗p��csv�`���ɕύX����B
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Input�`�� ���̃X�N���v�g�{��[ConvertTrajectoryData.m]�Ɠ����t�H���_���ɁA
%    [FileName].xlsx �Ƃ����ȉ��̌`����csv�t�@�C����ۑ��B
%    �t�@�C�����͔C�ӂŁA�ȉ��� fileName �p�����[�^�̒l�ɓ��͂���B
%
% 2) �ȉ��̃p�����[�^�[����͂���B
%     pitch [um/pix]
%     interval [s]
%
% ----- Output �̃t�@�C���`��
%   TrjNo    1    1    2    2  ...
%   x_or_y   x    y    x    y
%   0.00     0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2
% -----------------------------------
%
% �� ������ xlsx �t�@�C������܂Ƃ߂ĉ�͂������ꍇ�́A
%    �o�͂����P��csv�t�@�C���ɒP�ɗ���R�s�[���Ēǉ����Ă����΂悢�B
%    TrjNo�̓_�u���Ă悢�B
% ----- ���ƂŃR�s�[���� Output �̃t�@�C���`��     
%            |> 1�t�@�C���ڂ̃f�[�^         |> 2�t�@�C���ڂ̃f�[�^
%   TrjNo    1    1    2    2    3    3    1    1 
%   x_or_y   x    y    x    y    x    y    x    y
%   0.00     0.0  0.0  1.3  2.3  0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1  0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2  0.4  0.3  1.2  1.2
% -----------------------------------
%
% �� ��͂��钷����Z���������ꍇ�́A
%�@  ��͂������Ȃ����Ԉȍ~�̍s���폜����csv��p�ӂ��Ă��������B
%    ��Feu2-1733.csv > 31frame�ڈȍ~���폜���āAeu2-1733_frm1-30.csv


% <<�s�N�Z���T�C�Y�A���ԊԊu�A�t�@�C��������͂��� Cntrol + Enter �������Ď��s>>
pitch = 0.0338; % [um/pix]
interval = 0.551; % [s]

fileName = 'eu2-1733_488';

filePath = [pwd filesep fileName '.xlsx'];
t = readtable(filePath);
trjNoList = unique(t.ID);

% �e�[�u���Ǎ�
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

% Cell�z��ɐ��`
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

% csv �f�[�^�����o��
filePath = [pwd filesep fileName '.csv'];
writecell(Cell,filePath)
disp('CSV�t�@�C�������o������')