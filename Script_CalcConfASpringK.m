% �O�Ղ�xy���W�̃��X�g���� Area of Confinement �� Spring Constant ���v�Z����B
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Input�`�� 
%    ���̃X�N���v�g�{��[Calc_Fit_Anomalous_MSD.m]�Ɠ����t�H���_���ɁA
%    FileName.csv�Ƃ����ȉ��̌`����csv�t�@�C����ۑ��B
%    �t�@�C�����͔C�ӂŁA�ȉ��� fileName �p�����[�^�̒l�ɓ��͂���B
%
% ----- Input csv �̃t�@�C���`�� -----
%   TrjNo    1    1    2    2  ...
%   x_or_y   x    y    x    y
%   0.00     0.0  0.0  1.3  2.3
%   0.03     0.2  0.2  0.3  3.1
%   0.06     0.4  0.3  1.2  1.2
% -----------------------------------
% �� 1��ڂ�3�s�ڂ���́A�X�e�b�v���Ԃ���͂��Ă��������B
% �� 2��ڈȍ~�́A�O�Ղ��ƂɁA���̎��Ԃɂ�����x,y���W����͂��Ă��������B
% �� �Z���O�Ղ�x,y���W�͋󗓂ł����܂�Ȃ��B
% �� �P�ʂ͎��� s, ���� um �ɂ��낦�Ă���csv�ɂ��Ă��������B�O�Ղ̊J�n�� 0 s �Ƃ��܂��B
% �� �e�O�Ղ̎n�_�� (0,0) �łȂ��ėǂ��B�v���O������(0,0)�ɕ��s�ړ�����B
 
% ----- Output �̃t�@�C���`��
%   FileName_ConfA_SpringK.csv
%   TrjNo  Rlong Rshort ConfA SpringK         
%   1      0.0   0.0    1.3   2.3
%   2      0.2   0.2    0.3   3.1
%   3      0.4   0.3    1.2   1.2
%   
%   TrjNo: �O�Քԍ��̘A�ԁBInput csv�t�@�C���̔ԍ��ɂ͑Ή������A�P�ɏ���
%   Rlong: Trajectory Point �� 95%���܂܂��ȉ~�̒����a (�P�ʂ� um)
%   Rshort: �ȉ~�̒Z���a (�P�ʂ� um)
%   ConfA: Trajectory Point �� 95%���܂܂��ȉ~�̖ʐ� (�P�ʂ� um^2)
%   SpringK: �΂˒萔�B(�P�ʂ� kbT/um^2)

% Reference
% Confinment area :
%   Germier, T. et al. Biophys J 113, 1383?1394 (2017).
%   http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix
%
% Spring coefficient : 
%   Shukron, O., Seeber, A., Amitai, A. & Holcman, D.
%   Trends Genet 35, 685?705 (2019).

% <<�t�@�C��������͂��� Cntrol + Enter �������Ď��s>>
fileName = 'eu2-1733_frm1-30'
CalcConfASpringK(fileName);


% ----- �ȉ��͊֐� -----
function CalcConfASpringK(fileName)
    % �S�Ă̏��������s����֐�
    [timeList,trjList] = GetData(fileName);
    CheckTime(timeList)
    [rlList, rsList, confAList] = CalcConfA(trjList);
    kList = CalcSpringK(timeList, trjList);
    SaveConfASpringK(rlList,rsList,confAList,kList,fileName);
end

function [timeList,trjList] = GetData(fileName)
    % Csv�t�@�C����Ǎ��݁ACell�z��ɂ���
    filePath = [pwd filesep fileName '.csv']
    fid = fopen(filePath);
    rows = textscan(fid,'%s');
    fclose(fid);
    rowCell = cellfun(@(x) split(x,','), rows,'UniformOutput',false);
    rowCell = rowCell{1};
    
    % ���ԍs�����o��
    timeCell = rowCell(3:end,1)';
    timeList = cellfun(@str2double,timeCell);
    
    % �e�O�Ղ�xy���W�����o��
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
        
        %�����l����
        x = rmmissing(x);
        y = rmmissing(y);
        
        % �O�Ղ̍ŏ��̍��W��(0,0)�ɂ���
        x = x - x(1);
        y = y - y(1); 
        trjList{trjNo} = [x y];
    end
end

function CheckTime(timeList)
    % ���Ԃ����Ԋu���m�F����
    %  ���������_�덷��10^-15���x����̂ŁA���̍�(=0)��10^10�{���Ďl�̌ܓ����A
    %  0�ɂȂ邩�ǂ����Ŕ��f����
    if length(unique(round(diff(diff(timeList))*10^10))) > 1
        disp('���ԊԊu���ψ�ł͂���܂���B������SpringK���v�Z�ł��܂���B')
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
    %   Shukron et al Trends Genet 2019 �� Box 2 �Q��
    % - D * k * x + b �ƂȂ�̂� �X�� a = -D*k ����� k = -a/D
    
    kList = {};
    for trjNo = 1:length(trjList)
        trj = trjList{trjNo};
        xTrj = trj(:,1); yTrj = trj(:,2);
        dx = diff(xTrj); dy = diff(yTrj);
        xc = mean(xTrj); yc = mean(yTrj);
        
        % x = Distance from centroid (1-dimentional)
        % Displace�ƑΉ�����̂ōŌ�̂P�_������
        x = [(xTrj(1:end-1)-xc); (yTrj(1:end-1)-yc)];
        
        % y = velocity = displacement/interavl (1-dimentional)
        dt = mean(diff(timeList));
        y = [dx; dy]/dt;
        
        % ���`��A
        X = [ones(length(x),1) x];
        b = X\y;
        
        % k = -�X��/D, (1step displacement^2 = 4*D*t)
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
    
    % Table�ɂ��ďo�͂���
    t = array2table(dataList);
    Header = {'TrjNo','Rlong','Rshort','ConfA','SpringK'};
    t.Properties.VariableNames = Header;
    writetable(t,filePath)
    disp('ConfA SpringK exported')
end
