% �O�Ղ�xy���W�̃��X�g���畽��MSD���Z�o����Anomalous diffusion model �� fitting����B
%                                          Ver.1 by Yuma Ito 2021.9.13
%
% 1) Curve Fitting Toolbox ���K�v�B[�A�v��] > [����ɃA�v�����擾] ����C���X�g�[���B
% 2) Input�`�� ���̃X�N���v�g�{��[Script_CalcFitAnomalousMSD.m]�Ɠ����t�H���_���ɁA
%    [FileName].csv�Ƃ����ȉ��̌`����csv�t�@�C����ۑ��B
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
% �� �Z���O�Ղ�x,y���W�͋󗓂ł����܂�Ȃ����A����MSD��fitting�͍Œ��Ōv�Z�����B
% �@�@Fitting���钷����Z���������ꍇ�́AFitting�������Ȃ����Ԉȍ~�̍s���폜����csv��p�ӂ��Ă��������B
% �� �P�ʂ͎��� s, ���� um �ɂ��낦�Ă���csv�ɂ��Ă��������B�O�Ղ̊J�n�� 0 s �Ƃ��܂��B
% �� �e�O�Ղ̎n�_�� (0,0) �łȂ��ėǂ��B�v���O������(0,0)�ɕ��s�ړ�����B
 
% ----- Output �̃t�@�C���`��
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
%   ���ԊԊu�� 1/precision �ɂ��� Fit���f���̋Ȑ��B

% <<�t�@�C��������͂��� Cntrol + Enter �������Ď��s>>
fileName = 'eu2-1733_frm1-20'
CalcFitAnnomMSD(fileName);


% ----- �ȉ��͊֐� -----
function CalcFitAnnomMSD(fileName)
    % �S�Ă̏��������s����֐�
    
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
        disp('���ԊԊu���ψ�ł͂���܂���B������MSD���v�Z�ł��܂���B')
    end
end

function msdList = CalcMSD(trjList)
    % ������MSD�v�Z�@
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
    %�@���v�ʂ̌v�Z
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
    
    % ���v�ʂ��܂Ƃ߂�
    dataList = [timeList' msdAve' msdSD' msdSEM' msdN'];
    
    % �X��MSD��ǉ�����B
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
    
    % Table�ɂ��ďo�͂���
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
