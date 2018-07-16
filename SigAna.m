%-------------------------------------------------
% GNU Radio����󂯎�����o�C�i���t�@�C������C���p���X���������߂�v���O����
% ���M�M����SigGen.m�ō��ꂽ���̂̂ݑΉ��D
% ��M�M����m�n������ɐ؂�o���ʒu�Z�o��Ctsp�M���̍ŏ��ƍŌ�������ē������Z����D
%-------------------------------------------------

%% main�֐�
function SigAna(  )
%% ������
clear
close all
format compact
%% �t�H���g
set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','���S�V�b�N');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','���S�V�b�N');
%% ��������
% fc = input('���S���g������͂��Ă�������:');%���S���g��
% samplerate = input('�T���v�����O���g������͂��Ă�������:');
fc = 927e6;    samplerate = 20e6;
disp(['���S���g��',num2str(fc*1e-6),'MHz / �T���v�����O���g��',num2str(samplerate*1e-6),'MHz  �ɐݒ�'])
%% �t�@�C���ǂݍ���

disp('�����t�@�C��');
[data,datname] = readBinary;
    function [ data, fname ] = readBinary(  )
        %readBinary
        %   �o�C�i���t�@�C��(.dat)��ǂݍ���
        %   �����F�Ȃ�
        %   �߂�l�F�e�L�X�g�t�@�C�����̃f�[�^�D
        %   �e�L�X�g�t�@�C���̃f�[�^�`���F1�s�ڂɃf�[�^���C2�s�ڈȍ~�Ƀf�[�^�D
        
        % �t�@�C���_�C�A���O�\��
        [fname, dpath] = uigetfile({'*.dat','dat�t�@�C��(*.dat)';
            '*.*','���ׂẴt�@�C��(*.*)'}...
            ,'�t�@�C���̑I��');
        display(fname);
        fprintf('\n');
        % �t���p�X���𓾂�
        fullname = fullfile(dpath, fname);
        [~,fname,~]=fileparts(fullname);
        
        fp = fopen(fullname,'r');
        data = fread(fp,Inf,'float');  % �t�@�C���I�[�܂Ő��l�f�[�^��ǂݍ���
        fclose('all');                       % �t�@�C����S������
    end

%     TSP_deg = input('TSP�M���̎�������́F');
TSP_deg = 18;
disp(['TSP�M���̎�����',num2str(TSP_deg),'�ɐݒ�'])

%% M�n��̑��ݑ��ֈʒu������TSP�M���𓯊����Z

hakei = m_tsp_anly; %���̓����֐���
%     disp('�P�[�u�������̐M���ɂ���')     %�P�[�u�����������������鏈�����K�v�ȂƂ���Enable
%     cable_rxdata = m_tsp_anly;

    function  hakei = m_tsp_anly
        %m_tsp_anly ��M�M����m�n������Ɉ��TSP�M����M�g�`�ɂ���֐�
        %   out�ɓ������Z����TSP�M�����o�Ă���D
        %   �����F�Ȃ�
        %   �߂�l�F�������Z����tsp�M��
        
        % MLS_deg = input('M�n��̎�������́F');
        MLS_deg = 18;
        disp(['M�n��̎�����',num2str(MLS_deg),'�ɐݒ�'])
        fp = fopen(['M:\Research\SigGen\MLS\mls',num2str(MLS_deg),'.txt'],'r');
        mdata = fscanf(fp,'%f\n');  % �t�@�C���I�[�܂Ő��l�f�[�^��ǂݍ���
        fclose('all');              % �t�@�C����S������
        
        % MLS_link = input('M�n��̘A��������́F');
        % TSP_link = input('M�n��̘A��������́F');
        MLS_link = 6;
        TSP_link = 16;
        disp(['M�n��̘A������',num2str(MLS_link),'�R�CTSP�M����',num2str(TSP_link),'�R�ɐݒ�'])
        disp('|')
        %% �؂�o��
        % m�n��f�[�^�̓_��
        mdata_n = length(mdata);
        % ��M�f�[�^�̓_��
        % data_n = length(data);
        % ���ݑ���
        sogo_start = mdata_n; % while���[�v�ōX�V���Ă����D
        [f,lag] = xcorr(data(sogo_start:sogo_start*1+mdata_n*13),mdata);%M�n��̌��f�[�^
        lag_n = length(lag);
        
        % ���[�v�O����(do while�ɂ��悤�Ƃ���)
        search_f = 'a';
        cond = '0';
        LR = round(lag_n/2);
        lag1 = 1;
        lag2 = LR;
        cut_num = 3;    % MLS_link��菬�����l������D
        % y�����͂����܂ŉ�
        while cond ~= 'y'
            if cond == 'n'      %�����͈�
                disp([num2str(cut_num),'�Ԗڂ̎w��'])
                fprintf('���݂͊J�n:%d  �I��:%d �ɐݒ肳��Ă��܂�\n',lag1,lag2);
                disp('lag1�̓��Om�̌��� �J�n �ʒu���w�� �f�t�H���g:1');
                lag1 = input('lag1:');
                disp('lag2�̓��Om�̌��� �I�� �ʒu���w�� �f�t�H���g:LR');
                lag2 = input('lag2:');
                fprintf('���݂͊J�n:%d  �I��:%d �ɐݒ肳��Ă��܂�\n',lag1,lag2);
            end
            if search_f == 'a'  %�ő�l����
                disp('�ő�l����');
                disp([num2str(cut_num),'�Ԗڂ̎w��'])
                [f_value,f_index] = max(f(LR+lag1 :LR+lag2 -1 ));
                % lag1��lag2�̒����𔽉f
                f_index = f_index + lag1 -1 ;
                fprintf('�ő�l:%d\n',f_value);
                fprintf('�擪����%d�Ԗ�\n',f_index);
                
            elseif search_f ==  'b'     %�ŏ��l����
                disp('�ŏ��l����');
                disp([num2str(cut_num),'�Ԗڂ̎w��'])
                [f_value,f_index] = min(f(LR+lag1 :LR+lag2 -1 ));
                % lag1��lag2�̒����𔽉f
                f_index = f_index + lag1 -1 ;
                fprintf('�ŏ��l:%d\n',f_value);
                fprintf('�擪����%d�Ԗ�\n',f_index);
                
            end
            
            figure(1)
            %     plot(lag,f,'',f_max_index,f_max,'*')
            plot(lag,f,'',f_index,f_value,'*')
            title('18�� M�n��')
            xlabel('���Om')
            ylabel('���ݑ��֒l')
            grid on
            
            cond = input('[y:OK n:�����͈͂�ύX / a:�ő�l�����ɕύX b:�ŏ��l�����ɕύX p:�w��ӏ��ύX]:','s');
            if cond == 'a'
                search_f = 'a';
            elseif cond == 'b'
                search_f = 'b';
            elseif cond == 'p'
                % �A����������ő�l/�ŏ��l���牽�Ԗڂ����w�肷��D
                cut_num = input(['�w��ӏ��ύX[1�`',num2str(MLS_link),']:']);
            end
            disp('|')
        end
        %---------------------------------------------------------------�P�O�O�_�O��
        % cut_num�Ԗڂ̐擪��������������D
        xstart = sogo_start + f_index + mdata_n*(MLS_link-cut_num+1);% - 100; %����for���̉񐔂��킩��₷�����邽�߂�+mdata����
        figure
        % plot(data(xstart-2^17:xstart+2^18*35))
        plot(data(xstart-2^17:xstart+2^18*16-1))
        
        % �������Z
        adddata=zeros(2^TSP_deg,TSP_link-2);
        for i=2:TSP_link-1  %�ŏ��ƍŌ�̓m�C�Y���₷��
            adddata(:,i-1) = data((2^TSP_deg)*(i)+xstart : 1 : (2^TSP_deg)*(i+1)+xstart-1 , 1);% for����i���ڂɂ�������data�̎w�肪���G�ɂȂ����D
            if i==2 || i==3 || i==TSP_link-2 || i==TSP_link-1
                figure
                plot(adddata(:,i-1))
                title([num2str(i),'�Ԗڂɐ؂�o������M�M��']);
                xlabel('�f�[�^��[�Ԗ�]')
                ylabel('�傫��')
            end
        end
        hakei=mean(adddata,2); %�e�s�̕��ϒl
        figure
        plot(hakei)
        title('�������Z����������M�M��');
        xlabel('�f�[�^��[�Ԗ�]')
        ylabel('�傫��')
        % �������������炱���Œ�~�D
        input('�����𑱍s���܂��D[push:ENTER]�F','s');
        
        %     figure
        %     plot(out)
        %     title('�������Z�g�`�̕\��');
        
        %writeText(out);
        
    end

clear data

%% TSP�M������ݍ��݂��ăC���p���X���������߂�D
disp('TSP�M���ǂݍ��ݒ�')
fp = fopen(['M:\Research\SigGen\TSP\tsp',num2str(TSP_deg),'N.txt'],'r');%���K���ς݃f�[�^
fscanf(fp,'%d',1);%�_��
tsp = textscan(fp,'%f',Inf);%data
tsp = tsp{1,1}(:,1);
fclose('all');

disp('��ݍ��݉��Z��')
imp = tspim(tsp,hakei);

    function [ imp ] = tspim( tsp , rec )
        %tspim tsp�M���Ƙ^���M������C���p���X�����𐄒肷��֐�
        %   �����Ftsp�M���C�^�������g�`
        %   �߂�l�F����C���p���X����
        out = conv(flipud(tsp), rec);% ���Ԏ����]����TSP�M���Ƙ^���g�`����ݍ���
        imp = out(length(tsp) : 1 : length(out)); % ��ݍ��ݒx���ʂ�␳���Đ؂�o��
        %     figure
        %     %writeText(imp);
        %     plot(imp)
        %     title('TSP')      
    end

Xaxis=0:1:length(imp)-1;
figure
subplot(2,1,1);
plot(Xaxis,imp(:,1));
title('�C���p���X�����S�̂̕\��');
xlabel('�f�[�^��[�Ԗ�]')
ylabel('�傫��')

Mlength = 1000;
xstart = 1;
xend=xstart+Mlength-1;
Main=imp(xstart:1:xend,1);

subplot(2,1,2)
plot(Main)
title('�؂�o�����C���p���X�����̕\��');
xlabel('�f�[�^��[�Ԗ�]')
ylabel('�傫��')

%% �C���p���X�����̏�������
disp('�C���p���X�����t�@�C��imp')
%�t�@�C���o��
[fname, dpath] = uiputfile({'*.txt','txt�t�@�C��(*.txt)';
    '*.*','���ׂẴt�@�C��(*.*)'}...
    ,'�����o�̓t�@�C���̑I��',...
    ['imp_',datname,'.txt']);
% �t���p�X���𓾂�
fname = fullfile(dpath, fname);
display(fname);

fp = fopen(fname,'w');
% �������� - ���ƃf�[�^
fprintf(fp,'%d\n',numel(imp));
fprintf(fp,'%23.15e\n',imp);
fclose('all');
clear tsp

%% �C���p���X��������U���C�ʑ��C�Q�x���������Z�o
disp('���g�������Z�o��')
[AM1,f1,GD1] = chractor(Main);
%     imp2 = tspim(tsp,cable_rxdata);�@%�����Ƃ̍�����p����p
%     [~,~,GD2] = chractor(h,imp2);
    function [ AM, f, GD ] = chractor(Main)
        %chractor �C���p���X����������g���������v�Z����֐�
        %   �����F�C���p���X����
        %   �߂�l�F�Ȃ�
        h = 65536; %��͐��x:���ݕ�(2^12)
        disp(['���g�������̉�͐��x�F',num2str(h)])
        
        %% ���g�������v�Z
        [AM,f] =freqz(Main,1,h,1);
        GD = grpdelay(Main,1,h,1);
        
        %% �`��
        figure
        subplot(3,1,1)
        plot(f,20*log10(abs(AM)))
        title('Amplitude')
        xlabel('Normalized frequency')
        ylabel('Amplitude[dB]')
        
        subplot(3,1,2);
        plot(f,360/(2*pi)*angle(AM))
        title('Phase')
        xlabel('Normalized frequency')
        ylabel('Phase[degree]')
        
        subplot(3,1,3)
        plot(f,GD)
        title('GroupDelay')
        xlabel('Normalized frequency')
        ylabel('GroupDelay[sample]')
        
    end
%% ���g����������������
disp('�����t�@�C��Chr')
%�t�@�C���o��
[fname, dpath] = uiputfile({'*.txt','txt�t�@�C��(*.txt)';
    '*.*','���ׂẴt�@�C��(*.*)'}...
    ,'�����o�̓t�@�C���̑I��',...
    ['chr_',datname,'.txt']);
% �t���p�X���𓾂�
fname = fullfile(dpath, fname);
display(fname);

%fp = fopen(fname,'w');

%�f�[�^�Z�b�g
f_real = f1*samplerate+fc; % �T���v�����O���g��->�����g��
AM_dB = 20*log10(abs(AM1));% �U��������dB�\����
PH = 360/(2*pi)*angle(AM1);% �ʑ�������degree�\����
%     GD_smp = GD1-GD2;    % �������̕����Ђ��Ƃ��ɗp����
GD_smp = GD1;
GD_sec = GD_smp/samplerate;% �Q�x�������̒P��sample->sec

T=table(f1,f_real,AM_dB,PH,GD_smp,GD_sec);% �e�[�u���Z�b�g
writetable(T, fname,'Delimiter','\t'); % ��������

fclose('all');

end


%% �����m�F���ɊO���l�𗘗p���邽��
% function [ samples , data ] = readText(  )
% %readText �e�L�X�g�t�@�C���ǂݍ��݊֐�
% %   �e�L�X�g�t�@�C��(.txt)��ǂݍ���
% %   �����F�Ȃ�
% %   �߂�l�F�e�L�X�g�t�@�C�����̃f�[�^�D
% %   �e�L�X�g�t�@�C���̃f�[�^�`���F1�s�ڂɃf�[�^���C2�s�ڈȍ~�Ƀf�[�^�D
%     % �t�@�C���_�C�A���O�\��
%     [fname, dpath] = uigetfile({'*.txt','txt�t�@�C��(*.txt)';
%                                 '*.*','���ׂẴt�@�C��(*.*)'}...
%                                 ,'���̓t�@�C���̑I��');
%     % �t���p�X���𓾂�
%     fname = fullfile(dpath, fname);
%     fp = fopen(fname,'r');
%     display(fname);
%     samples = fscanf(fp,'%d',1);
%     data = textscan(fp,'%f',samples);
%
%     % �t�@�C������
%     fclose('all');
%
% end
%% �����m�F���ɊO���o�͂��邽��
% function [  ] = writeText( data )
% %writeText ���̊֐��̊T�v�������ɋL�q
% %   �ڍא����������ɋL�q
%     disp('�������݃t�@�C��');
%
%     % �t�@�C���_�C�A���O�\��
%     [fname, dpath] = uiputfile({'*.txt','TEXT�t�@�C��(*.txt)';
%                                 '*.*','���ׂẴt�@�C��(*.*)'}...
%                                 ,'�t�@�C���̑I��');
%     display(fname);
%     % �t���p�X���𓾂�
%     fname = fullfile(dpath, fname);
%     fp = fopen(fname,'w');
%
%     % �������� - ���ƃf�[�^
%     fprintf(fp,'%d\n',numel(data));
%     fprintf(fp,'%23.15e\n',data);
%     close('all');
% end
