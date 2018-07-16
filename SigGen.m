% M sequence
% mseq18.m
% ��������M�n���C�ӌƐ�������TSP�M����C�ӌȂ���D
clear;
format compact;

%% m�n�񐶐�
%https://astamuse.com/ja/drawing/JP/2009/124/499/A/000013.png
% m = 18; %����
% a = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % �����l
% h = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1]; % ���n������
%m = 20; %����:1-26�͈̔͂Ȃ�gfprimdf����������ȊO�͎�ł�
MLS_deg = input('M�n��̎���(2��n��)��1�`26�͈͓̔��Ŏw��:');
MLS_link = input('M�n��̘A�������w��F');
TSP_deg = input('TSP�̎���(2��n��)���w�肵�Ă�������:');
TSP_link = input('TSP�M���̘A�������w��F');

h = gfprimdf(MLS_deg);h = fliplr(h);
h = h(2:MLS_deg+1);
a = zeros(size(h));a(MLS_deg) = 1;

mdata = zeros(1,(2^MLS_deg-1)); % for���[�v�������̂��ߏ�����

for n=1:(2^MLS_deg-1)
    ms = rem(sum(a.*h),2);
	a = [ms a(1:MLS_deg-1)];
	% 0/1�f�[�^�� -1/1�f�[�^�ɕς��Ċi�[
	if ms == 0
        mdata(n) = -1;
    else
        mdata(n) = ms;
	end
end
% ���ݑ��֋��߂�Ƃ��Ɏg������t�@�C���ԒʐM
fileID = fopen(['M:\Research\SigGen\MLS\mls',num2str(MLS_deg),'.txt'],'w');
fprintf(fileID,'%23.15e\n',mdata);
fclose(fileID);
% disp('mdata�M�����������݂܂��D.');
% fileID = fopen('m20.dat','w');
% fwrite(fileID,mdata,'float');
% fclose(fileID);
%% �A�b�v�T���v��
% �݌v�d�l
% ���[���I�t�t�B���^�̎�ށ@	
% ���[���I�t�����@	0.5
% �I�[�o�[�T���v����	 5
% �ł��؂萔    �@    3
% ���K���̗L���@	�Ȃ�
% 
% disp('�R�T�C�����[���I�t�t�B���^�f�[�^��ǂݍ��݂܂�.');
% fileID = fopen('cosfilter.txt','r');
% filter = fscanf(fileID,'%e\n',Inf);
% fclose(fileID);
% 
% oversample = 5;
% mdata_up = upsample(mdata,oversample);
% mdata_filter = conv(mdata_up,filter);
% 
% fileID = fopen('mf18.txt','w');
% fprintf(fileID,'%e\n',mdata_filter);
% fclose(fileID);
%% TSP�M��
% �w�b�_�[���Ƃ����^�C�v�̃f�[�^������tsp�M���e�L�X�g�t�@�C��
% TSP�M���̒���

TSP_n = 2^TSP_deg;
%N = 262144;
m_tsp = TSP_n / 4;

% �ʑ������
k1 = 0 : TSP_n / 2;
omega1 = 4 * pi * m_tsp * (k1 / TSP_n) .^ 2;
k2 = TSP_n / 2 + 1 : TSP_n - 1;
omega2 = -4 * pi * m_tsp * ((TSP_n - k2) / TSP_n) .^ 2;
omega = [omega1 omega2];

% �tFFT   % �X�y�N�g�������
tspdata = ifft(exp(0 + 1i * omega));

tspdata = real(tspdata); % �v�Z�덷�ŋ������c��̂Ŏ��������
tspdata = circshift(tspdata, [0, -m_tsp]); % �g�`�̎n�܂�����炷
% �v���}�C�P�ɐ��K��
tsp_nw=fix(1/max(abs(tspdata)));
tspdata = tspdata * tsp_nw; 

fileID = fopen(['M:\Research\SigGen\TSP\tsp',num2str(TSP_deg),'N.txt'],'w');
fprintf(fileID,'%23.15e\n',tspdata);
fclose(fileID);

%% �����������ăo�C�i����������
% m�n��f�[�^������

mls_data = mdata;%�P��
for n=1:MLS_link-1%�Q����
    mls_data = horzcat(mls_data,mdata);
end
% tsp�M��������

tsp_data = tspdata;%�P��
for n=1:TSP_link-1%2����
    tsp_data = horzcat(tsp_data,tspdata);
end
data = horzcat(mls_data,tsp_data);
% disp('tspdata�M�����������݂܂��D.');
% fileID = fopen('tsp32.dat','w');
% fwrite(fileID,tsp_data,'float');
% fclose(fileID);
plot(data)

disp('���M�p�M�����������݂܂��D.');
fileID = fopen(['M:\Research\SigGen\result\mls',num2str(MLS_deg),'-',num2str(MLS_link),'tsp',num2str(TSP_deg),'-',num2str(TSP_link),'N.dat'],'w');
fwrite(fileID,data,'float');
fclose(fileID);

fclose('all');
disp('�I��');