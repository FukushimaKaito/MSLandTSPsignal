% M sequence
% mseq18.m
% 生成したM系列を任意個と生成したTSP信号を任意個つなげる．
clear;
format compact;

%% m系列生成
%https://astamuse.com/ja/drawing/JP/2009/124/499/A/000013.png
% m = 18; %次数
% a = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % 初期値
% h = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1]; % 原始多項式
%m = 20; %次数:1-26の範囲ならgfprimdfが動くそれ以外は手打ち
MLS_deg = input('M系列の次数(2のn乗)を1～26の範囲内で指定:');
MLS_link = input('M系列の連結個数を指定：');
TSP_deg = input('TSPの次数(2のn乗)を指定してください:');
TSP_link = input('TSP信号の連結個数を指定：');

h = gfprimdf(MLS_deg);h = fliplr(h);
h = h(2:MLS_deg+1);
a = zeros(size(h));a(MLS_deg) = 1;

mdata = zeros(1,(2^MLS_deg-1)); % forループ高速化のため初期化

for n=1:(2^MLS_deg-1)
    ms = rem(sum(a.*h),2);
	a = [ms a(1:MLS_deg-1)];
	% 0/1データを -1/1データに変えて格納
	if ms == 0
        mdata(n) = -1;
    else
        mdata(n) = ms;
	end
end
% 相互相関求めるときに使うやつをファイル間通信
fileID = fopen(['M:\Research\SigGen\MLS\mls',num2str(MLS_deg),'.txt'],'w');
fprintf(fileID,'%23.15e\n',mdata);
fclose(fileID);
% disp('mdata信号を書き込みます．.');
% fileID = fopen('m20.dat','w');
% fwrite(fileID,mdata,'float');
% fclose(fileID);
%% アップサンプル
% 設計仕様
% ロールオフフィルタの種類　	
% ロールオフ率α　	0.5
% オーバーサンプル比	 5
% 打ち切り数    　    3
% 正規化の有無　	なし
% 
% disp('コサインロールオフフィルタデータを読み込みます.');
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
%% TSP信号
% ヘッダーをとったタイプのデータだけのtsp信号テキストファイル
% TSP信号の長さ

TSP_n = 2^TSP_deg;
%N = 262144;
m_tsp = TSP_n / 4;

% 位相を作る
k1 = 0 : TSP_n / 2;
omega1 = 4 * pi * m_tsp * (k1 / TSP_n) .^ 2;
k2 = TSP_n / 2 + 1 : TSP_n - 1;
omega2 = -4 * pi * m_tsp * ((TSP_n - k2) / TSP_n) .^ 2;
omega = [omega1 omega2];

% 逆FFT   % スペクトルを作る
tspdata = ifft(exp(0 + 1i * omega));

tspdata = real(tspdata); % 計算誤差で虚部が残るので実部を取る
tspdata = circshift(tspdata, [0, -m_tsp]); % 波形の始まりをずらす
% プラマイ１に正規化
tsp_nw=fix(1/max(abs(tspdata)));
tspdata = tspdata * tsp_nw; 

fileID = fopen(['M:\Research\SigGen\TSP\tsp',num2str(TSP_deg),'N.txt'],'w');
fprintf(fileID,'%23.15e\n',tspdata);
fclose(fileID);

%% 水平結合してバイナリ書き込み
% m系列データを結合

mls_data = mdata;%１個目
for n=1:MLS_link-1%２から
    mls_data = horzcat(mls_data,mdata);
end
% tsp信号を結合

tsp_data = tspdata;%１個目
for n=1:TSP_link-1%2から
    tsp_data = horzcat(tsp_data,tspdata);
end
data = horzcat(mls_data,tsp_data);
% disp('tspdata信号を書き込みます．.');
% fileID = fopen('tsp32.dat','w');
% fwrite(fileID,tsp_data,'float');
% fclose(fileID);
plot(data)

disp('送信用信号を書き込みます．.');
fileID = fopen(['M:\Research\SigGen\result\mls',num2str(MLS_deg),'-',num2str(MLS_link),'tsp',num2str(TSP_deg),'-',num2str(TSP_link),'N.dat'],'w');
fwrite(fileID,data,'float');
fclose(fileID);

fclose('all');
disp('終了');