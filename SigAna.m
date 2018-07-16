%-------------------------------------------------
% GNU Radioから受け取ったバイナリファイルからインパルス応答を求めるプログラム
% 送信信号はSigGen.mで作られたもののみ対応．
% 受信信号のm系列を元に切り出し位置算出後，tsp信号の最初と最後を除いて同期加算する．
%-------------------------------------------------

%% main関数
function SigAna(  )
%% 初期化
clear
close all
format compact
%% フォント
set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','游ゴシック');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','游ゴシック');
%% 初期条件
% fc = input('中心周波数を入力してください:');%中心周波数
% samplerate = input('サンプリング周波数を入力してください:');
fc = 927e6;    samplerate = 20e6;
disp(['中心周波数',num2str(fc*1e-6),'MHz / サンプリング周波数',num2str(samplerate*1e-6),'MHz  に設定'])
%% ファイル読み込み

disp('復調ファイル');
[data,datname] = readBinary;
    function [ data, fname ] = readBinary(  )
        %readBinary
        %   バイナリファイル(.dat)を読み込む
        %   引数：なし
        %   戻り値：テキストファイル内のデータ．
        %   テキストファイルのデータ形式：1行目にデータ個数，2行目以降にデータ．
        
        % ファイルダイアログ表示
        [fname, dpath] = uigetfile({'*.dat','datファイル(*.dat)';
            '*.*','すべてのファイル(*.*)'}...
            ,'ファイルの選択');
        display(fname);
        fprintf('\n');
        % フルパス名を得る
        fullname = fullfile(dpath, fname);
        [~,fname,~]=fileparts(fullname);
        
        fp = fopen(fullname,'r');
        data = fread(fp,Inf,'float');  % ファイル終端まで数値データを読み込む
        fclose('all');                       % ファイルを全部閉じる
    end

%     TSP_deg = input('TSP信号の次数を入力：');
TSP_deg = 18;
disp(['TSP信号の次数は',num2str(TSP_deg),'に設定'])

%% M系列の相互相関位置を求めTSP信号を同期加算

hakei = m_tsp_anly; %次の内部関数へ
%     disp('ケーブル直結の信号について')     %ケーブル直結分を引き去る処理が必要なときにEnable
%     cable_rxdata = m_tsp_anly;

    function  hakei = m_tsp_anly
        %m_tsp_anly 受信信号をm系列を元に一つのTSP信号受信波形にする関数
        %   outに同期加算したTSP信号が出てくる．
        %   引数：なし
        %   戻り値：同期加算したtsp信号
        
        % MLS_deg = input('M系列の次数を入力：');
        MLS_deg = 18;
        disp(['M系列の次数は',num2str(MLS_deg),'に設定'])
        fp = fopen(['M:\Research\SigGen\MLS\mls',num2str(MLS_deg),'.txt'],'r');
        mdata = fscanf(fp,'%f\n');  % ファイル終端まで数値データを読み込む
        fclose('all');              % ファイルを全部閉じる
        
        % MLS_link = input('M系列の連結数を入力：');
        % TSP_link = input('M系列の連結数を入力：');
        MLS_link = 6;
        TSP_link = 16;
        disp(['M系列の連結個数は',num2str(MLS_link),'コ，TSP信号は',num2str(TSP_link),'コに設定'])
        disp('|')
        %% 切り出し
        % m系列データの点数
        mdata_n = length(mdata);
        % 受信データの点数
        % data_n = length(data);
        % 相互相関
        sogo_start = mdata_n; % whileループで更新していく．
        [f,lag] = xcorr(data(sogo_start:sogo_start*1+mdata_n*13),mdata);%M系列の元データ
        lag_n = length(lag);
        
        % ループ前処理(do whileにしようとした)
        search_f = 'a';
        cond = '0';
        LR = round(lag_n/2);
        lag1 = 1;
        lag2 = LR;
        cut_num = 3;    % MLS_linkより小さい値を入れる．
        % yが入力されるまで回す
        while cond ~= 'y'
            if cond == 'n'      %検索範囲
                disp([num2str(cut_num),'番目の指定'])
                fprintf('現在は開始:%d  終了:%d に設定されています\n',lag1,lag2);
                disp('lag1はラグmの検索 開始 位置を指定 デフォルト:1');
                lag1 = input('lag1:');
                disp('lag2はラグmの検索 終了 位置を指定 デフォルト:LR');
                lag2 = input('lag2:');
                fprintf('現在は開始:%d  終了:%d に設定されています\n',lag1,lag2);
            end
            if search_f == 'a'  %最大値検索
                disp('最大値検索');
                disp([num2str(cut_num),'番目の指定'])
                [f_value,f_index] = max(f(LR+lag1 :LR+lag2 -1 ));
                % lag1とlag2の調整を反映
                f_index = f_index + lag1 -1 ;
                fprintf('最大値:%d\n',f_value);
                fprintf('先頭から%d番目\n',f_index);
                
            elseif search_f ==  'b'     %最小値検索
                disp('最小値検索');
                disp([num2str(cut_num),'番目の指定'])
                [f_value,f_index] = min(f(LR+lag1 :LR+lag2 -1 ));
                % lag1とlag2の調整を反映
                f_index = f_index + lag1 -1 ;
                fprintf('最小値:%d\n',f_value);
                fprintf('先頭から%d番目\n',f_index);
                
            end
            
            figure(1)
            %     plot(lag,f,'',f_max_index,f_max,'*')
            plot(lag,f,'',f_index,f_value,'*')
            title('18次 M系列')
            xlabel('ラグm')
            ylabel('相互相関値')
            grid on
            
            cond = input('[y:OK n:検索範囲を変更 / a:最大値検索に変更 b:最小値検索に変更 p:指定箇所変更]:','s');
            if cond == 'a'
                search_f = 'a';
            elseif cond == 'b'
                search_f = 'b';
            elseif cond == 'p'
                % 連結個数分ある最大値/最小値から何番目かを指定する．
                cut_num = input(['指定箇所変更[1～',num2str(MLS_link),']:']);
            end
            disp('|')
        end
        %---------------------------------------------------------------１００点前へ
        % cut_num番目の先頭を検索したから．
        xstart = sogo_start + f_index + mdata_n*(MLS_link-cut_num+1);% - 100; %次にfor文の回数をわかりやすくするために+mdata処理
        figure
        % plot(data(xstart-2^17:xstart+2^18*35))
        plot(data(xstart-2^17:xstart+2^18*16-1))
        
        % 同期加算
        adddata=zeros(2^TSP_deg,TSP_link-2);
        for i=2:TSP_link-1  %最初と最後はノイズ乗りやすい
            adddata(:,i-1) = data((2^TSP_deg)*(i)+xstart : 1 : (2^TSP_deg)*(i+1)+xstart-1 , 1);% for文のiを個目にしたためdataの指定が複雑になった．
            if i==2 || i==3 || i==TSP_link-2 || i==TSP_link-1
                figure
                plot(adddata(:,i-1))
                title([num2str(i),'番目に切り出した受信信号']);
                xlabel('データ数[番目]')
                ylabel('大きさ')
            end
        end
        hakei=mean(adddata,2); %各行の平均値
        figure
        plot(hakei)
        title('同期加算したした受信信号');
        xlabel('データ数[番目]')
        ylabel('大きさ')
        % おかしかったらここで停止．
        input('処理を続行します．[push:ENTER]：','s');
        
        %     figure
        %     plot(out)
        %     title('同期加算波形の表示');
        
        %writeText(out);
        
    end

clear data

%% TSP信号を畳み込みしてインパルス応答を求める．
disp('TSP信号読み込み中')
fp = fopen(['M:\Research\SigGen\TSP\tsp',num2str(TSP_deg),'N.txt'],'r');%正規化済みデータ
fscanf(fp,'%d',1);%点数
tsp = textscan(fp,'%f',Inf);%data
tsp = tsp{1,1}(:,1);
fclose('all');

disp('畳み込み演算中')
imp = tspim(tsp,hakei);

    function [ imp ] = tspim( tsp , rec )
        %tspim tsp信号と録音信号からインパルス応答を推定する関数
        %   引数：tsp信号，録音した波形
        %   戻り値：推定インパルス応答
        out = conv(flipud(tsp), rec);% 時間軸反転したTSP信号と録音波形を畳み込み
        imp = out(length(tsp) : 1 : length(out)); % 畳み込み遅延量を補正して切り出す
        %     figure
        %     %writeText(imp);
        %     plot(imp)
        %     title('TSP')      
    end

Xaxis=0:1:length(imp)-1;
figure
subplot(2,1,1);
plot(Xaxis,imp(:,1));
title('インパルス応答全体の表示');
xlabel('データ数[番目]')
ylabel('大きさ')

Mlength = 1000;
xstart = 1;
xend=xstart+Mlength-1;
Main=imp(xstart:1:xend,1);

subplot(2,1,2)
plot(Main)
title('切り出したインパルス応答の表示');
xlabel('データ数[番目]')
ylabel('大きさ')

%% インパルス応答の書き込み
disp('インパルス応答ファイルimp')
%ファイル出力
[fname, dpath] = uiputfile({'*.txt','txtファイル(*.txt)';
    '*.*','すべてのファイル(*.*)'}...
    ,'特性出力ファイルの選択',...
    ['imp_',datname,'.txt']);
% フルパス名を得る
fname = fullfile(dpath, fname);
display(fname);

fp = fopen(fname,'w');
% 書き込み - 個数とデータ
fprintf(fp,'%d\n',numel(imp));
fprintf(fp,'%23.15e\n',imp);
fclose('all');
clear tsp

%% インパルス応答から振幅，位相，群遅延特性を算出
disp('周波数特性算出中')
[AM1,f1,GD1] = chractor(Main);
%     imp2 = tspim(tsp,cable_rxdata);　%直結との差分を用いる用
%     [~,~,GD2] = chractor(h,imp2);
    function [ AM, f, GD ] = chractor(Main)
        %chractor インパルス応答から周波数特性を計算する関数
        %   引数：インパルス応答
        %   戻り値：なし
        h = 65536; %解析精度:刻み幅(2^12)
        disp(['周波数特性の解析精度：',num2str(h)])
        
        %% 周波数特性計算
        [AM,f] =freqz(Main,1,h,1);
        GD = grpdelay(Main,1,h,1);
        
        %% 描画
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
%% 周波数特性を書き込み
disp('特性ファイルChr')
%ファイル出力
[fname, dpath] = uiputfile({'*.txt','txtファイル(*.txt)';
    '*.*','すべてのファイル(*.*)'}...
    ,'特性出力ファイルの選択',...
    ['chr_',datname,'.txt']);
% フルパス名を得る
fname = fullfile(dpath, fname);
display(fname);

%fp = fopen(fname,'w');

%データセット
f_real = f1*samplerate+fc; % サンプリング周波数->実周波数
AM_dB = 20*log10(abs(AM1));% 振幅特性をdB表示に
PH = 360/(2*pi)*angle(AM1);% 位相特性をdegree表示で
%     GD_smp = GD1-GD2;    % 直結時の分をひくときに用いる
GD_smp = GD1;
GD_sec = GD_smp/samplerate;% 群遅延特性の単位sample->sec

T=table(f1,f_real,AM_dB,PH,GD_smp,GD_sec);% テーブルセット
writetable(T, fname,'Delimiter','\t'); % 書き込み

fclose('all');

end


%% 処理確認時に外部値を利用するため
% function [ samples , data ] = readText(  )
% %readText テキストファイル読み込み関数
% %   テキストファイル(.txt)を読み込む
% %   引数：なし
% %   戻り値：テキストファイル内のデータ．
% %   テキストファイルのデータ形式：1行目にデータ個数，2行目以降にデータ．
%     % ファイルダイアログ表示
%     [fname, dpath] = uigetfile({'*.txt','txtファイル(*.txt)';
%                                 '*.*','すべてのファイル(*.*)'}...
%                                 ,'入力ファイルの選択');
%     % フルパス名を得る
%     fname = fullfile(dpath, fname);
%     fp = fopen(fname,'r');
%     display(fname);
%     samples = fscanf(fp,'%d',1);
%     data = textscan(fp,'%f',samples);
%
%     % ファイル閉じる
%     fclose('all');
%
% end
%% 処理確認時に外部出力するため
% function [  ] = writeText( data )
% %writeText この関数の概要をここに記述
% %   詳細説明をここに記述
%     disp('書き込みファイル');
%
%     % ファイルダイアログ表示
%     [fname, dpath] = uiputfile({'*.txt','TEXTファイル(*.txt)';
%                                 '*.*','すべてのファイル(*.*)'}...
%                                 ,'ファイルの選択');
%     display(fname);
%     % フルパス名を得る
%     fname = fullfile(dpath, fname);
%     fp = fopen(fname,'w');
%
%     % 書き込み - 個数とデータ
%     fprintf(fp,'%d\n',numel(data));
%     fprintf(fp,'%23.15e\n',data);
%     close('all');
% end
