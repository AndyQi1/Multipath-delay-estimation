%% LoS identification, Ding, Qian 2023/4/28

clear all
close all
load pilot.mat % transmit SRS signal frequency-domain data
load example_64Tc.mat % LoS path with no noise

Tc = 1/(480 * 1000 * 4096); % minimum time slot (0.509ns)
fc = 2565e6; % carrier frequency (2565MHz)
B = 100e6; % bandwidth (100MHz)
scs = 30e3; % gap between subcarriers (30kHz)
dx = 5e-2;dy = 5e-2; % antenna size (5cm)
comb = 4; % number of combs
T = 64*Tc/3; % symbol time length according to example_64Tc

addpath '..\data'
addpath '.\raw data'
namelist = dir('.\raw data\*.mat');
%
len = length(namelist);
NewNameLen = 8;
SortNum_1 = 16;
SortNum_4 = 18;
x = pilot;
% x  = x.';
x = pilot.';
y = example_64Tc.';
y = data.y.';
H = y./x;
a = H*H';
% a = ifft(a);
[U,D] = eig(a);
D = diag(D)';
[D, I] = sort(D);
U = fliplr(U(:,I));
Lp = 4;
Res = Tc/3;
for kk = 1:256*3
    V = exp(-1j*2*pi*[0:815]'*scs*comb*kk*Res);
    P_MUSIC(kk) = 1/abs((V'*U(:,Lp + 1:end)*(V'*U(:,Lp + 1:end))'));
end
[v,id] = max(P_MUSIC);
id*Res/Tc
[local, pks, w, p] = findpeaks(P_MUSIC)
plot(P_MUSIC)

for i = 1:len
    
    file_name = namelist(i).name;
    %     Index = file_name(10:end);
    %     IndexLen = length(Index);
    %     for n = 1:NewNameLen-IndexLen
    %         Index = [num2str(0),Index];
    %     end
    data = load(file_name);
    y = data.y;
    %     if i <= 400
    %         y = data.ant1_data;
    %     else
    %         y = data.ant4_data;
    %     end
    %     save(Index,'y');
    %     y = example_64Tc;
    if i<=400
        
%         y = y*exp(-1j*2*pi*fc);
        
%         y = fft(filter(Lowpassfilter,ifft(y))); % Low-pass filter
        sxy = y.*conj(x);
%         sxy = sxy./abs(sxy); % GCC path
        sxy = (ifft(sxy));
        
        cxy = sxy;
        %     [val, Ind] = max(abs(corr));
        [~, max_ind] = max(abs(cxy));
        delay_max(i) = (max_ind - 1)*T/Tc; % 56 is the delay of FIR filter
        
        [sort_val, sort_ind] = sort(abs(cxy(2:end)), 'descend');
        delay_sort(i) = (min(sort_ind(sort_val >= mean(sort_val(1:SortNum_1)))))*T/Tc;
%         delay_sort(i) = (min(sort_ind(sort_val >= 0.2*sort_val(1))) - 1)*T/Tc;
        
        %     plot(lag, abs(corr));
%             plot(abs(cxy(1:50)))
    else
%         y = y*exp(-1j*(2*pi*(fc - 2*pi*dx/2*fc/1e8)));
%         y = fft(filter(Lowpassfilter,ifft(y,length(y),2),2),length(y),2); % Low-pass filter
        sxy = y.*conj(x);
%         sxy = sxy./abs(sxy);
        sxy = (ifft(sxy,length(y),2));
        cxy = sxy(:,1:end);
        [~, max_ind] = max(sum(abs(cxy)));
        delay_max(i) = (max_ind - 1)*T/Tc; % 56 is the delay of FIR filter
        
        [sort_val, sort_ind] = sort(sum(abs(cxy(:,2:end))), 'descend');
        delay_sort(i) = (min(sort_ind(sort_val >= mean(sort_val(1:SortNum_4)))))*T/Tc;
%         delay_sort(i) = (min(sort_ind(sort_val >= 0.2*sort_val(1))) - 1)*T/Tc;
        %         plot(abs(sum(cxy)));
    end
end
% delay_sort - delay_max

%% write into txt file

fid = fopen(['..\data\','answer.txt'],'w');
for kk = 1:length(delay_sort)
    if kk < length(delay_sort)
        fprintf(fid,'%.2f,\n',delay_sort(kk));
    else
        fprintf(fid,'%.2f',delay_sort(kk));
    end
end
fclose(fid);

