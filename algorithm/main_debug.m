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
% x = pilot;
% x  = x.';
x = pilot.';

%     Index = file_name(10:end);
%     IndexLen = length(Index);
%     for n = 1:NewNameLen-IndexLen
%         Index = [num2str(0),Index];
%     end
[1.6900    1.7000    1.7600    1.9000    1.9900    2.1100    2.2500    2.4400...,
    2.5400    2.5700    2.7600    2.8400    2.9000    4.0800    5.6100,...
    8.2400    8.9000   10.8500   66.8800, 72.5100   75.0600];

[ 763    27   127   163   521   767   242   708   251,...
    45   737   455   480   435   655   248   537    91];
Datanum = 655; % choose the data label
file_name = namelist(Datanum).name;
N = 816;L = [408,512];
Res = 1;
SNR = 3;
eigThreshold = 1e-4;
H_Sample_gap = 1;

data = load(file_name);
y = example_64Tc.';
y = data.y.';


DelayLen = 1024/Res;
Time = [0:Res:(DelayLen - 1)*Res];

H = y./x;

% for HH = 1:N-2*L+1 % all non-overlapping
%     for hh = HH + L:N-L+1
%         H_ = H_ + H(HH:HH+L-1,:)*H(hh:hh+L-1,:)';
%         Scalar = Scalar + 1;
%     end
% end

% for hh = 1:1:N - 2*L + 1  % original non-overlapping
%     H_ = H_ + H(hh:hh+L-1,:)*H(hh+L:hh+2*L-1,:)';
%     Scalar = Scalar + 1;
% end
% 
H_ = zeros(L(1),1);
Scalar = 0;
H = H*H';
for HH = 1:H_Sample_gap:N + 1-L(1) % postive smoothing
%     eigThreshold = 3e-3;
        H_ = H_ + H(HH:HH+L(1)-1,HH:HH+L(1)-1);
        Scalar = Scalar + 1;
end
F = fliplr(eye(N));                             % transpose matrix
H = F*(conj(H))*F;
for HH = 1:H_Sample_gap:N + 1-L(1) % negative smoothing
%     eigThreshold = 3e-3;
        H_ = H_ + H(HH:HH+L(1)-1,HH:HH+L(1)-1);
        Scalar = Scalar + 1;
end
H_ = 1/Scalar*H_;
% for HH = 1:N - 2*L + 1 % partial overlapping
% %     eigThreshold = 3e-3;
%     for hh = 0:H_Sample_gap
%         H_ = H_ + H(HH:HH+L-1,:)*H(HH+L-hh:HH+2*L-1-hh,:)';
%         Scalar = Scalar + 1;
%     end
% end

% for hh = 1:N - L + 1  % all overlapping
%     eigThreshold = 0.01;
%     H_ = H_ + H(hh:hh+L-1,:)*H(hh:hh+L-1,:)';
%     Scalar = Scalar + 1;
% end

% H_ = 1/Scalar*H_;
a = H_;
% a = H*H';
% a = H_*H_';

[U,D] = eig(a);
% D = abs(D);
D = diag(D)';
[D, I] = sort(D);
U = fliplr(U(:,I));
% D_dB = 10*log10(D/max(D));
[Lp(1),MDL] = LS_MDL(D,Scalar,L(1));
% Lp = sum( abs(D) >= eigThreshold);
% Lp = 4;

for kk = 1:DelayLen
    V = exp(-1j*2*pi*[0:length(D)-1]'*scs*comb*kk*Res*Tc);
    P_MUSIC(kk) = 1/abs((V'*U(:,Lp(1) + 1:end)*(V'*U(:,Lp(1) + 1:end))'));
end

P_MUSIC = 10*log10((P_MUSIC)/max(P_MUSIC));
[pks, pksid, w, p] = findpeaks(P_MUSIC);
% IND = min(pksid(P_MUSIC(pksid) - min(P_MUSIC) >= SNR));
flag = 0;

for findlos = 1:length(pksid)
    if  P_MUSIC(pksid(findlos)) - min(P_MUSIC) > SNR
        IND = pksid(findlos);
        flag = 1;
        break;
    end
end
if flag == 0
    IND = pksid(1);
end
% IND = pksid(min(sortedpksid(1:min(Lp,length(sortedpksid)))));
delay(1) = IND*Res;
subplot 311
plot(Time,P_MUSIC,'b');hold on;
% plot(Time(PeaksIndex),P_MUSIC(PeaksIndex),'bo'); hold on;
plot(Time(pksid),P_MUSIC(pksid),'k*');hold on;
plot(Time(IND),P_MUSIC(IND),'rs','markerface','r');hold on;
% compare terms---------------------------------------------------------
H_ = zeros(L(2),1);
Scalar = 0;
for HH = 1:N + 1-L(2) - H_Sample_gap % postive smoothing
%     eigThreshold = 3e-3;
        H_ = H_ + H(HH:HH+L(2)-1,HH:HH+L(2)-1);
        Scalar = Scalar + 1;
end
F = fliplr(eye(N));                             % transpose matrix
H = F*(conj(H))*F;
for HH = 1:N + 1-L(2) - H_Sample_gap % negative smoothing
%     eigThreshold = 3e-3;
        H_ = H_ + H(HH:HH+L(2)-1,HH:HH+L(2)-1);
        Scalar = Scalar + 1;
end
H_ = 1/Scalar*H_;

a = H_;

[U,D] = eig(a);
% D = abs(D);
D = diag(D)';
[D, I] = sort(D);
U = fliplr(U(:,I));
% D_dB = 10*log10(D/max(D));
[Lp(2),MDL] = LS_MDL(D,Scalar,L(2));
% Lp = sum( abs(D) >= eigThreshold);
% Lp = 4;

for kk = 1:DelayLen
    V = exp(-1j*2*pi*[0:length(D)-1]'*scs*comb*kk*Res*Tc);
    P_MUSIC(kk) = 1/abs((V'*U(:,Lp(2) + 1:end)*(V'*U(:,Lp(2) + 1:end))'));
end


P_MUSIC = 10*log10((P_MUSIC)/max(P_MUSIC));

[pks, pksid, w, p] = findpeaks(P_MUSIC);
sxy = y.*conj(x);
%         sxy = sxy./abs(sxy);
sxy = (ifft(sxy,length(y)));
cxy = sxy(:,1:end);

% IND = min(pksid(P_MUSIC(pksid) - min(P_MUSIC) >= SNR));
flag = 0;
[~, pksid_corr, ~, ~] = findpeaks(abs(sum(cxy(1:20,:),2)));
peaksind = min(pksid_corr)*T/Tc;

for findlos = 1:length(pksid)
    if  P_MUSIC(pksid(findlos)) - min(P_MUSIC) > SNR
        IND = pksid(findlos);
        flag = 1;
        break;
    end
end
if flag == 0
    IND = pksid(1);
end
% IND = pksid(min(sortedpksid(1:min(Lp,length(sortedpksid)))));
delay(2) = IND*Res;
delay
Lp
% delay = min(pksid(p>= mean(p)))*Res
subplot 312
plot(Time,P_MUSIC,'b');hold on;
% plot(Time(PeaksIndex),P_MUSIC(PeaksIndex),'bo'); hold on;
plot(Time(pksid),P_MUSIC(pksid),'k*');hold on;
plot(Time(IND),P_MUSIC(IND),'rs','markerface','r');hold on;




Timecorr = 0:T/Tc:T*(length(sum(cxy,2)) - 1)/Tc;
subplot 313
plot(Timecorr(1:20),abs(sum(cxy(1:20,:),2)),'b');hold on;
[sort_val, sort_ind] = sort(abs(cxy(2:end,:)), 'descend');
delay_sort = (min(sort_ind(sort_val >= mean(sort_val(1:SortNum_1)))));
[~, pksid_corr, ~, ~] = findpeaks(abs(sum(cxy(1:20,:),2)));
peaksind = min(pksid_corr);
plot(Timecorr(delay_sort),abs(sum(cxy(delay_sort,:),2)),'bs','markerface','b');hold on;
plot(Timecorr(pksid_corr),abs(sum(cxy(pksid_corr,:),2)),'k*');hold on;
plot(Timecorr(peaksind),abs(sum(cxy(peaksind,:),2)),'rs','markerface','r');hold on;