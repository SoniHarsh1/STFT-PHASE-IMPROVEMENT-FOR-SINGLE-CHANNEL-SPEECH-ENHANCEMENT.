%ip audio signal
[ip, fs] = audioread('audio_file_path');

%stft of ip signal+baseband
[IP, ~, ~] = stft(ip, fs, 'Window', hamming(256, 'periodic'), 'OverlapLength', 32);
[p,q]=size(IP);
for k=1:p
    for l=1:q
        IP(k,l)=IP(k,l)*exp(k*l*-0.785i);
    end
end

%add noise to ip signal; n= ip+noise
ip1=ip;
signal_power = sum(ip.^2) / length(ip);
desired_snr_db = 5;
desired_snr_linear = 10^(desired_snr_db/10);
desired_noise_power = signal_power / desired_snr_linear;
noise = sqrt(desired_noise_power) * randn(size(ip));
n = ip + noise;
n1=n;

%convert n into stft
[N, F, T] = stft(n, fs, 'Window', hamming(256, 'periodic'), 'OverlapLength', 32);

%|n| matrix and phase IP matrix
absN=abs(N);
phaseIP= angle(IP);

%phase enhancement, no. of harmonics have been taken as 5
[f0, time] = pitch(ip, fs);
average_f0 = mean(f0);
f0=average_f0;

 for k=1:p
    min=0.02453125*k-(2*3.14*f0/fs);
    for h= 2:5
        t=min-(2*3.14*h*f0/fs);
        if min>abs(t)
            min=t;
        end
    end
        for l=2:q
            phaseIP(k,l)=phaseIP(k,l-1)+(32*min-0.785*k);
        end
end

%rev baseband, enhanced signal construction
for k=1:p
    for l=1:q
        phaseIP(k,l)=phaseIP(k,l)+(k*l*0.785);
    end
end
OP=zeros(size(phaseIP));
for k=1:p
    for l=1:q
        OP(k,l)=absN(k,l)*exp(phaseIP(k,l)*1i);
    end
end

%IDFT, overlap and add
[op, ~] = istft(OP, fs, 'Window', hamming(128, 'periodic'), 'OverlapLength', 32);
op1=abs(op);

%calculate SNRs
ensignal_power = sum(op1.^2) / length(ip1);
en_snr = ensignal_power / desired_snr_linear;
signal_power = sum(ip1.^2) / length(op1);

%display results
disp('SNR of noisy signal');
disp(desired_snr_db);
disp('Improvement in SNR');
disp(ensignal_power/signal_power);

%plot spectrograms
subplot(3,1,1);
spectrogram(ip1, hamming(256), 32, 128, fs, 'yaxis');
title('Original Signal');
subplot(3,1,2);
spectrogram(n1, hamming(256), 32, 128, fs, 'yaxis');
title('Noisy Signal');
subplot(3,1,3);
spectrogram(op1, hamming(256), 32, 128, fs, 'yaxis');
title('Enhanced Signal');