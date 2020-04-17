clc; clear; close all

g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];
k = 8;                  % длина сообщения
r = length(g_x) - 1;    % степень g(x)
num_messages = 2^k;     % колво сообщений
[H, G] = hammgen(4);    % 

messages = de2bi((0:num_messages-1)');     
codewords = zeros(num_messages, k+r);      % массив кодовых слов
mes_3_11 = zeros(num_messages, 3, 11);     % х3 по 11
                                           
ham_3_15 = zeros(num_messages, 3, 15);     
channel_3_12 = zeros(num_messages, 3, 12); %  BPSK (-1, 1)
                                          
                        
for j=1:num_messages
    [~, c] = gfdeconv([zeros(1, r), messages(j,:)], g_x);   % CRC 
    codewords(j,:) = xor([zeros(1, r), messages(j,:)], ... 
                         [c, zeros(1, k+r - length(c))]);   
    mes_3_11(j, :, :) = [reshape(codewords(j, :), 8, 3)' ...
                         [0 0 0; 0 0 0; 0 0 0]];            
    for i=1:3
        ham_3_15(j, i, :) = mod(reshape(mes_3_11(j, i, :), 1, 11) * G, 2); % кодирование по хэммингу
        channel_3_12(j, i, :) = ham_3_15(j, i, 1:end-3)*(-2) + 1; 
    end
end

all_codewords = de2bi((0:(2^8-1))');
h = zeros(2^8, 12);
for j=1:2^8
    temp = mod([all_codewords(j, :) 0 0 0] * G, 2);
    h(j, :) = temp(1, 1:end-3) .*(-2) + 1;
end

% Теоретические значения
SNRdB = -10:0;
SNR = 10.^(SNRdB./10);
Pe_bit_theor = qfunc(sqrt(2.*SNR));         
Ped_theor_up = ones(1, length(SNR)).*1/2^r; 
Ped_theor = zeros(1,length(SNRdB));
w = sum(codewords(2:end, :)');
d = min(w);
for i=1:length(SNR)
    Pe_theoretical = 0;
    for j=d:(k + r)
        Pe_theoretical = Pe_theoretical + sum(w == j) * ...
                (Pe_bit_theor(i)^j) * ((1 - Pe_bit_theor(i))^((k + r)-j));
    end
    Ped_theor(i) = Pe_theoretical;
end

 
Pe_bit = zeros(1,length(SNRdB));  % ошибка на бит
Pe_bit_after_Hamming = zeros(1,length(SNRdB));
Pe_bit_after_Hamming_soft = zeros(1,length(SNRdB));
Ped = zeros(1,length(SNRdB));     % вер ошибки декодера
Ped_soft = zeros(1,length(SNRdB));
T = zeros(1,length(SNRdB));       % пропускная способность канала
for i=1:length(SNR)
    disp(i);
    sigma = sqrt(1/(2*SNR(i)));
    nTests = 0; nSent = 0; %колво попыток передачи %сколько было отправлено
    nErrDecode = 0; nErrDecode_soft = 0;
    nErrBits = 0; nErrBits_H = 0; nErrBits_H_soft = 0;
    messages_to_send = 0;
    while messages_to_send < 10*(SNR(i)+10)
        messages_to_send = messages_to_send + 1;
        rand_ind = randi(num_messages, 1);
        c = reshape(channel_3_12(rand_ind, :, :), 3, 12);
        while 1
            nTests = nTests + 1;
            AWGN = c + sigma*randn(3, 12); %сообщение плюс шум
            
            [corrected, nErr, v] = correctError(AWGN, H, c);    %жесткий
            nErrBits = nErrBits + nErr;
            nErrBits_H = nErrBits_H + v;
            [soft_corrected, v_soft] = ...                      %мягкий
                correctSoft(AWGN, h, all_codewords, c);
            nErrBits_H_soft = nErrBits_H_soft + v_soft;
            
            [~, S] = gfdeconv(corrected, g_x); 
            v = sum(xor(codewords(rand_ind, :), corrected));
            
            [~, S_soft] = gfdeconv(soft_corrected, g_x);
            v_soft = sum(xor(codewords(rand_ind, :), soft_corrected));
            
            if (bi2de(S_soft) == 0) && (v_soft > 0)
                nErrDecode_soft = nErrDecode_soft + 1;
            end
            
            if (bi2de(S) == 0)
                if (v > 0)
                    nErrDecode = nErrDecode + 1;
                end
                nSent = nSent + 1;
                break;
            end
        end
    end
    Ped(i) = nErrDecode/nTests;
    Ped_soft(i) = nErrDecode_soft/nTests;
    Pe_bit(i) = nErrBits/(nTests*36);
    Pe_bit_after_Hamming(i) = nErrBits_H/(nTests*36);
    Pe_bit_after_Hamming_soft(i) = nErrBits_H_soft/(nTests*36);
    T(i) = k * nSent / (36 * nTests); %пропускная способность канала
end

figure;
subplot(2, 2, 1);
semilogy(SNRdB, Ped, 'ko', ...
         SNRdB, Ped_soft, 'bo', ...
         SNRdB, Ped_theor_up, 'r.-', ...
         SNRdB, Ped_theor, 'b.-')
legend('Ped', 'Ped soft', 'Ped theor up', 'Ped theor');
subplot(2, 2, 3);
semilogy(SNRdB, Pe_bit, 'k', ...
         SNRdB, Pe_bit_after_Hamming, 'r.--', ...
         SNRdB, Pe_bit_after_Hamming_soft, 'r.-', ...
         SNRdB, Pe_bit_theor, 'b.-');
legend('Pe bit', 'Pe bit after Hamming', ...
       'Pe bit after Hamming soft', 'Pe bit theor');
xlabel('E/N_0, dB')
subplot(1, 2, 2);
plot(SNRdB, T, 'b.-');
ylabel('T, Пропускная способность');
xlabel('E/N_0, dB')

% 
function [concatenated, nErrBits, v] = correctError(AWGN, H, c)
	unBPSK = AWGN < 0;
    nErrBits = sum(sum(xor(c < 0, unBPSK))); %колво ошибочных битов в слове
    plus_zeros = [unBPSK(:, :) [0 0 0; 0 0 0; 0 0 0]];
    corrected = zeros(3, 12);

            % [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
            %  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
            %  0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
            %  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
            %  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
            %  0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
            %  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
            %  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
            %  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
            %  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
            %  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
            %  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
            %  0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
            %  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];

    for part=1:3
        corr = plus_zeros(part, :); %
                                    
        S = bi2de(mod(corr*H',2)); % синдром
                                  
        switch S
            case 1
                corr = gfadd(corr,[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
            case 2
                corr = gfadd(corr,[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
            case 3
                corr = gfadd(corr,[0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]);
            case 4
                corr = gfadd(corr,[0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]);
            case 5
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]);
            case 6
                corr = gfadd(corr,[0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]);
            case 7
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]);
            case 8
                corr = gfadd(corr,[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]);
            case 9
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);
            case 10
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]);
            case 11
                corr = gfadd(corr,[0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);
            case 12
                corr = gfadd(corr,[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]);
            case 13
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);
            case 14
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]);
            case 15
                corr = gfadd(corr,[0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);
        end
        corrected(part, :) = corr(1, 1:end-3);
    end
    v = sum(sum(xor(corrected, c<0)));
    concatenated = [corrected(1, 5:end) ...
                    corrected(2, 5:end) ...
                    corrected(3, 5:end)];
end

% Soft
function [concatenated, v] = correctSoft(AWGN, h, codes, c)
    min_d = [100 100 100];
    corrected = codes(1, :);
    cor = zeros(3, 12);
    for j=1:256
        for AWGN_part=1:3
            min = sqrt(sum((AWGN(AWGN_part, :) - h(j, :)).^2));
            if min < min_d(AWGN_part)
                min_d(AWGN_part) = min;
                corrected(AWGN_part, :) = codes(j, :);
                cor(AWGN_part, :) = h(j, :);
            end
        end
    end
    v = sum(sum(xor(cor<0, c<0)));
    concatenated = [corrected(1, :) ...
                    corrected(2, :) ...
                    corrected(3, :)];
end
