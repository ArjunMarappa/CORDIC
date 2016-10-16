% Cordic Implementation
% Author: Arjun Marappa
% Description: MATLAB Implementation of CORDIC as Digitally Controlled Oscillator. 

clc
clear all
close all

nsamples = 4096;
fftsize = 65536*7;
fsample = 100e3;
fsignal = 10e3;
frange = (-0.5:1/fftsize:0.5-1/fftsize)*fsample;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate I and Q ADC samples 
[I, Q] = TestSigGen(fsignal,fsample,nsamples,0);

Isample = fix(I*32767);
Qsample = fix(Q*32767);
x_in = (I+i*Q)*32767;

% Store the Inputs to ease of signal generation in future development iteration
fileID = fopen('C:\Users\Nagarjun\Documents\MATLAB\Thesis\cordicI.txt', 'w+');%%%%%%%%%%
for i= 1:nsamples
	fprintf(fileID, '%d\n',int32(Isample));
end
fclose(fileID); 
fileID = fopen('C:\Users\Nagarjun\Documents\MATLAB\Thesis\cordicQ.txt', 'w+');%%%%%%%%%%
for i= 1:nsamples
     fprintf(fileID, '%d\n',int32(Qsample));
end
fclose(fileID);

figure('NumberTitle', 'off',...
         'Name', 'Received Signal');
plot(0:nsamples-1,x_in)
title('\bfReceived Signal')
xlabel('\bfTime')
ylabel('\bfMagnitude')
% 
SigSpec = fftshift(fft(x_in,fftsize));
figure
plot(frange,dB(psdg(SigSpec/max(SigSpec))));
title('\bfFrequency Specturm Received Signal')
xlabel('\bfFrequency')
ylabel('\bfMagnitude in dB')
ylim([-80,10])
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_str = sprintf('Please enter a freq (less than %dHz) you want to generate',fsample); % DCO tuning frequency
NCO_Freq = input(input_str);

% double values cannot represent all integers greater than 2^53 correctly
% convert PhaseInc to int64.
PhaseInc = 2*(((NCO_Freq * 2^(32-1)) / fsample))
% dec2hex(int32(PhaseInc),8)
% PhaseInc = 536870912;% 32bit phase increment corresponding to pi/4
% PhaseInc = 6039797;
% PhaseInc = 10737418;

STG = 20;   % Number of rotations 
K = 0.6073; % Cordic Gain 
K = 0.60725294104140;

for i = 0:1:(23)
    c = (round((atan(1.0/(2^i))/(2*pi)) * (2^24)));
    consts(i+1) = c;
end
consts = int32(consts);
% phase = (round(((angle)/(2*pi)) * (2^24)));
phasepast = 0;
init_0padding = 2;
% PhaseInc = PhaseInc/256; % conversion 32bit -> 24bit
x = int32(1)*(2^(16)-1);% 24 bit input to cordic
y = int32(0);
Is = int32([zeros(1,init_0padding) ones(1,nsamples-init_0padding)].*32767);
Qs = int32([zeros(1,nsamples+init_0padding)]);
% Is = Is*2^8;
% Qs = Qs*2^8;
Isample = Isample*2^8;
Qsample = Qsample*2^8;
% Is = int32([zeros(1,init_0padding) Isample(1:nsamples-init_0padding)]);
% Qs = int32([zeros(1,init_0padding) Qsample(1:nsamples-init_0padding)]);
phase = (0);
initial = 1;
for N = 1:nsamples
%######################################################
%           CORDIC PROCESSOR
%######################################################
   x = sign_ext(Is(N),24,1);% 25 bit Cordic processor
   y = sign_ext(Qs(N),24,1);
% Phase accumulator
phase = phase + PhaseInc;
% phasepast = phase;
% ############------[32 bit Integer Wrap Around]------############
if (phase > 2147483647)% if > 2^31 - 1 roll it back to -ve's
    phase = phase - 4294967295; %(2^31-1 + 2^31)
else if(phase < -2147483648) % -2^31
        phase = phase + 4294967295;
    end 
end

zin_24 = floor(phase/256);% right shift by 8 since zwidth is 24
zin = int32(zin_24);

% Phase pre rotation since cordic is limited to +pi/2 to -pi/2
if (zin>=2^22 && zin<2^23)% interval between 90 to 180 degrees
    xpast = x;
    x = -y;
    y = xpast;
    temp = bitand(zin, 4194303);
    zin = bitand(zin,0);
    zin = bitor(zin,temp);
 else if(zin>=-2^23 && zin<-2^22)% interval between -180 to -90 degrees
        xpast = x;
        x = y;
        y = -xpast;
        zin = bitor(zin,-12582912);
     end
end

% 20 stage, 27 bit cordic pipeline 
x = sign_ext(x,25,2); % sign extend to accomodate bit growth
y = sign_ext(y,25,2);

j = int32(0);
d=1;
P_vec(N) = phase;
% z = phase;
while j < STG
    if (zin > 0)
        d = 1;
    else
        d = -1;
    end 
    xpast = x;
    x = xpast - (d*bitshift(y,-j));
    y = y + (d*bitshift(xpast,-j));
	j = j+1;
    zin = zin - (d*consts(j));
% ############------[27bit OverFlow Compensation]------############    
    if (x > 134217727)% if > 2^27 -1 roll it back to -ve's
        x = x - 268435455;
    else if(x < -134217728)
            x = x + 268435455;
        end 
    end
    if (y > 134217727)% if > 2^27 -1 roll it back to -ve's
        y = y - 268435455;
    else if(y < -134217728)
            y = y + 268435455;
        end 
    end
end	

% Multiply the cordic gain
% x = K*x;
% y = K*y;
% Clip the cordic output to 24 bit
X(N) = bitshift(x,-2);
Y(N) = bitshift(y,-2);
Z(N) = zin;
end
%% Save Real and Imaginary values to test.mat file.
X_fft = fftshift(fft(double(X), fftsize));
save('C:\Users\Nagarjun\Documents\MATLAB\Thesis\cordic_out.mat', 'X', 'Y','NCO_Freq', 'fsample','nsamples');
%% Time Domain Plots.
figure('NumberTitle', 'off',...
        'Name', 'CORDIC output');
plot(X,'b')
hold 
plot(Y,'g')
grid on
title('CORDIC time domain')
legend('Real Samples', 'Imag Samples')
xlabel('Time')
ylabel('Magnitute')
%% Frequency Plots
figure('NumberTitle','off',...
        'Name','Power Spectral Desity plot')
plot(frange, dB(psdg(X_fft/max(X_fft))))
axis([-fsample/2 fsample/2 -80 10])
grid on;
title('Frequency domain plot of NCO')
xlabel('Frequency')
ylabel('Magnitude in dB')

% fileID = fopen('C:\Users\Nagarjun\Desktop\MathWorks\Matfiles\cordicI2cic_TB.txt', 'w+');%%%%%%%%%%
% for i= 1:nsamples
%     fprintf(fileID, '%d => %d,\n',i,int32(X(i)));
% end
% fclose(fileID);
cordic_delay  = 21;
X = [zeros(1,cordic_delay) X];
fileID = fopen('C:\Users\Nagarjun\Documents\MATLAB\Thesis\cordicI2cic_TB.txt', 'w+');%%%%%%%%%%
for i= 1:nsamples
    fprintf(fileID, '%d\n',(X(i)));
end
fclose(fileID);

fileID = fopen('C:\Users\Nagarjun\Documents\MATLAB\Thesis\cordicQ2cic_TB.txt', 'w+');%%%%%%%%%%
for i= 1:nsamples
    fprintf(fileID, '%d\n',(Y(i)));
end
fclose(fileID);

