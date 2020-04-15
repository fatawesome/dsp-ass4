clear()
clf()
// Read signal
//[s, Fs, _] = wavread("signal_with_low_freq_noise_2.wav")

load("signal_with_noise_and_filtered.sod");
// Take the first channel
s = signal_with_noise(1, :)
// Plot
//plot(s)
//xlabel("Time, n", 'fontsize', 2)
//ylabel("Amplitude", 'fontsize', 2)
//title("Signal with noise in time domain", 'fontsize', 3) 
//xs2png(0, "signal_with_noise_time.png")
 
// Plot spectrum
// calculate frequencies

//s_len = length(s)
//frequencies = (0:s_len-1)/s_len * Fs;
//
//// plot
//
//plot2d("nl", frequencies, abs(fft(s)), color("blue")) 
//xlabel("Frequency component, n", 'fontsize', 2)
//ylabel("Freq. Amplitude", 'fontsize', 2)
//title("Signal with noise in frequency domain", 'fontsize', 3) 
//xs2png(0, "signal_with_noise_freq.png")


// N: integer length of FIR filter
// cutoff: fraction of Fs, at which frequencies are stopped 
// stop_value: the value for frequencies in the stop band 
//   (after cutoff frequency)
// return: frequency representation of an ideal
//   low pass FIR filter of length N+1 if N is even
//   or N if N is odd
function H = ideal_lowpass(N, cutoff, stop_value)
    N = (N - modulo(N, 2)) / 2
    cutoff = floor(2 * N * cutoff)
    H = ones(1, N) * stop_value
    H(1, 1:cutoff) = 1
    printf("%d, %d\n", N, cutoff)
    // need to make N odd
    H = [1. H flipdim(H, 2)]
endfunction

function H = highpass(N, cutoff, stop_value)
    N = (N - modulo(N, 2)) / 2
    cutoff = floor(2 * N * cutoff)
    H = ones(1, N)
    H(1, 1:cutoff) = stop_value
    
    printf("%d, %d\n", N, cutoff)
    // need to make N odd
    H = [1. H flipdim(H, 2)]
endfunction

function H = task1_filter(N, cutoff_low, cutoff_high, stop_low, stop_high)
    N = (N - modulo(N, 2)) / 2
    cutoff_low = floor(2 * N * cutoff_low)
    cutoff_high = floor(2 * N * cutoff_high)
    printf("%d, %d, %d\n", N, cutoff_low, cutoff_high)
    H = ones(1, N)
    H(1, 1:cutoff_high) = stop_high
    H(1, cutoff_low:N) = stop_low
    H(1, 1:30) = 0
    
    // need to make N odd
    H = [1. H flipdim(H, 2)]    
endfunction

function sv = shift(v,n)
if (size(v,'r')<>1 & size(v,'c')<>1) then
    error("1st argument must be column vector or row vector")
end
if size(v,'r')<>1 then
    v = shift(v',n)
    sv = v'
else
    n = modulo(n,size(v,'c'))
    sv = v($-n+1:$)
    sv = [sv v(1:$-n)]
end
endfunction

H_l = ideal_lowpass(256, 0.15, 0.);
h_len = length(H_l)
frequencies = (0:h_len-1)/h_len * Fs;
// Compute impulse response
// project into temporal domain
// imaginary part should be close to 0
h_l = real(ifft(H_l))
mid = floor(h_len / 2)
//h_l_shifted = [h_l(:, mid + 1:h_len), h_l(:, 1:mid)]
h_l_shifted = shift(h_l, int(h_len / 2))
h_l_windowed = h_l_shifted .* window('kr', length(h_l), 8)

disp(h_l_windowed)

//plot2d('nn', 0:length(h_l_windowed)-1, h_l_windowed, color("blue"))
//xlabel("Time, n", 'fontsize', 2)
//ylabel("Amplitude", 'fontsize', 2)
//title("Impulse response of ideal low-pass filter", 'fontsize', 3)
//xs2png(gcf(), "ideal_lowpass_time.png")

////plot2d("nl", frequencies, abs(fft(h_l_windowed)), color("blue"))
////xlabel("Frequency, Nz", 'fontsize', 2)
////ylabel("Freq amplitude", 'fontsize', 2)
////title("Frequency response of the final FIR filter", 'fontsize', 3)
////xs2png(gcf(), "frequency_response_final_fir_filter.png")

H_h = highpass(256, 0.001, 0.);
h_len = length(H_l)
frequencies = (0:h_len-1)/h_len * Fs;
h_h = real(ifft(H_h))
mid = floor(h_len / 2)
h_h_shifted = shift(h_h, int(h_len / 2))
h_h_windowed = h_h_shifted .* window('kr', length(h_h), 8)

task1_fir = task1_filter(8192, 0.15, 30 / 44100, 0., 0.)
task1_fir_len = length(task1_fir)

frequencies = (0:task1_fir_len-1)/task1_fir_len * Fs;

// квадтратный график
//plot2d("nn", frequencies, task1_fir, color("blue"))
//xlabel("Frequency, Hz", 'fontsize', 2)
//ylabel("Freq amplitude", 'fontsize', 2)
//title("Frequency response of band-pass task1 filter", 'fontsize', 3)
//xs2png(0, "task1-band-pass-freq.png")

task1_fir = real(ifft(task1_fir))
task1_fir_shifted = shift(task1_fir, int(task1_fir_len / 2))
task1_fir_windowed = task1_fir_shifted .* window('kr', length(task1_fir_shifted), 8)

//plot2d('nl', frequencies, abs(fft(task1_fir_windowed)), color("blue"))
//xlabel("Time, n", 'fontsize', 2)
//ylabel("Amplitude", 'fontsize', 2)
//title("Impulse response of ideal low-pass filter", 'fontsize', 3)
//xs2png(gcf(), "ideal_lowpass_time.png")

//zopa = convol(h_h_windowed, convol(h_l_windowed, s))
//savewave('task1-piped-filter.wav', zopa, 44100)
//zopa = convol(ffilt("bp", 256, 9000, 10000), s)       
zopa = convol(task1_fir_windowed, s)
savewave('task1-combined-filter.wav', zopa, 44100)

frequencies = (0:length(zopa)-1)/length(zopa) * Fs;

//plot2d("nl", frequencies, abs(fft(zopa)), color("blue"))
//xlabel("Frequency, Nz", 'fontsize', 2)
//ylabel("Freq amplitude", 'fontsize', 2)
//title("Frequency response of the final FIR filter", 'fontsize', 3)
//xs2png(gcf(), "frequency_response_final_fir_filter.png")


//--------------

h_len = length(task1_fir);
//frequencies = (0:h_len-1)/h_len * Fs;
//plot2d("nn", frequencies, task1_fir, color("blue"));
//xlabel("Frequency, Hz", 'fontsize', 2);
//ylabel("Freq amplitude", 'fontsize', 2);
//title("Frequency response of ideal low-pass filter", 'fontsize', 3);

//h_l = real(ifft(task1_fir));
//plot2d('nn', 0:length(h_l)-1, h_l, color("blue"));
//xlabel("Time, n", 'fontsize', 2);
//ylabel("Amplitude", 'fontsize', 2);
//title("Impulse response of ideal low-pass filter", 'fontsize', 3);

task1_fir_shifted = shift(task1_fir, int(task1_fir_len / 2))
//plot2d('nn', 0:length(task1_fir_shifted)-1, task1_fir_shifted , color("blue"));
//xlabel("Time, n", 'fontsize', 2);
//ylabel("Amplitude", 'fontsize', 2);
//title("Impulse response of shifted task1 band-pass filter", 'fontsize', 3);

task1_fir_windowed = task1_fir_shifted .* window('kr', length(task1_fir_shifted), 8)
plot2d('nn', 0:length(task1_fir_windowed)-1, task1_fir_windowed, color("blue"));
xlabel("Time, n", 'fontsize', 2);
ylabel("Amplitude", 'fontsize', 2);
title("Impulse response of shiffted task1 filter after window", 'fontsize', 3);

//frequencies = (0:task1_fir_len-1)/task1_fir_len * Fs;
//plot2d("nl", frequencies, abs(fft(task1_fir_windowed)), color("blue"));
//xlabel("Frequency, Hz", 'fontsize', 2);
//ylabel("Freq amplitude", 'fontsize', 2);
//title("Frequency response of the final FIR filter", 'fontsize', 3);
