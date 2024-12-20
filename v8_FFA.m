% Fundamental Frequency Analyzer (MATLAB R2019a)
% This script analyzes all audio files in the current directory, plots
% the time-domain and frequency-domain for each, calculates the fundamental
% frequency using a manual FFT implementation, maps it to the corresponding
% musical note, and saves the figure with the same name as each audio file.

% Step 1: Get a list of all audio files in the current folder
audioFiles = dir('*.wav');  % Supported extensions are: '*.flac' '*.m4a' '*.mp4' '*.oga' '*.ogg' '*.wav'. replace if needed
if isempty(audioFiles)
    disp('No audio files found in the directory.');
    return;
end

% Loop through each audio file in the directory
for k = 1:length(audioFiles)
    % Get the audio file name and path
    audioFile = audioFiles(k).name;
    [~, audioName, ~] = fileparts(audioFile);  % Get the name without extension

    % Step 2: Read the audio file
    [audioData, fs] = audioread(audioFile);

    % Ensure mono audio (if stereo, average the channels)
    if size(audioData, 2) == 2
        audioData = mean(audioData, 2);
    end

    % Check if the audio duration is at least 3 seconds
    if length(audioData) / fs < 3
        disp(['Skipping ', audioFile, ': audio is too short (< 3 seconds).']);
        continue;
    end

    % Optional: Cut the audio to the first X seconds (default: disabled)
    % Uncomment the following lines to enable truncation
    % maxDuration = 3;  % Truncate to X seconds
    % audioData = audioData(1:min(length(audioData), fs * maxDuration));

    % Zero-pad the signal to the next power of 2
    N = 2^nextpow2(length(audioData));
    audioData = [audioData; zeros(N - length(audioData), 1)];

    % Step 3: Plot the time-domain representation
    t = (0:length(audioData)-1) / fs;  % Time vector

    % Create the figure and set the name of the window to the audio file name
    figure;
    set(gcf, 'Name', audioName, 'NumberTitle', 'off');

    % Time-domain plot
    subplot(2, 1, 1);
    plot(t, audioData);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Domain Representation');
    grid on;

    % Step 4: Compute the FFT using the manual implementation
    audioDataFFT = abs(my_fft_iterative(audioData));
    f = (0:N-1) * (fs / N);  % Frequency vector

    % Frequency-domain plot
    subplot(2, 1, 2);
    % Apply color coding for frequencies out of human hearing range
    if all(f(1:N/2) < 20 | f(1:N/2) > 20000)
        plot(f(1:N/2), audioDataFFT(1:N/2), 'Color', [0.5, 0, 0.5]);  % Purple
    else
        plot(f(1:N/2), audioDataFFT(1:N/2), 'b');  % Default color
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Frequency Domain (Manual FFT)');
    grid on;

    % Step 5: Find the fundamental frequency
    [~, idx] = max(audioDataFFT(1:N/2));  % Index of the peak
    fundamentalFreq = f(idx);  % The fundamental frequency (in Hz)

    % Step 6: Map the fundamental frequency to a musical note
    note = frequencyToNote(fundamentalFreq);

    % Annotate the plot with a red dashed line for the fundamental frequency
    hold on;
    line([fundamentalFreq fundamentalFreq], [0 max(audioDataFFT)], ...
        'Color', 'r', 'LineWidth', 2);  % Dashed red line 'LineStyle', '--',

    % Add legend with red line and note information
    legend(plot(fundamentalFreq, max(audioDataFFT), 'r-'), ...
           ['Fundamental: ', num2str(fundamentalFreq), ' Hz [', note, ']'], ...
           'Location', 'northeast', 'FontSize', 10);

    % Save the figure with the same name as the audio file
    saveas(gcf, [audioName, '.png']);  % Save the figure as a PNG file
end

% Custom FFT Implementation (Iterative)
function X = my_fft_iterative(x)
    N = length(x);
    if mod(N, 2) ~= 0
        error('Input length must be a power of 2.');
    end
    X = bitrevorder(x); % Rearrange the input in bit-reversed order
    stages = log2(N);
    for s = 1:stages
        m = 2^s;
        half_m = m / 2;
        W_m = exp(-2j * pi * (0:half_m-1) / m);
        for k = 1:m:N
            for j = 0:half_m-1
                t = W_m(j+1) * X(k + j + half_m);
                u = X(k + j);
                X(k + j) = u + t;
                X(k + j + half_m) = u - t;
            end
        end
    end
end

% Function to map frequency to note
function note = frequencyToNote(frequency)
    % Define the standard pitch A4 = 440 Hz
    A4 = 440;  
    notes = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};

    % Calculate the number of semitones from A4
    semitonesFromA4 = 12 * log2(frequency / A4);

    % Round the semitones to the nearest integer
    semitoneIndex = round(semitonesFromA4);

    % Map semitoneIndex to the corresponding note
    noteIndex = mod(semitoneIndex, 12) + 1;  % Wrap around to get a valid note index

    % Handle negative noteIndex correctly
    if noteIndex == 0
        noteIndex = 12;  % Ensure the note index is valid
    end

    % Return the note corresponding to the frequency
    note = notes{noteIndex};
end
