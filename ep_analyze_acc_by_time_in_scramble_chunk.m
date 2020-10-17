%%Analysis instructions:

%Goal: test whether performance accuracy relates to predictability. 
%The 1B, 2B, and 8B scrambled conditions have some (quarter-note) beats that follow a beat within 
%the same musical chunk (e.g., beat 2 in a 1B measure is predicted by beat 1, but beat 1 isn’t 
%predicted by the previous beat 2). 
%Plot beat accuracy as a function of the length of the preceding, predictive context (e.g., 
%how many beats preceded this beat that were part of the same musical chunk). 
%Edit this script to do this. 

%Tips: 
%The “stimuli” folder (in the “Piano study analysis project” folder) contains audio files, 
%sheet music, and scramble orders for all conditions. 
%Important: The “scramble orders” files contain 2 columns each. Each row is a musical chunk. 
%Column 1 = the starting measure of a chunk in the scrambled piece and column 2 = 
%the starting measure in the original Intact piece that that chunk was pulled from. 
%You’ll notice that not all chunks in the 8B version are exactly 8 measures (this is not a big deal; 
%we had to do some fudging to keep musical phrases intact).

%Navigate to the folder on your computer where this file lives
cd('/Users/eap/Dropbox/fMRI_music/keyboard_main/behavior/MIDI');

%Load the beat-by-beat accuracy data
load('analyzed/prop_corr_all_conditions_beatpadding=0.mat');

%Note: each 3D matrix in prop_corr_all_conditions is 13x240x3, where the
%1st dim = subjects (N = 13), 2nd dim = beats (all stimuli have 240 beats,
%or 120 measures containing 2 quarter note beats each), 3rd dim = rep (each
%subject completed 3 repetitions of each condition). The conditions = 1B,
%2B, 8B, and I.

%Now, load the scramble order files for each condition and for each beat, 
%compute how many beats separate that beat from the previous scramble
%boundary. For example, for the first chunk of music in the 2B stimulus,
%the first beat has had 0 beats since a scramble boundary, the 2nd beat has
%had 1, the 3rd beat has had 2, the 4th beat has had 3, and the 5th beat
%has had 0 again...and so on.

%Then, for each condition, compute the accuracy according to how many beats 
%it's been since a scramble boundary. You should have some kind of bar or line 
%plot with the x axis = distance from the last scramble boundary (for 1B,
%this axis will have values of 0 and 1, for 2B it will have values of 0 1 2
%3, etc.) Y axis = average accuracy across all musical chunks. You can make
%separate plots for different reps if you want.

%I've started an example for you for 2B below:
load('../stimuli/scramble_orders/2B_scramble_order.mat');

%Create a vector that looks like this....
2B_beat_by_boundary = [1 2 3 4 1 2 3 4 1 2 3 4...]; 

acc_data_2B = prop_corr_all_conditions{2};

%...