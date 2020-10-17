%ep_analyzeMIDIcsv
%Note: in 1B, 2B, in beat 1 (time 0), there's stuff

clear;
% cd('/Users/elise/Dropbox/fMRI_music/keyboard_main/behavior/MIDI/');
cd('/Users/eap/Dropbox/fMRI_music/keyboard_main/behavior/MIDI/');

subjects = [2 3 4 5 8 9 10 15 17 20 21 23 22]; nSubs = length(subjects); 
groups = {'AM', 'AM', 'M', 'M', 'M', 'M', 'AM', 'AM', 'M', 'AM', 'M', 'AM', 'M'};

conditions = {'1B', '2B', '8B', 'I'}; nCond = length(conditions);

completed_run1 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1];

completed_run2 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1];

completed_run3 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1];

all_completed_runs = cat(3, completed_run1, completed_run2, completed_run3); nRuns = 3;

n_samples_in_beat = 480; %For 120 measures, 2 beats/measure, this should be 480 (115200/120/2)
beat_padding = 0; %accept events that occur +/- this value outside of a beat (0.25 = 16th-note precision)


for c = 1:nCond
    
    condition = conditions{c};
    
    load(['M_' condition '_answer_key.mat']);
    beat_starts = correct_times(1):n_samples_in_beat:correct_times(end); nBeats = length(beat_starts);
    prop_corr = zeros(nSubs,length(beat_starts),2);

    for s = 1:length(subjects)
        
        subject = subjects(s);
        
        for run = 1:nRuns
            
            if all_completed_runs(c,s,run) == 1
                
                filename = ['s' num2str(subject) '_' condition '_' num2str(run) '.csv'];
                
                [data_tbl] = ep_clean_MIDIcsv(['csv_files/' filename]);
                
                %Extract column vectors
                last_time_slot = find(~isnan(data_tbl{:,5}), 1, 'last' ); %(cut off last few meaningless rows)
                times = data_tbl{1:last_time_slot,2}; %vector of event times
                MIDInums = data_tbl{1:last_time_slot,5}; %vector of MIDI numbers
                on_offs = data_tbl{1:last_time_slot,6}; on_offs = on_offs>1; %vector of on/offs
                
                
                %In each of b beats, check what prop. of events in the answer key also occurred in the subject
                for b = 1:nBeats
                    
                    if b == nBeats
                        curr_rows_answer = correct_times >= beat_starts(b);
                        curr_rows_subject = times >= beat_starts(b);
                    else
                        curr_rows_answer = correct_times >= beat_starts(b) & correct_times < beat_starts(b+1);
                        curr_rows_subject = times >= (beat_starts(b)-(beat_padding*n_samples_in_beat)) & times < (beat_starts(b+1)+(beat_padding*n_samples_in_beat));
                    end
                    
                    curr_data_answer = horzcat(correct_MIDInums(curr_rows_answer), correct_on_offs(curr_rows_answer));
                    curr_data_subject = horzcat(MIDInums(curr_rows_subject), on_offs(curr_rows_subject));
                    
                    %Which rows in "answer" are also in "subject"?
                    indices = find(ismember(curr_data_answer, curr_data_subject, 'rows'));
                    prop_corr(s,b,run) = length(indices)/size(curr_data_answer,1);
                end
                
                figsize = [100 100 200*length(subjects) 250];
                figure((run-1)*nCond+c); set(gcf, 'Position', figsize);
                subplot(1,nSubs,s); plot(1:nBeats, prop_corr(s,:,run), 'o', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 .4 1]); ylim([0 1]); xlim([1 nBeats]); title(['Subject ' num2str(subjects(s))]); set(gca, 'FontSize', 16); ylabel('Prop. Correct Events'); xlabel('Beat');
                print(gcf, '-dtiff', ['figures/Acc by beat_' condition '_run' num2str(run) '_beatpadding=' num2str(beat_padding) '.tif']);
            
            else
                prop_corr(s,:,run) = NaN;                  
            end
            
        end
    end
    
    prop_corr_all_conditions{c} = prop_corr;
end

save(['analyzed/prop_corr_all_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'prop_corr_all_conditions');


%Analyze group data
load(['analyzed/prop_corr_all_conditions_beatpadding=' num2str(beat_padding) '.mat']);

%Reshape the cell into an s x c x run matrix of average beat values. 
for c = 1:nCond
    data = prop_corr_all_conditions{c}; 
    for run = 1:nRuns
        data_avgbeat(:,c,run) = nanmean(data(:,:,run),2);
        nSubs_by_cond(:,run) = sum(all_completed_runs(:,:,run),2);
        
        y_compareCond(c,run) = nanmean(data_avgbeat(:,c,run),1);
        errors_compareCond(c,run) = nanstd(data_avgbeat(:,c,run),1)/sqrt(nSubs_by_cond(c,run));
        
    end
end

y_compareSubs = nanmean(data_avgbeat(:,:,:),2);
errors_compareSubs = nanstd(data_avgbeat(:,:,:),[],2)/sqrt(nCond);

%Analyze the group across conditions
for run = 1:nRuns
    figsize = [100 100 450 375]; barwidth = .6; barcolor = [.9 0 .3];
    figure('Units', 'pixels', 'Position', figsize);
    
    x = 1:nCond;
    bar(x,y_compareCond(:,run),barwidth,'facecolor', barcolor); hold on;
    errorbar(x,y_compareCond(:,run),errors_compareCond(:,run),'k.', 'LineWidth', 1)
    xlabel('Scramble condition'); ylabel('Mean Accuracy Across Beats'); ylim([0 1]); title(['Run ' num2str(run)]); set(gca, 'XTickLabel', conditions, 'FontSize', 16, 'FontName', 'Helvetica');
    print(gcf, '-dtiff', ['figures/Acc (avgd across beats), per condition, run ' num2str(run) '_beatpadding=' num2str(beat_padding) '.tif']);  
end

%Analyze data (collapsed across conditions) across subjects
for run = 1:nRuns
    figsize = [100 100 700 375]; barwidth = .6; barcolor = [.2 .4 .9];
    figure('Units', 'pixels', 'Position', figsize);
    
    x = 1:nSubs;
    bar(x,y_compareSubs(:,run),barwidth,'facecolor', barcolor); hold on;
    errorbar(x,y_compareSubs(:,run),errors_compareSubs(:,run),'k.', 'LineWidth', 1)
    xlabel('Subject'); ylabel('Mean Accuracy Across Beats'); ylim([0 1]); title(['Run ' num2str(run)]); set(gca, 'FontSize', 16, 'FontName', 'Helvetica');
    print(gcf, '-dtiff', ['figures/Acc (avgd across beats), per subject, run ' num2str(run) '_beatpadding=' num2str(beat_padding) '.tif']);
end


% all_data = vertcat(mean_acc_across_beats_I, mean_acc_across_beats_10perc, mean_acc_across_beats_20perc, mean_acc_across_beats_8B, mean_acc_across_beats_2B, mean_acc_across_beats_1B);
% x = 1:nSubs;
% y = mean(all_data);
% errors = std(all_data)/sqrt(size(all_data,1));
% 
% figsize = [100 100 325 375]; barwidth = .6; barcolor = [.2 .4 .9];
% figure('Units', 'pixels', 'Position', figsize);
% bar(x,y,barwidth,'facecolor', barcolor); hold on;
% errorbar(x,y,errors,'k.', 'LineWidth', 1)
% xlabel('Subject'); ylabel('Mean Accuracy Across Beats'); ylim([0 1]); set(gca, 'XTickLabel', num2cell(subjects), 'FontSize', 16, 'FontName', 'Helvetica');
% print(gcf, '-dtiff', ['figures/' session '/run' num2str(run) '/Acc (avgd across beats) per subject, beatpadding=' num2str(beat_padding) '.tif']);



%Note: Intact is only 151 beats because the last measure is held and the
%last beat isn't counted as a new event

%Other code
% headers = {'Track', 'TimeStamp', 'OnOff', 'Unknown', 'MIDINum', 'OnOff'};

%     c = intersect(curr_data_subject, curr_data_answer, 'rows');
%     indices = find(ismember(curr_data_subject, c, 'rows'));
%     beat_starts = times(1):n_samples_in_beat:times(end); nBeats = length(beat_starts);

% %Save the ideal data
% [data_tbl] = ep_clean_MIDIcsv(['csv_files/Pretest_8B.csv']);
% last_time_slot = find(~isnan(data_tbl{:,5}), 1, 'last' ); %(cut off last few meaningless rows)
% times = data_tbl{1:last_time_slot,2}; %vector of event times
% MIDInums = data_tbl{1:last_time_slot,5}; %vector of MIDI numbers
% on_offs = data_tbl{1:last_time_slot,6}; on_offs = on_offs>1; %vector of on/offs
%
% correct_times = times; correct_MIDInums = MIDInums; correct_on_offs = on_offs; correct_data_tbl = data_tbl;
% save('Pretest_8B_answer_key.mat', 'correct_times', 'correct_MIDInums', 'correct_on_offs', 'correct_data_tbl');
%


% Alternative way:
% -Interpolate the matrices so they have the same # of samples, then correlate/RSA each subject?s with the original
% possible_MIDI_nums = 41:77;
% note_labels = {'F2', 'F#2', 'G2', 'G#2', 'A2', 'A#2', 'B2', 'C3', 'C#3', 'D3', 'D#3', 'E3', 'F3', 'F#3', 'G3', 'G#3', 'A3', 'A#3', 'B3', ...
%     'C4', 'C#4', 'D4', 'D#4', 'E4', 'F4', 'F#4', 'G4', 'G#4', 'A4', 'A#4', 'B4', 'C', 'C#', 'D', 'D#', 'E', 'F5'};

%
% % MIDInum_range = max(MIDInums) - min(MIDInums)
% MIDI_to_row_conv = min(possible_MIDI_nums)- 1;
% % plot(1:length(times)-1,diff(times),'ro');
%
%
% %Save notes x timestamps (for both ons and offs) in 2 forms (table, vector form)
% MIDI_table = zeros(length(possible_MIDI_nums), length(times));
% ons = NaN(1,length(times)); offs = NaN(1,length(times));
%
% for i = 1:length(times)
%     if on_offs(i) == 1
%         MIDI_table(MIDInums(i)-MIDI_to_row_conv,i) = 1;
%         ons(i) = MIDInums(i)-MIDI_to_row_conv;
%     elseif on_offs(i) == 0
%         MIDI_table(MIDInums(i)-MIDI_to_row_conv,i) = -1;
%         offs(i) = MIDInums(i)-MIDI_to_row_conv;
%     end
% end
%
% %Plot
% data = round(2*(rand(20)-0.5));
% figure; hAxes = gca;
% imagesc(hAxes, flipud(MIDI_table)); set(gca, 'FontSize', 16, 'ytick', 1:length(note_labels), 'yticklabels', flip(note_labels)); colormap(hAxes, [0 1 0; 1 1 1; 0 0 1])
%
% figure; plot(1:length(times),ons, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 .3 1]); hold on; plot(1:length(times),offs,'o', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); ylim([1 length(note_labels)]); set(gca, 'FontSize', 16, 'ytick', 1:length(note_labels), 'yticklabels', flip(note_labels));
%


% if strcmp(condition, 'I')
%     mean_acc_across_beats_I = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_I', '-append');
%     
% elseif strcmp(condition, '.05')
%     mean_acc_across_beats_5perc = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_5perc', '-append');
%     
% elseif strcmp(condition, '.1')
%     mean_acc_across_beats_10perc = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_10perc', '-append');
%     
% elseif strcmp(condition, '.2')
%     mean_acc_across_beats_20perc = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_20perc', '-append');
%     
% elseif strcmp(condition, '8B')
%     mean_acc_across_beats_8B = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_8B', '-append');
%     
% elseif strcmp(condition, '2B')
%     mean_acc_across_beats_2B = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_2B', '-append');
%     
% elseif strcmp(condition, '1B')
%     mean_acc_across_beats_1B = nanmean(prop_corr_all_subs);
%     save(['analyzed/' session '/run' num2str(run) '/mean_acc_across_conditions_beatpadding=' num2str(beat_padding) '.mat'], 'mean_acc_across_beats_1B', '-append');
% end

% nSubs_by_cond = horzcat(sum(all_completed_runs(:,:,1),2),sum(all_completed_runs(:,:,2),2)); %nSubs completed each condition for each run

%Run data including .05 condition
completed_run1 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 0 1 1 1 1 1 1 1 1 1 1];

completed_run2 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 0 1 1 1 1 1 1 1 1 1 1];

completed_run3 = [
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 1 1 1 1 1 1 1 1
    1 0 0 0 0 0 0 0 0 0 0 0 0];
