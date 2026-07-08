function [syncMethod] = identifySyncMethod


answer = questdlg('Which method do you want to use to sync the data?', ...
	'Sync Method', ...
	'nVista SYNC LED OR AnyMaze','Random LED pulses', 'UCLA timestamps', 'Maybe');

% Handle response
switch answer
    case 'nVista SYNC LED OR AnyMaze'
        disp(['We will apply a mask on the tracking data, so the beginning of the tracking and imaging data is the same.'])
        syncMethod = 1;
    case 'Random LED pulses'
        disp([answer 'We will compute the aligner and convert the imaging timestamps'])
        syncMethod = 2;
    case 'UCLA timestamps'
        syncMethod = 3
end


