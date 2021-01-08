function x = CommandDiary(filename,command)
% operate the command in the system terminal
% and save the output that was written in the screen
dlmwrite(filename,evalc('system(command);'),'delimiter','');
% importdata into a temperory array.
x = importdata(filename);