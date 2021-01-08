function [structure,energy] = A03_Binding_Complex_Analysis(Input1,Input2,nupackhome,i)
    str = ['2','\n',Input1,'\n',Input2,'\n','1 2'];
    str = compose(str);
    dlmwrite(['pair',num2str(i),'.in'],str,'delimiter','');
    command = [nupackhome 'pfunc -multi pair',num2str(i)];
    filename = ['pair',num2str(i),'.pfunc'];
    partition = CommandDiary(filename,command);
    energy = str2double(partition{3});
    delete(['pair',num2str(i),'.pfunc'])
    
    command = [nupackhome 'mfe -multi pair',num2str(i)];
    system(command);
    filename = ['pair',num2str(i),'.mfe'];
    if isfile(filename)
        fileID = fopen(filename);
        structure = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',15);
        structure = cell2mat(structure{1});
        fclose(fileID);
        delete(['pair',num2str(i),'.mfe'])
    else
        structure = {};
    end
    delete(['pair',num2str(i),'.in'])
    
end