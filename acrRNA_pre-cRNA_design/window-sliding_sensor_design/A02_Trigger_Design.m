function [Trigger] = A02_Trigger_Design(Sensors,terminator,nupackhome)

seqrc = cellfun(@(x) seqrcomplement(x([1:9,end-8:end])),Sensors,'UniformOutput',false);

Stem_5 =    'GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';
structure = '...((((((((((......))))))))))..............................(.((((((((((((((......)))))))))))))).)...........';

N = numel(seqrc);
Trigger = cell(N,1);

for i = 1:N
    constraint_seq = [Stem_5,dna2rna(seqrc{i}),terminator];
    str = [structure,'\n',constraint_seq];
    str = compose(str);
    dlmwrite('trigger.fold',str,'delimiter','');
    command = [nupackhome 'complexdesign -T 37 -prevent trigger.prevent trigger'];
    [~] = CommandDiary('screen.summary',command);
    if isfile('trigger.summary')
    x = importdata('trigger.summary');
    Trigger(i) = {rna2dna(x{3})};
    delete trigger.summary
    else
        warning('Terminator is not compatiable to sensor %s',i)
    end
    delete trigger.fold
    delete screen.summary
end
end