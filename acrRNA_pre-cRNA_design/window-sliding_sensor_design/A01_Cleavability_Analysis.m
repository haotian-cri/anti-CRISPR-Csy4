function [Prediction,scores,stdevs] = A01_Cleavability_Analysis(Mdl,SeqAll,nupackhome,Viennahome)
    % count the library size
    smplsize = numel(SeqAll);
    SeqRNA = dna2rna(SeqAll);

    %% 1.Nupack analysis
    ComplexEnergy = zeros(smplsize,1);
    Pfunc = zeros(smplsize,1);
    Mfe = zeros(smplsize,1);
    ComplexCount = zeros(smplsize,1);
    MfeProb = zeros(smplsize,1);
    PPair = cell(smplsize,1);

    parfor i = 1:smplsize
        % create the input file for nupack analysis
        dlmwrite(['target',num2str(i),'.in'],SeqRNA(i,:),'delimiter','');

        % 1. calculate the partition function and complex energy
        command = [nupackhome 'pfunc target' num2str(i)];
        filename = ['target',num2str(i),'.pfunc'];
        x = CommandDiary(filename,command);
        ComplexEnergy(i) = str2double(x{3}); % free energy (kcal/mol)
        Pfunc(i) = str2double(x{4});

        % 2. calculate the base-pairing observables
        command = [nupackhome 'pairs target' num2str(i) ' -cutoff 0'];
        system(command);
        x = importdata(['target',num2str(i),'.ppairs']);
        pairs = x.data;% first line of the data is the number of bases
        pairs(1) = [];pairs = reshape(pairs,3,[]);% reshape the data into x,y,probability
        PPair(i) = {pairs(3,:)};

        % 3. count the number of secondary structures in the ensemble
        command = [nupackhome 'count target' num2str(i)];
        filename = ['target',num2str(i),'.count'];
        x = CommandDiary(filename,command);
        ComplexCount(i) = str2double(x{3});

        % 4. find the minimum free energy (MFE) secondary structure(s)
        command = [nupackhome 'mfe target' num2str(i)];
        system(command);
        x = importdata(['target',num2str(i),'.mfe']);
        Mfe(i) = str2double(x{4});

        % create the input file to analyse the mfe structure
        dlmwrite(['target',num2str(i),'.in'],x{5},'-append','delimiter','');
        % 5. calculate the equilibrium probability of a secondary structure
        command = [nupackhome 'prob target' num2str(i)];
        filename = ['target',num2str(i),'.prob'];
        x = CommandDiary(filename,command);
        MfeProb(i) = str2double(x{3});
        % delete the target files and prepare for the next round calculation
        delete(['target' num2str(i) '.mfe'])
        delete(['target' num2str(i) '.prob'])
        delete(['target' num2str(i) '.count'])
        delete(['target' num2str(i) '.pfunc'])
        delete(['target' num2str(i) '.in'])
        delete(['target' num2str(i) '.ppairs'])
    end

    %% 2. kinwalker − predicts RNA folding trajectories
    intmdt_No = zeros(smplsize,1);
    barrier_mean = zeros(smplsize,1);
    barrier_median = zeros(smplsize,1);
    barrier_std = zeros(smplsize,1);
    barrier_max = zeros(smplsize,1);
    parfor i = 1:smplsize
        % create the input file for kinwalker
        dlmwrite(['target',num2str(i),'.seq'],SeqRNA(i,:),'delimiter','');

        command = [Viennahome 'kinwalker --transcribed=3',... % arbitrary setting 
            ' --transcription_rate=55',... % nt/s, Ref: U Vogel and K F Jensen 1994, LB media
            ' <target',num2str(i),'.seq > target',num2str(i),'.kw'];

        system(command);

        fid = fopen(['target',num2str(i),'.kw']);
        x = textscan(fid, '%s');x = x{1};% extract all the data from .kw file
        fclose(fid);

        x(end-3:end) = [];% Kinwalker run time: x seconds
        x(1:4) = [];% sequence mfe structre mfe TRAJECTORY
        x = reshape(x,6,[])';
        intmdt = x(:,1); % all intermediate structrues
        data = x(:,2:6); % energy, time, barrier,threshold,length
        data = cell2mat(cellfun(@(x) str2double(x),data,'UniformOutput',false));
        intmdt_No(i) = numel(intmdt);% number of intermediate structures

        % statiss of barrier energy
        barrier_mean(i) = mean(data(:,3));
        barrier_median(i) = median(data(:,3));
        barrier_std(i) = std(data(:,3));
        barrier_max(i) = max(data(:,3));
        delete(['target',num2str(i),'.kw'])
        delete(['target',num2str(i),'.seq'])
    end
    %% 3. Sequence bioinformas
    % remove the sequence that is invariant, and extract only the 9 bp upstream
    % and downstream flanking regions which is randomized in the library
    rSeq = cellfun(@(x) x([1:9,25:end]),SeqAll,'UniformOutput',false);
    % count nucleotides in sequences
    base = cellfun(@(x) basecount(x),rSeq,'UniformOutput',false);
        % turn the data into a better format.
        BaseFreq = table(cellfun(@(x) x.A,base),cellfun(@(x) x.C,base),cellfun(@(x) x.G,base),cellfun(@(x) x.T,base));
        BaseFreq.Properties.VariableNames = {'A','C','G','T'};
        clearvars base
    % generate the integer representation of random flanking regions, turn ACGT to 1234
    rSeq_int = double(cell2mat(cellfun(@nt2int,rSeq,'UniformOutput',false)));
    % generate the prime representation of random flanking regions
    p = [2,3,5,7]; q = [11,13,17,19];
    rSeq_p = zeros(size(rSeq_int));rSeq_q = zeros(size(rSeq_int));
    for i = 1:4
        rSeq_p(rSeq_int == i) = p(i);
        rSeq_q(rSeq_int == i) = q(i);
    end
    rSeq_p = mat2cell(rSeq_p,ones(numel(rSeq),1),18);
    rSeq_q = mat2cell(rSeq_q,ones(numel(rSeq),1),18);
    % generate the pairing combination
    % create a lower triangular matrix, use 2 primer products to represent a
    % unique pairing 
    rPair_int = cellfun(@(x,y) tril(x'*y,-1),rSeq_p,rSeq_q,'UniformOutput',false);
    % turn the pairing composite number representation from matrix into vector,
    % in an order of:
    % (1,2),(1,3),...,(1,18),(2,1),(2,2),...,(16,17),(16,18),(17,18).
    % which also matches the order of ppair file from nupack calculation.
    rPair_int = cellfun(@(x) x(tril(true(size(x)),-1)),rPair_int,'UniformOutput',false);
    % count base pairing in sequences
    rPairFreq = cellfun(@(x) hist(unique(x),unique(p'*q)),rPair_int,'UniformOutput',false);
    % p = [2,3,5,7]; q = [11,13,17,19];
    rPairFreq = array2table(cell2mat(rPairFreq),...
        'VariableNames',{'pairAA','pairAC','pairCA','pairAG',...
        'pairAT','pairCC','pairCG','pairGA',...
        'pairCT','pairGC','pairTA','pairGG',...
        'pairTC','pairGT','pairTG','pairTT'});
    clearvars rSeq_p rSeq_q
    
    % count n-mer numbers in sequences
    mer2 = {'TT';'TG';'GT';'AT';'CA';'TA';'CC';'GG';'TC';'AC';'CT';'AA';'CG';'GC';'GA';'AG'};
    mer3 = {'TTT';'TGT';'TTG';'TGG';'TAT';'ATT';'GTT';'GTG';'CCC';'TAC';'TTA';'GGG';'TTC';'TCA';'ACA';'CAA';'CAT';'ATA';'CTT';'GGT';'CAC';'ACT';'AAT';'CGT';'CCA';'TCT';'AAA';'ATC';'TCC';'ACC';'CCT';'CTA';'GTC';'AAC';'GGC';'ATG';'GCG';'TGC';'GAT';'CTC';'TAA';'CTG';'TCG';'CCG';'GCC';'CGG';'CAG';'GCA';'GTA';'AGT';'GCT';'GGA';'CGC';'ACG';'AGC';'CGA';'TGA';'GAC';'AAG';'GAA';'TAG';'GAG';'AGG';'AGA'};
    mer4 = {'TTTT';'GTGT';'TTTA';'TTGT';'TGTT';'TGTG';'GTTT';'CCCC';'TTGG';'TGGG';'TGGT';'GTGG';'ATTG';'ATTT';'GTTG';'TATT';'TTAC';'TATA';'TTAT';'GATT';'CTAT';'TTTG';'TACA';'CAAT';'TCTT';'TGTC';'CCCA';'CTTT';'ACAA';'TTCA';'CACA';'ATGG';'TTTC';'CGTT';'GGTG';'AATA';'TCGT';'GGGT';'TACT';'ATTC';'TCAT';'TCCC';'ACAT';'CATC';'ATAC';'TTGC';'TCAA';'GTTA';'GTAT';'ATAT';'TGGC';'CATT';'TCAC';'AAAC';'ATCC';'CTGT';'TATG';'TACC';'ACTA';'AATT';'GGGG';'TATC';'AACA';'GCGG';'GGCG';'ACAC';'GGTT';'GGGC';'TAAA';'CTTG';'GGCC';'GTAC';'ACTT';'TTCC';'TTCT';'ATGT';'CTCA';'ATAA';'TGTA';'CCAC';'AAAA';'AGTT';'GCGT';'ATCA';'ACCA';'CACT';'TCCT';'TAAT';'CAAC';'AGTG';'TTAA';'CCCG';'CTTC';'GGGA';'GCCC';'CCAT';'ACCC';'CATG';'TCTA';'GGTC';'CAAA';'TGCA';'CGTG';'GTCC';'TACG';'CTAC';'TGGA';'AGCA';'TTGA';'ACCT';'CCGC';'ATCT';'ATTA';'CAAG';'CGAT';'TTCG';'CCTT';'GTCT';'CCTC';'GTCA';'AACT';'CAGT';'CGGG';'TCGG';'GTGC';'CCGA';'AACC';'CCAA';'CGTC';'CACC';'TGAT';'ACTC';'CCCT';'CAGC';'CCTA';'TGCT';'CATA';'GGAT';'CCTG';'CGCG';'CGGT';'ACGA';'CCAG';'CGAC';'CTCG';'CTCT';'GCTG';'GTTC';'AGGG';'GCAA';'AATC';'CTGG';'GGAC';'TCTC';'TAAC';'TGCG';'AAAG';'GCTT';'ACGT';'ACTG';'GCAT';'GAAT';'TCCA';'AAGC';'CACG';'CGGC';'GTCG';'GCAG';'GATG';'GGTA';'TTAG';'GACC';'GCCG';'TCAG';'CGCT';'GCCT';'AAAT';'TGCC';'ATAG';'TAGT';'CCGT';'GGCT';'CTTA';'GACT';'ACCG';'CAGA';'CGTA';'CTAA';'ACAG';'TGAA';'CGCC';'GCTC';'AGCG';'TCCG';'GCGC';'ATGC';'GCAC';'GCCA';'TCTG';'AAGT';'GAGT';'ACGC';'AGAT';'AATG';'GAAA';'CTGC';'AGCC';'GTGA';'AAGG';'GCGA';'GGAG';'ATCG';'TCGC';'AACG';'CGGA';'GGAA';'TGAG';'AGAG';'GGCA';'CCGG';'GAGC';'GATC';'GCTA';'CTGA';'TAGC';'CGCA';'CAGG';'AGCT';'CGAA';'GATA';'GAGG';'GACA';'GAAC';'TAGA';'CTCC';'ACGG';'AGAA';'GACG';'CTAG';'AGTA';'GTAA';'ATGA';'AGGC';'TGAC';'AGGT';'GAGA';'CGAG';'AAGA';'TAGG';'GAAG';'TCGA';'AGAC';'AGGA';'AGTC';'TAAG';'GTAG'};
    % count for each possible n-mer
    Mer2 = array2table(cell2mat(cellfun(@(x) count(rSeq,x)',mer2,'UniformOutput',false))',...
        'VariableNames',mer2');
    Mer3 = array2table(cell2mat(cellfun(@(x) count(rSeq,x)',mer3,'UniformOutput',false))',...
        'VariableNames',mer3');
    Mer4 = array2table(cell2mat(cellfun(@(x) count(rSeq,x)',mer4,'UniformOutput',false))',...
        'VariableNames',mer4');
    clearvars mer2 mer3 mer4
    
    
    % translation outcomes  
    SeqAA = nt2aa(SeqAll);
    % position of the early stops
    posStop = zeros(numel(SeqAA),1);
    posStop2 = zeros(numel(SeqAA),1);
    for i = 1
    % strfind outputs the location of stop codon
    % numel()==i determine if there are i stop codons
    % find output all the seq having i stop codons
    index = find(cellfun(@(x) numel(x)==i,strfind(SeqAA,'*')));
    posStop(index,1:i) = cell2mat(strfind(SeqAA(index),'*'));
    end
    for i = 1
        index = find(cellfun(@(x) numel(x)==i,strfind(SeqAA,'**')));
        posStop2(index,1:i) = cell2mat(strfind(SeqAA(index),'**'));
    end

    % Codon usage
    % codon usage reference: https://www.genscript.com/tools/codon-frequency-table
    CodonUsage = {"TTT",0.580000000000000,22.1000000000000;"TTC",0.420000000000000,16;"TTA",0.140000000000000,14.3000000000000;"TTG",0.130000000000000,13;"TAT",0.590000000000000,17.5000000000000;"TAC",0.410000000000000,12.2000000000000;"TAA",0.610000000000000,2;"TAG",0.0900000000000000,0.300000000000000;"CTT",0.120000000000000,11.9000000000000;"CTC",0.100000000000000,10.2000000000000;"CTA",0.0400000000000000,4.20000000000000;"CTG",0.470000000000000,48.4000000000000;"CAT",0.570000000000000,12.5000000000000;"CAC",0.430000000000000,9.30000000000000;"CAA",0.340000000000000,14.6000000000000;"CAG",0.660000000000000,28.4000000000000;"ATT",0.490000000000000,29.8000000000000;"ATC",0.390000000000000,23.7000000000000;"ATA",0.110000000000000,6.80000000000000;"ATG",1,26.4000000000000;"AAT",0.490000000000000,20.6000000000000;"AAC",0.510000000000000,21.4000000000000;"AAA",0.740000000000000,35.3000000000000;"AAG",0.260000000000000,12.4000000000000;"GTT",0.280000000000000,19.8000000000000;"GTC",0.200000000000000,14.3000000000000;"GTA",0.170000000000000,11.6000000000000;"GTG",0.350000000000000,24.4000000000000;"GAT",0.630000000000000,32.7000000000000;"GAC",0.370000000000000,19.2000000000000;"GAA",0.680000000000000,39.1000000000000;"GAG",0.320000000000000,18.7000000000000;"TCT",0.170000000000000,10.4000000000000;"TCC",0.150000000000000,9.10000000000000;"TCA",0.140000000000000,8.90000000000000;"TCG",0.140000000000000,8.50000000000000;"TGT",0.460000000000000,5.20000000000000;"TGC",0.540000000000000,6.10000000000000;"TGA",0.300000000000000,1;"TGG",1,13.9000000000000;"CCT",0.180000000000000,7.50000000000000;"CCC",0.130000000000000,5.40000000000000;"CCA",0.200000000000000,8.60000000000000;"CCG",0.490000000000000,20.9000000000000;"CGT",0.360000000000000,20;"CGC",0.360000000000000,19.7000000000000;"CGA",0.0700000000000000,3.80000000000000;"CGG",0.110000000000000,5.90000000000000;"ACT",0.190000000000000,10.3000000000000;"ACC",0.400000000000000,22;"ACA",0.170000000000000,9.30000000000000;"ACG",0.250000000000000,13.7000000000000;"AGT",0.160000000000000,9.90000000000000;"AGC",0.250000000000000,15.2000000000000;"AGA",0.0700000000000000,3.60000000000000;"AGG",0.0400000000000000,2.10000000000000;"GCT",0.180000000000000,17.1000000000000;"GCC",0.260000000000000,24.2000000000000;"GCA",0.230000000000000,21.2000000000000;"GCG",0.330000000000000,30.1000000000000;"GGT",0.350000000000000,25.5000000000000;"GGC",0.370000000000000,27.1000000000000;"GGA",0.130000000000000,9.50000000000000;"GGG",0.150000000000000,11.3000000000000};
    CodonFraction = codon_fraction(SeqAll,CodonUsage);
    CodonFreq1k = codon_freq(SeqAll,CodonUsage);
    CodonFraction_mean = mean(CodonFraction,2);CodonFraction_var = std(CodonFraction,[],2);
    CodonFreq1k_mean = mean(CodonFreq1k,2);CodonFreq1k_var = std(CodonFreq1k,[],2);   
    %% 4. Output features
    % data clearning
    % remove some of the dataset without uniformed format

    % reshape the format for some of the cell array dataset
    PPair = cell2mat(PPair);PPair = array2table(PPair);
    rPair_int = cell2mat(rPair_int')';rPair_int = array2table(rPair_int);
    CodonFraction = array2table(CodonFraction);
    CodonFreq1k = array2table(CodonFreq1k);
    rSeq_int = array2table(rSeq_int);
    clearvars rSeq SeqAll % rSeq_int as the alternative
    
    % Convert aa translation from letter to integer representation
    SeqAA = cell2mat(cellfun(@(x) double(aa2int(x)),SeqAA,'UniformOutput',false));
    SeqAA = array2table(SeqAA);
    % create the feature table.
    fRawdata = [rSeq_int,SeqAA];
    fEnzemble = table(Pfunc,ComplexCount,ComplexEnergy);
    fKinwalker = table(barrier_max,barrier_mean,barrier_median,barrier_std,intmdt_No);                     
    fMfe = table(Mfe,MfeProb);
    fPairWise = [PPair,rPair_int,rPairFreq];
    fSeq = [BaseFreq,Mer2,Mer3,Mer4];
    fTranslation = [table(posStop(:,1),posStop2(:,1),CodonFraction_mean,CodonFraction_var,...
        CodonFreq1k_mean,CodonFreq1k_var),CodonFraction,CodonFreq1k];
    Features = [fRawdata,fEnzemble,fKinwalker,fMfe,fPairWise,fSeq,fTranslation];
    clearvars -except Features Mdl
    %% 5. Predict the output
    [Countfit,scores,stdevs] = predict(Mdl,Features);
    Prediction = cellfun(@str2double,Countfit);
end

function a = codon_fraction(SeqAll,CodonUsage)
    smplsize = numel(SeqAll);
    rSeq = cellfun(@(x) x([1:9,25:end]),SeqAll,'UniformOutput',false);

    a = zeros(numel(rSeq),6);
    for i = 1:6
        % find strcmp functions locate the triplet in the codon usage table
        codon = cellfun(@(x) x((3*i-2):3*i),rSeq,'UniformOutput',false);
        codon_index = cellfun(@(x) find(contains([CodonUsage{:,1}],x)),codon,'UniformOutput',false);
        for j = 1:smplsize
            a(j,i) = CodonUsage{codon_index{j},2};     
        end
    end
end

function a = codon_freq(SeqAll,CodonUsage)
    smplsize = numel(SeqAll);
    rSeq = cellfun(@(x) x([1:9,25:end]),SeqAll,'UniformOutput',false);
    a = zeros(numel(rSeq),6);
    for i = 1:6
        % find strcmp functions locate the triplet in the codon usage table
        codon = cellfun(@(x) x((3*i-2):3*i),rSeq,'UniformOutput',false);
        codon_index = cellfun(@(x) find(contains([CodonUsage{:,1}],x)),codon,'UniformOutput',false);
        for j = 1:smplsize
            a(j,i) = CodonUsage{codon_index{j},3};     
        end
    end
end
