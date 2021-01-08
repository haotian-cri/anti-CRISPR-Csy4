function potential_sensors = A00_Window_Sliding(ROI)
length_ROI = length(ROI);
potential_sensors = cell(length_ROI-17,1);
rcROI = upper(seqrcomplement(ROI));
for i = 1:(length_ROI-17)
    potential_sensors(i) = {[rcROI(i:i+8),'CTGCCGTATAGGCAG',rcROI(i+9:i+17)]};
end
end