% This program count the number of RNA and the number of
% mitochondrial-colocalized mRNA for each cell

tdfread('RNA_mito_distance.txt'); % open the RNA_mito_distance.txt file
output_file_name = sprintf ('number_of_mito_colocalizing_RNA.txt'); % outpput file name
fileID_write = fopen (output_file_name, 'w');
fprintf(fileID_write, 'cell_number\tnumber_of_RNA\tnumber_of_mito_RNA\tratio_of_mito_RNA\n'); % write headers
fclose (fileID_write);

distance_threshold = 500; %define the threshold distance between RNA and mito, I use 500 nm
number_of_cell = length (unique (cell_number));

for i = 0 : (number_of_cell-1)
    dlmwrite(output_file_name, [i, sum(cell_number==i), sum(cell_number == i & RNA0x2Dmito_distance < distance_threshold), sum(cell_number == i & RNA0x2Dmito_distance < distance_threshold) / sum(cell_number==i)], '-append','Delimiter', '\t'); %Write in the output file: cell_number, number of RNA, number of mito RNA, ratio of mito RNA
    
end
 

