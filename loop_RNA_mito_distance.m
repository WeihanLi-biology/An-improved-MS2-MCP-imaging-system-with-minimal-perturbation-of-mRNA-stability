% This program measure the distance between RNA and mito distance in
% multiple cells. Input files are re-formated FISH Quant output files and
% mitoGraph output files. 

FISH_folder ='/Users/weihanli/Dropbox (EinsteinMed)/shared files/NMD-resistant MS2 Weihan-Anna/092921_FISH_ATP2_Su9-GFP/092921_W303WT_ATP2Q570_Su9-GFP/RNA-mito distance/RNA coordinates';%the folder where RNA coordinates are
mito_folder = '/Users/weihanli/Dropbox (EinsteinMed)/shared files/NMD-resistant MS2 Weihan-Anna/092921_FISH_ATP2_Su9-GFP/092921_W303WT_ATP2Q570_Su9-GFP/RNA-mito distance/mito coordinates (use)';%the folder where mito coordinates are

FISH_files = dir(append(FISH_folder,'/FISH_*.txt')); % create a directory that contains all the FISH_*.txt files in the FISH folder
mito_files = dir(append(mito_folder,'/mito*.txt')); % create a directory that contains all the ATP2*.txt files in the mito folder

pixel_size = 64.5; % This defines the xy pixel size in nm. It is 64.5nm for 100x on DS1
xy_dim = 200; % This is the size of the crop from MitoGraph (the unit is in pixel). It is 200 by default in mitoGraph. 

number_of_file = length (FISH_files); % count the number of files

output_file_name = sprintf ('RNA_mito_distance.txt'); % creat a string to annotate the CY5 file name (to write in)
fileID_write = fopen (output_file_name, 'w');
fprintf(fileID_write, 'cell_number\tRNA-mito distance\tRNA_X\tRNA_Y\tRNA_Z\tRNA-projection_X\tRNA-projection_Y\tRNA-projection_Z\tRNA intensity (raw)\tRNA intensity (filtered)\n'); % write headers
fclose (fileID_write);

for i = 1: number_of_file
    
    %FISH_file_name = sprintf (FISH_files(i).name); % creat a string to FISH file name (to read from)
    %mito_file_name = sprintf (mito_files(i).name); % creat a string to mito file name (to read from)
    
    tdfread(append(FISH_folder,'/',FISH_files(i).name),'\t'); % read the FISH files, which is the re-formated FISH quant output files
    RNA = [Pos_X, Pos_Y, Pos_Z]; % Pos_X, Pos_Y, Pos_Z are in nm. They are the three columns from the .txt file. They document the coordinates of RNA. Combining them give the RNA coodinates. RNA is a nx3 matrix 
    
    tdfread(append(mito_folder,'/',mito_files(i).name),'\t'); % read the mito file, which is the output file from mitograph
    mito=[x,y,z]; % x,y,z are in um. They are the three columns from the .txt file. They document the coordinates of mitochondria. Combining them give the mitochondrial coodinates. mito is a nx3 matrix
    mito = mito * 1000; % converts coordinates into nm from um
    mito(:,2) = xy_dim * pixel_size - mito(:,2); % This line needs to be done because Mitograph flips the y axis

    distance = pdist2(mito, RNA); % measure the distance between every RNA and every mito skeleton dot
    [min_distance, IndexOfMinimum] = min (distance); % min_distance is the minimum distance between RNA and mito. IndexOfMinimum is the index number of the closest mito skeleton dot 
    mito_projection = mito(IndexOfMinimum,:); % mito_projection is the cloest mito dot for each RNA
   
    cell_number = zeros (length(RNA),1);
    cell_number(:) = i-1; % These two lines creat a column that shows which cell the corresponding RNA comes from
    dlmwrite(output_file_name, [cell_number, transpose(min_distance), RNA, mito_projection, INT_raw, INT_filt], '-append','Delimiter', '\t'); %Write in the output file: cell_number, RNA-mito distance, RNA_X, RNA_Y, RNA_Z, RNA-projection_X, RNA-projection_Y, RNA-projection_Z, RNA intensity (raw), RNA intensity (filtered)

end
