function []=savesignature(sigfile,outfile)
  disp(sigfile);
  disp(outfile);
  load(sigfile);
  type = input.types;
  subtype = input.subtypes;
  totalProcesses = size(processes, 2);
  totalMutationTypes = size(processes, 1);
  [sortedType sortedIndex] = sort(type);
  sortedSubType = subtype(sortedIndex);

  for i = 1 : totalProcesses
    processes(:, i) = processes(sortedIndex, i);
  end
  
  filecon=fopen(outfile,'wt');
  for i=1:totalMutationTypes
    for j=1:totalProcesses
      fprintf(filecon,'%f\t',processes(i,j));
	end
	fprintf(filecon,'\n');
  end
  fclose(filecon);
end