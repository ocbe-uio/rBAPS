comparePartitions <- function(data, c.rows, partitionCompare.partitions, ninds, rowsFromInd, noalle, adjprior) {
  stop("Comparing partitions not yet implemented") # TODO: implement
  # nsamplingunits = size(c.rows,1);
  # partitions = partitionCompare.partitions;
  # npartitions = size(partitions,2);
  # partitionLogml = zeros(1,npartitions);
  # for i = 1:npartitions
  #     % number of unique partition lables
  #     npops = length(unique(partitions(:,i)));

  #     partitionInd = zeros(ninds*rowsFromInd,1);
  #     partitionSample = partitions(:,i);
  #     for j = 1:nsamplingunits
  #         partitionInd([c.rows(j,1):c.rows(j,2)]) = partitionSample(j);
  #     end
  #     partitionLogml(i) = initialCounts(partitionInd, data(:,1:end-1), npops, c.rows, noalle, adjprior);

  # end
  # % return the logml result
  # partitionCompare.logmls = partitionLogml;
  # set(h1, 'userdata', partitionCompare);
  # return
}
