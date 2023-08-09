#' @title Clustering of individuals
#' @param data data file
#' @param format Data format. Format supported: "FASTA", "VCF" ,"BAM", "GenePop"
#' @param verbose if \code{TRUE}, prints extra output information
#' @importFrom utils read.delim
#' @importFrom vcfR read.vcfR
#' @importFrom Rsamtools scanBam
#' @importFrom adegenet read.genepop .readExt
#' @references Samtools: a suite of programs for interacting
#' with high-throughput sequencing data. <http://www.htslib.org/>
#' @export
#' @examples
#' data <- system.file("extdata", "FASTA_clustering_haploid.fasta", package = "rBAPS")
#' greedyMix(data)
greedyMix <- function(data, format, fixedK = FALSE, partition_compare = FALSE, verbose = TRUE) {
  data <- importFile(data, format, verbose)
  if (partition_compare) {
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
    #     partitionLogml(i) = ...
    #         initialCounts(partitionInd, data(:,1:end-1), npops, c.rows, noalle, adjprior);

    # end
    # % return the logml result
    # partitionCompare.logmls = partitionLogml;
    # set(h1, 'userdata', partitionCompare);
    # return
  }

  if (fixedK) {
      # [logml, npops, partitionSummary]=indMix_fixK(c);
  } else {
      # [logml, npops, partitionSummary]=indMix(c);
  }
  stopifnot(logml != 1)

  # data = data(:,1:end-1);

  # h0 = findobj('Tag','filename1_text');  inp = get(h0,'String');
  # h0 = findobj('Tag','filename2_text');
  # outp = get(h0,'String');
  # changesInLogml = writeMixtureInfo(logml, rowsFromInd, data, adjprior, priorTerm, ...
  #     outp,inp,partitionSummary, popnames, fixedK);

  # viewMixPartition(PARTITION, popnames);

  talle <- questdlg(
    c(
      'Do you want to save the mixture populations ',
      'so that you can use them later in admixture analysis?'
    ),
    'Save results?',
    defbtn = 'n'
  );
  if (talle == "y") {
  #     [filename, pathname] = uiputfile('*.mat','Save results as');
  #     if (sum(filename)==0) || (sum(pathname)==0)
  #         % Cancel was pressed
  #         return;
  #     else
  #         % copy 'baps4_output.baps' into the text file with the same name.
  #         if exist('baps4_output.baps','file')
  #         copyfile('baps4_output.baps',[pathname filename '.txt'])
  #         delete('baps4_output.baps')
  #         end
  #     end;
  #     c.PARTITION = PARTITION; c.COUNTS = COUNTS; c.SUMCOUNTS = SUMCOUNTS;
  #     c.alleleCodes = alleleCodes; c.adjprior = adjprior; c.popnames = popnames;
  #     c.rowsFromInd = rowsFromInd; c.data = data; c.npops = npops;
  #     c.noalle = noalle; c.mixtureType = 'mix';
  #     c.logml = logml; c.changesInLogml = changesInLogml;
  #     save([pathname filename], 'c', '-v7.3'); % added by Lu Cheng, 08.06.2012
  } else {
  #     if exist('baps4_output.baps','file')
  #     delete('baps4_output.baps')
  #     end
  }
}
