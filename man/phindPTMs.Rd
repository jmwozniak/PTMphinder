\name{phindPTMs}
\alias{phindPTMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Localize PTMs from Proteomic Data
}
\description{
phindPTMs() uses modified peptide information to localize the PTMs in the full-length protein and returns flanking sequences for subsequent analyses.
}
\usage{
phindPTMs(data_table, reference_table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_table}{
a data table containing your experimentally acquired proteomic data. See "phindPTMs_Input_Example.csv" for an example.
}
  \item{reference_table}{
a data table containing full protein sequence information for the peptides in the data_table. This input can be created from any proteome database using the parseDB() function. See "Human_Uniprot_Parsed_Example.txt" for an example.
}
}
\details{
input files must contain the following 6 columns (exactly as below - see "phindPTMs_Input_Example.csv"):

  1) Identifier - unique identifier for each modified peptide

  2) Protein_ID - protein accession ID

  3) Peptide_Seq - peptide sequence detected in experiment

  4) Total_Sites - total number of modified sites on peptide

  5) PTM_Loc - (potential) locations of PTM sites. Multiple PTMs should be separated by a semi-colon  (eg. "S9;S18")

  6) PTM_Score - confidence scores of PTM localizations. Multiple scores should be separated by a semi-colon  (eg. "98.56;99.84")

}
\value{
a data table with 8 columns:

  1) Identifier - see above

  2) Protein ID - see above

  3) Pep_Loc - location of the PTMs in the identified peptide

  4) Prot_Loc - location of the PTMs the full-length protein

  5) PTM_Score - see above

  6) Flank_Seq - flanking sequences extracted from PTMs

  7) Ambiguity - ambiguity of PTMs based on input data

  8) Prot_Seq - full-length protein from which motifs were extracted

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jacob M. Wozniak (jakewozniak@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Locate example files
examples.path = system.file("extdata", package = "PTMphinder")
ptm_data.path = paste(examples.path, "/phindPTMs_Input_Example.csv", sep="")
parsed_ref.path = paste(examples.path, "/Human_Uniprot_Parsed_Example.txt", sep="")

# Read in proteomic data example (or your own data)
ptm_data <- read.csv(ptm_data.path, header=TRUE, row.names = NULL, stringsAsFactors = FALSE)

# Read in parsed Uniprot database example (or one you created with parseDB above)
parsed_ref <- read.table(parsed_ref.path, header=FALSE, row.names=NULL, sep="\t")

# Locate PTMs within full-length proteins and extract neighboring motifs
phindPTMs_Example <- phindPTMs(ptm_data, parsed_ref)

# Create file name for phindPTMs output
fileName3 <- paste(examples.path, "/PTMs_phound_", Sys.Date(), ".txt", sep="")

# Write phindPTMs output to new file (will be found in the PTMphinder/extdata/ directory within R.framework)
write.table(phindPTMs_Example, fileName3, row.names=FALSE, append = FALSE,  sep = "\t")
}

