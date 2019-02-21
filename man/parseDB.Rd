\name{parseDB}
\alias{parseDB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse a .fasta Database
}
\description{
parseDB() parses a .fasta proteome database into a data table containing each protein accession ID in one column and the protein sequence in another column. Can also filter out reverse and contaminant sequences.
}
\usage{
parseDB(database, db_source, filt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{
a .fasta proteome database.
}
  \item{db_source}{
a string denoting the source of the input database. Key: "UP" - Uniprot; "RS" - RefSeq; "HM" - Homemade.
}
  \item{filt}{
a boolean variable (TRUE/FALSE) specifying whether or not to filter out reverse and contaminant sequences.
}
}
\details{
File extension of database should be changed from ".fasta" to ".txt" prior to import into R. This function was built to organize .fasta databases for easier manipulation in R (such as prior to input into the phindPTMs() function) or other analysis software.
}
\value{
a data table with two columns: protein accession ID and protein sequence
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
examples.path = system.file("/extdata", package = "PTMphinder")
uniprot_ref.path = paste(examples.path, "/Human_Uniprot_Example.txt", sep="")

# Read in proteome database example (or your own database)
uniprot_ref = readLines(uniprot_ref.path)

# Parse proteome database into 2 columns: protein accession and protein sequence (may take a while depending on database size)
parseDB_Example <- parseDB(uniprot_ref, "UP", FALSE)

# Create file name for parsed database
filename1 <- paste(examples.path, "/Human_Uniprot_Parsed_", Sys.Date(), ".txt", sep="")

# Write parsed  database to new file (will be found in the PTMphinder/extdata/ directory within R.framework)
write.table(parseDB_Example, filename1, row.names=FALSE, col.names=FALSE, sep="\t")
}
