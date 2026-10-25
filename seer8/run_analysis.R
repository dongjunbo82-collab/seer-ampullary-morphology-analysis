# Run only with an independently authorized local SEER export.
args <- commandArgs(TRUE)
if(length(args)<2) stop('Usage: Rscript run_analysis.R INPUT.txt OUTPUT_DIRECTORY [full|smoke] [comma-separated module numbers]')
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
code <- normalizePath(dirname(script),winslash='/',mustWork=TRUE)
input <- normalizePath(args[1],winslash='/',mustWork=TRUE)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
output <- normalizePath(args[2],winslash='/',mustWork=TRUE)
inside <- function(p) identical(tolower(p),tolower(code))||startsWith(tolower(p),paste0(tolower(code),'/'))
if(inside(input)||inside(output))stop('Keep ALL data and generated outputs outside the repository.')
mode <- if(length(args)>=3) args[3] else 'full'
stopifnot(mode %in% c('full','smoke'))
modules <- c('01_import_describe.R','02_secondary_hazards.R','03_original_overlap_risks.R',
 '04_grade_comparator_ipcw.R','05_reduced_ipcw.R','06_primary_and_bounded.R',
 '07_composition_and_subcodes.R','08_coding_selection.R','09_verify_selection.R',
 '10_figure.R','11_followup.R','12_characteristics.R')
steps <- if(length(args)>=4)as.integer(strsplit(args[4],',',fixed=TRUE)[[1]])else seq_along(modules)
stopifnot(all(steps %in% seq_along(modules)))
if(mode=='smoke')warning('SMOKE TEST ONLY: five resamples; intervals are not publication estimates.')
rs <- file.path(R.home('bin'),'Rscript')
for(i in steps){
 cat('Running module ',i,': ',modules[i],' [',mode,']\n',sep='');flush.console()
 cmd <- c(shQuote(file.path(code,'R','worker.R')),shQuote(code),shQuote(input),shQuote(output),mode,modules[i])
 status <- system2(rs,cmd)
 if(status!=0)stop('Module failed: ',modules[i],'; no subsequent module executed.')
}
cat('Completed requested modules. Outputs remain local; do not upload them.\n')
