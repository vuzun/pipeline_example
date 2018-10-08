from ruffus import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.Sra as Sra

PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


#v23 transcriptome version
@transform("gencode.v23.transcripts.fa.gz",
	   suffix(".fa.gz"),
	   ".idx")
def buildKallistoIndex(infile, outfile):

    job_memory="15G"
    statement = '''kallisto index -i %(outfile)s %(infile)s > kallisto_index.log'''
    P.run()

#fetching BAMs via .remote
@follows(mkdir("kalisto_quant"))
@transform(["input_files.dir/*.remote", "input_files.dir/*.bam"],
           regex(".+/(.+).(?:remote|bam)"),
           add_inputs(buildKallistoIndex),
           r"kalisto_quant/\1/abundance.tsv.gz")
def runKalistoOnRemoteBAM(infiles, outfile):
    '''running kalisto on .bam or .remote files'''

    job_memory="6G"
    infile, kallisto_index = infiles

    outdir = os.path.dirname(outfile)
    outfile = P.snip(outfile, ".gz")

    statement = []
    

    tempfastq1 = P.getTempFilename()
    tempfastq2 = P.getTempFilename()
    rm_files = [tempfastq1, tempfastq2]
    
    
    
    if infile.endswith(".remote"):
            tempbam = P.getTempFilename()
	    s, infiles = Sra.process_remote_BAM(infile, outdir=tempbam)
    	    infile = " ".join(infiles)
            statement.append(s)
            rm_files.append(tempbam)


    statement.append('''samtools fastq %(infile)s 
                                      -1 %(tempfastq1)s 
                                      -2 %(tempfastq2)s''')

    statement.append('''kallisto quant -i %(kallisto_index)s
                                       -o %(outdir)s
                                       %(tempfastq1)s %(tempfastq2)s''')
    rm_files = " ".join(rm_files)
    statement.append('''rm -R %(rm_files)s''')
    statement.append('''gzip %(outfile)s''')
    statement = "; \n checkpoint;\n".join(statement)

    P.run()
    

@merge(runKalistoOnRemoteBAM, "counts_table.tsv.gz")
def mergeCounts(infiles, outfile):

	statement = '''Rscript test.R --stdout %(outfile)s '''
	"""
	statement = '''python %(scriptsdir)s/combine_tables.py 
                       -c 1
                       -k 5
                       --regex-filename 'kalisto_quant/(.+)/abundance.tsv.gz'
                       --use-file-prefix
                       --stdin %(infiles)s
                       --stdout %(outfile)s

#                       --stdout %(outfile)s.log'''
	""" 
	P.run()
#	infiles = " ".join(infiles)
#	"""
#	statement = '''python %(scriptsdir)s/combine_tables.py 
#                       -c 1
#                       -k 5
#                       --regex-filename 'kalisto_quant/(.+)/abundance.tsv.gz'
##                       --use-file-prefix
#                       --stdin %(infiles)s
#                       --stdout %(outfile)s'''
#                       --stdout %(outfile)s.log'''
#	""" 
#	P.run()


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
