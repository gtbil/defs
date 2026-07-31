#!/usr/bin/env perl

########################################################################################################
# Copyright Notice and Disclaimer for BatchPrimer3
#
# Copyright (c) 2007, Depatrtment of Plant Sciences, University of California,
# and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#
# Authors: Frank M. You, and Dave Hane 
#
# The BatchPrimer3 web application is developed based on the Primer3 core program and Primer3Web
# of the Whitehead Institute (see Copyright notice and Disclaimer below).  This product itself must
# include the Primer3 core program developed by the Whitehead Institute for Biomedical Research.
#
# Redistribution and use in source and binary forms of the BatchPrimer3 CGI scripts, with or without
# modification,are permitted provided that the following conditions are met:
#
#   1. Redistributions must reproduce the above copyright notice, this list of conditions and the
#      following disclaimer in the documentation and/or other materials provided with the distribution.
#      Redistributions of source code must also reproduce this information in the source code itself.
#   2. If the program is modified, redistributions must include a notice (in the same places as above)
#      indicating that the redistributed program is not identical to the version distributed by Whitehead
#      Institute.
#   3. All advertising materials mentioning features or use of this software must display the following
#      acknowledgment: This product includes software developed by Department of Plant Sciences, UC Davis
#      and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#   4. The name of the UC Davis and USDA-ARS-WRRC may not be used to endorse or promote products derived
#      from this software without specific prior written permission. 
#
#
# THE BATCHPRIMER3 APPLICATION IS PROVIDED BY THE ERASMUS MC ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ERASMUS MC BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL,EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
# Copyright Notice and Disclaimer for Primer3
# Copyright (c) 1996,1997,1998 Whitehead Institute for Biomedical Research. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
#   1. Redistributions must reproduce the above copyright notice, this list of conditions and the
#      following disclaimer in the documentation and/or other materials provided with the distribution.
#      Redistributions of source code must also reproduce this information in the source code itself.
#   2. If the program is modified, redistributions must include a notice (in the same places as above)
#      indicating that the redistributed program is not identical to the version distributed by
#      Whitehead Institute.
#   3. All advertising materials mentioning features or use of this software must display the following
#      acknowledgment: This product includes software developed by the Whitehead Institute for Biomedical
#      Research.
#   4. The name of the Whitehead Institute may not be used to endorse or promote products derived from
#      this software without specific prior written permission. 

# We also request that use of this software be cited in publications as
# Steve Rozen and Helen J. Skaletsky (2000) Primer3 on the WWW for general users and for biologist
# programmers. In: Krawetz S, Misener S (eds) Bioinformatics Methods and Protocols: Methods in Molecular
# Biology. Humana Press, Totowa, NJ, pp 365-386.
# (Code available at http://www-genome.wi.mit.edu/genome_software/other/primer3.html.)
# THIS SOFTWARE IS PROVIDED BY THE WHITEHEAD INSTITUTE ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE WHITEHEAD INSTITUTE BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
########################################################################################################

# THE FOLLOWING PACKAGES MUST BE AVAILABLE IN YOUR SERVER MACHINE. PLEASE CHECK YOUR SERVER MACHINE TO
# SEE IF THEY ARE AVAILABLE. IF NOT, YOU NEED INSTALL THEN BEFOR RUN THE APPLICATION

use FileHandle;
use IPC::Open3;
use Carp;
use CGI qw(:standard);
use Socket;

use Primer3Output;
use PrimerPair;
use Primer;
use OligoTM;
use QuickSort;
use BStats;

use GD::Graph::bars;
use GD::Graph::colour;
use GD::Text;
use POSIX qw(ceil floor strftime);

use Email::Valid;
#use Thread;
  
use Fcntl qw ( :flock);

use warnings;
use strict;

########################################################################################################
# You may wish to modify the following variables to suit your installation.

# the end user will complain to:
my $MAINTAINER
    ='<a href="mailto:frankyou@pw.usda.gov">frankyou@pw.usda.gov</a>';

my $MAIL_FROM = 'frankyou@pw.usda.gov'; 
my $SENDMAIL_CMD = "|/usr/sbin/sendmail";

# The location of the primer3_core executable program. This is a platform dependant program. Users need
# to get or compile the correct version of the program
my $PRIMER_BIN =  '__BP3_CGI__/primer3_core';

# YOU MUST SET CORRECT DIRECTORIES OR URL FOR BATCHPRIMER3 SERVER. ALL DIRECTORIES MUST BE WRITABLE
#############################################################################
my $BATCHPRIMER3_CGI_HOME = "__BP3_CGI__";
my $BATCHPRIMER3_HTDOC_DIRECTORY = "__BP3_HTDOCS__";
my $BATCHPRIMER3_HTDOC_HOME = "http://127.0.0.1:__BP3_PORT__";

# The URL for help regarding this screen (which will normally
# be in the same directory as the this script)
my $HTDOC = $BATCHPRIMER3_HTDOC_HOME;
my $CGI_URL = "http://127.0.0.1:__BP3_PORT__/cgi-bin";

#Tempory result and parameter directory. chmod 777 to the directory
my $RESULT_DIR = "$BATCHPRIMER3_HTDOC_DIRECTORY/batch_primers";
my $RELATIVE_RESULT_DIR = "$HTDOC/batch_primers";
my $PARAM_DIR = "$BATCHPRIMER3_CGI_HOME/parameters";
#############################################################################

# Web application version
my $CGI_VERSION = "BatchPrimer3 v1.0";

# batchprimer3 HTML style file
my $HTML_STYLE_FILE = "style_batchprimer3.css";

# Add mispriming / mishybing libraries;  make corresponding changes
# to the $SELECT_SEQ_LIBRARY variable in batchprimer3.cgi
my %SEQ_LIBRARY=
    ('NONE' => '',
     # Put more repeat libraries here, e.g.
     'Arabidopsis' => 'repeat_libs/TIGR_Arabidopsis_Repeats.fasta',
     'Brachypodium' => 'repeat_libs/brachy_whole_genome_repeats.fasta',
     'Brassica' => 'repeat_libs/TIGR_Brassica_GSS_Repeats.fasta',
     'Glycine' => 'repeat_libs/TIGR_Glycine_Repeats.fasta',
     'Gramineae' => 'repeat_libs/TIGR_Gramineae_Repeats.v3.fasta',
     'Hordeum' => 'repeat_libs/TIGR_Hordeum_Repeats.v3.fasta',
     'Lotus' => 'repeat_libs/TIGR_Lotus_GSS_Repeats.fasta',
     'Oryza' => 'repeat_libs/TIGR_Oryza_Repeats.v2.fasta',
     'Sorghum' => 'repeat_libs/TIGR_Sorghum_GSS_Repeats.fasta',
     'Solanum' => 'repeat_libs/TIGR_Solanum_Repeats.v3.fasta',
     'Triticum' => 'repeat_libs/TIGR_Triticum_Repeats.v3.fasta',
     'Zea' => 'repeat_libs/TIGR_Zea_GSS_Repeats.fasta',
     'Wheat' => 'repeat_libs/trep.fasta'
     );


my @primer_type_descriptions = (
    "Generic primers",
    "SSR screening and primers",
    "Hybridization oligo primers",
    "SNP (Allele) flanking primers",
    "Single base extension (SBE) primers and SNP (allele) flanking primers",
    "Single base extension (SBE) primers",
    "Allele-specific primers and allele-flanking primer pairs",
    "Allele-specific primers",
    "Tetra-primer ARMS PCR primers",
    "Sequencing primers");

########################################################################################################


BEGIN{
    print "Content-type: text/html\n\n";

    # Ensure that errors will go to the web browser.
    open(STDERR, ">&STDOUT");
    $| = 1;   
    print '';
}

# remove expired directories
&remove_expired_results();

my $query = new CGI;

my $primer_type = $query->param('PRIMER_TYPE');
$primer_type = 1 if (!$primer_type);

my $email = $query->param('EMAIL');


my $tmpurl = $query->url;
my $DO_NOT_PICK;

my $ip = $ENV{'REMOTE_ADDR'};
my $user_id = $ip . "_" . time();


my $wrapup = "<hr>\n<pre>$CGI_VERSION</pre></div>".$query->end_html;
print $query->start_html(
    -title => "Report of BatchPrimer3 Primer Design ($CGI_VERSION)",
    -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
);
  
# Primer design results
# file name for table format

# make a user result directory 
my $user_result_dir = "$RESULT_DIR/$user_id";
my $dir_name = $user_id;

if (! -e $RESULT_DIR) {
    if (!mkdir($RESULT_DIR)) {
        print "<p>Can't create the result directory. Please check if your htdoc directory is correct. Probably you need to change the mode of the directory: chmod 777</p>"; 
        exit;
    }
    chmod (0777, $RESULT_DIR); 
}
if (! -e $user_result_dir) {
    if (!mkdir($user_result_dir)) {
        print "<p>Can't create the user's result directory. Please check if your htdoc directory is correct. Probably you need to change the mode of the directory: chmod 777</p>"; 
        exit;
    }
    chmod (0777, $user_result_dir); 
}

my $relative_user_result_dir = "$RELATIVE_RESULT_DIR/$user_id";


#pre-analysis of input sequences
if ($query->param('pre-analysis')) {
    &pre_analysis_of_input_sequences();
    exit;
}


# primer design 
my $PR_DEFAULT_PRODUCT_MIN_SIZE   = 100;
my $PR_DEFAULT_PRODUCT_MAX_SIZE   = 1000;
my $PRIMER_SALT_CONC              = 50.0;
my $PRIMER_DNA_CONC               = 50.0;

my $PRIMER_MAX_DIFF_TM = 5;

# for SNP primer or allele-specific primers
my $SNP_PRIMER_MIN_GC             = 40.0;
my $SNP_PRIMER_MAX_GC             = 70.0;
my $SNP_PRIMER_MIN_TM             = 50.0;
my $SNP_PRIMER_OPT_TM             = 60.0;
my $SNP_PRIMER_MAX_TM             = 65.0;
my $SNP_PRIMER_MIN_SIZE           = 14;
my $SNP_PRIMER_OPT_SIZE           = 20;
my $SNP_PRIMER_MAX_SIZE           = 30;
my $SNP_PRIMER_MAX_N              = 0;
my $SNP_PRIMER_SALT_CONC          = 50;
my $SNP_PRIMER_DNA_CONC           = 50;
my $SNP_MAX_SELF_COMPLEMENTARITY  = 8;
my $SNP_MAX_3_SELF_COMPLEMENTARITY= 3;

my $SECOND_MISMATCH = 0;
my $SECOND_MISMATCH_POS = 3;


# for tetra-primer ARMS PCR only
my $SNP_INNER_PRODUCT_MIN_SIZE = 100;
my $SNP_INNER_PRODUCT_OPT_SIZE = 200;
my $SNP_INNER_PRODUCT_MAX_SIZE = 300;
my $SNP_INNER_PRODUCT_MIN_DIFF = 1.1;
my $SNP_INNER_PRODUCT_MAX_DIFF = 1.5;

# for sequencing primer only
my $SEQUENCING_DIRECTION = "toward3";
my $PRIMER_NUM_RETURN = 1;

# for complementarity check
my $MATCH     =    1;
my $MISMATCH  =   -1;
my $GAP       =   -2;
my $N_SCORE   =   -0.25;
my $FORWARD   =    0;
my $BACKWARD  =    1;
my $REVERSE   =    2;


my $table_file_name = "$user_result_dir/".$user_id."tmp.txt";
my $table_file_url = "$relative_user_result_dir/".$user_id."tmp.txt";
if (!open (TABLE_FILE, ">$table_file_name")) {
    print "<p>can\'t create the file $table_file_name. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777</p>";
}

#html file for table format
my $html_file_name = ">$user_result_dir/".$user_id."html_table.html";
my $html_table_file_url = "$relative_user_result_dir/".$user_id."html_table.html";
if(!open (HTML_TABLE_FILE, $html_file_name)) {
    print "<p>can\'t create the file $html_file_name. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777</p>";
}
print HTML_TABLE_FILE $query->start_html(
    -title => "Primer List",
    -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
);

# html file for primer design report. this file will be sent to user if an e-mail address is provided
my $html_report_file_name = ">$user_result_dir/".$user_id."_primer_report.html";
my $html_report_file_url = "$relative_user_result_dir/".$user_id."_primer_report.html";
if(!open (HTML_REPORT_FILE, $html_report_file_name)) {
    print "<p>can\'t create the file $html_report_file_name. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777</p>";
}
print HTML_REPORT_FILE $query->start_html(
    -title => "Primer design report",
    -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
);

# SSR screening results
my $ssr_screening_text_file_name;
my $ssr_screening_text_file_url;

my $ssr_screening_html_file_name;
my $ssr_screening_html_file_url;

if ($primer_type eq '2') {  # SSR screening and primer design
    $ssr_screening_text_file_name = "$user_result_dir/".$user_id."ssr.txt";
    $ssr_screening_text_file_url = "$relative_user_result_dir/".$user_id."ssr.txt";
    if (!open (SSR_TXT, ">$ssr_screening_text_file_name")) {
    	print "<p>can\'t create the file $ssr_screening_text_file_name. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777</p>";
    }

    #html file for table format
    $ssr_screening_html_file_name = ">$user_result_dir/$user_id" . "ssr.html";
    $ssr_screening_html_file_url = "$relative_user_result_dir/$user_id" . "ssr.html";
    if(!open (SSR_HTML, $ssr_screening_html_file_name)) {
        print "<p>can\'t create the file $ssr_screening_html_file_name. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777</p>";
    }

    print SSR_HTML $query->start_html(
        -title => "SSR Screening Results",
        -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
    );

}

#SSR screening
#motif-repeat parameters:
#specify motif length, minimum number of repeats.
#modify according to researcher's preferences
my @repeats = (
   [2,10],  #dinucl. with >= 10 repeats
   [3,7],   #trinucl. with >= 7 repeats
   [4,5],   #tetranucl. with >= 5 repeats
   [5,4],   #pentanucl. with >= 4 repeats
   [6,4]);  #hexanucl. with >= 4 repeats
my @selected = (0, 0, 0, 0, 0);
my @detected_ssrs = (0, 0, 0, 0, 0);
my $total_ssrs = 0;

my $total_seqs = 0;
my $total_primer_pairs = 0;
my $sucess_seqs = 0;
my $failed_seqs = 0;

my $snp_sucess_seqs = 0;
my $snp_failed_seqs = 0;
my $snp_total_primers = 0;

my $unique_id = 0;
my $table_row_count = 0;
my $ssr_table_row_count = 0;

my $first_base_index = 1;


# process input sequences and pick primers
if ($query->param('Pick Primers')) {
    print qq {
        <table class="standard">
             <tbody><tr><td>
             
             <h2>BatchPrimer3 is working on primer picking.....</h2>\n
        	<p>This will take several seconds to minutes depending on numbers of sequences a user submited. Please click the following
            link to check the primer design report. If you provide an email address, the results will be sent to your email account.
            You can just leave it and then check your email. </p>
            <p><a target="_blank" href = "$html_report_file_url">$html_report_file_url</a></p>
       

    };
    &process_input_and_pick_primers($query);

    print "<h2><font color=green>The job is finished. Please click the above link to check the results.</font></h2>";

    # redirect to the primer design result page using javascript. Perl's solution would not work
    print "<script type=\"text/javascript\">window.location.href='";
    print $html_report_file_url;
    print "';</script>";
 
    print qq {
        </td></tr></tbody></table>\n
        $wrapup\n
    };

}     


################# Main program ends ##################################


sub check_server_side_configuration {
    my ($query) = @_;

    unless (-e $PRIMER_BIN) {
	print qq{Please contact webmaster: cannot find $PRIMER_BIN executable
		     $wrapup};
	exit;
    }
    unless (-x $PRIMER_BIN) {
	print qq{Please contact webmaster: wrong permissions for $PRIMER_BIN
		     $wrapup};
	exit;
    }

    # Check mispriming / mishyb library setup.
    my @names = $query->param;
    for (@names) {
	if (/^PRIMER_(MISPRIMING|INTERNAL_OLIGO_MISHYB)_LIBRARY$/) {
	    my $v = $query->param($_);
	    my $v1 = $SEQ_LIBRARY{'$v'};
	    if (!defined($v)) {
		print qq{
		    <h3>There is a configuration error at $tmpurl;
		    cannot find a library file name for  "$v1".  Please clip and
		    mail this page to $MAINTAINER.$wrapup</h3>
			};
		exit;
	    }
	    if ($v1 && ! -r $v1) {
		print qq{
		    <h3>There is a configuration error at $tmpurl;
		    library file  $v1 cannot be read.  Please clip and
		    mail this page to $MAINTAINER.$wrapup</h3>
		    };
		exit;

	    }
	}
    }
}

sub print_output_header() {
    # $primer_type is a global variable
    print HTML_REPORT_FILE qq{
         <table class="standard">
         <tbody><tr><td>
         
         <h2>Report of BatchPrimer3 Primer Design</h2>\n
	     <p><p>
             <h4><a target=_blank href="$HTDOC/batchprimer_results_help.html">Help</a> | 
             <a target=_blank href=$html_table_file_url>View primers in HTML table format</a> |
             <a target=_blank href=$table_file_url>View primers in tab-delimited table format</a> <br>
             <a href="$CGI_URL/save_results.cgi?TYPE=txt&FILE=$table_file_name">Save as a tab-delimited text file</a> |
             <a href="$CGI_URL/save_results.cgi?TYPE=excel&FILE=$table_file_name">Save as an Excel file</a> |
             <a name=STATISTICS href="#STATISTICS">Primer Report Statistics</a> |  
             <a href="$relative_user_result_dir.zip">Download entire results (a zip file)</a> </h4>
             <p>
    };
    
    if ($primer_type eq '2') {
        print HTML_REPORT_FILE qq{
             <h4>
             <a target=_blank href=$ssr_screening_html_file_url>View SSR screening results in HTML table format</a> |
             <a target=_blank href=$ssr_screening_text_file_url>View SSR screening results in tab-delimited table format</a> <br>
             <a href="$CGI_URL/save_results.cgi?TYPE=txt&FILE=$ssr_screening_text_file_name">Save SSR screening results as a tab-delimited text file</a> |
             <a href="$CGI_URL/save_results.cgi?TYPE=excel&FILE=$ssr_screening_text_file_name"save_results.cgi">Save SSR screening resultsas an Excel file</a>
             </h4>
             <p>
        };
    
    }
    
    print HTML_REPORT_FILE "<hr>\n";
    
    # get the primer ype
    # 1. General primers
    # 2. SSR screening and primers
    # 3. Hybridization oligo primers
    print HTML_REPORT_FILE "Primer type: $primer_type_descriptions[$primer_type-1]<p>\n";

    # print the header of result table
    # print header of html table file
    my $title = $primer_type_descriptions[$primer_type-1];
    
    print HTML_TABLE_FILE "<h1>$title</h1>\n";
    print HTML_TABLE_FILE "<h3><a href=\"$CGI_URL/save_results.cgi?TYPE=txt&FILE=$table_file_name\">Save as a tab-delimited text file</a>&nbsp;&nbsp;|";
    print HTML_TABLE_FILE "<a href=\"$CGI_URL/save_results.cgi?TYPE=excel&FILE=$table_file_name\">&nbsp;&nbsp;Save as an Excel file</a></h3>\n";
    print HTML_TABLE_FILE "<table class='standard'><thead><tr>\n";
    
    
    if ($primer_type eq '2') {  #SSR  
        print SSR_HTML "<h1>SSR Screening Results</h1>\n";
        print SSR_HTML "<h3><a href=\"$CGI_URL/save_results.cgi?TYPE=txt&FILE=$ssr_screening_text_file_name\">Save as a tab-delimited text file</a>&nbsp;&nbsp;|";
        print SSR_HTML "<a href=\"$CGI_URL/save_results.cgi?TYPE=excel&FILE=$ssr_screening_text_file_name\">&nbsp;&nbsp;Save as an Excel file</a></h3>\n";
        print SSR_HTML "<table class='standard'><thead><tr>\n";

   	print TABLE_FILE "Index\tSeq ID\tCount\tPrimer type\tOrientation\tStart\tLen\ttm\tGC%\tAny compl\t3' compl\tSeq\tProd size\tSeq len\tIncluded len\tPair any compl\tPair 3' compl\tMotif\tMotif Len\tSSR\tSSRLen\n";
        print HTML_TABLE_FILE "<th>Index</th><th>Seq ID</th><th>Count</th><th>Primer type</th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Seq</th><th>Prod Size</th><th>Seq Length</th><th>Included Length</th><th>Pair any compl</th><th>Pair 3' compl</th><th>Motif</th><th>Motif len</th><th>SSR</th><th>SSR len</th>\n";
        
   	print SSR_TXT "Index\tSeq ID\tCount\tSSR start\tSSR end\tMortif\tMortif len\tSSR\tSSR len\n";
        print SSR_HTML "<th>Index</th><th>Seq ID</th><th>Count</th><th>SSR start</th><th>SSR end</th><th>Motif</th><th>Mortif len</th><th>SSR</th><th>SSR len</th>\n";
    }
    elsif ($primer_type eq '3') {   #oligo design
    	print TABLE_FILE "Index\tSeq ID\tCount\tPrimer_type\tstart\tlen\tTm\tGC\tany\t3'\tSeq\tSeq Length\tIncluded Length\n";
        print HTML_TABLE_FILE "<th>Index</th><th>Seq ID</th><th>Count</th><th>Primer type</th><th>start</th><th>len</th><th>Tm</th><th>GC</th><th>any</th><th>3'</th><th>Seq</th><th>Seq Length</th><th>Included Length</th>\n";
    }
    elsif ($primer_type >= 5 && $primer_type <= 9) {   #SBE primer only or allele-specifc primer only, #flanking primer and SBE/allele-specific primer/tetra-primer ARMS PCR
    	print TABLE_FILE "Index\tSeq ID\tCount\tPrimer_type\tOrientation\tStart\tLength\tTm\tGC%\tAny compl\t3' compl\tQ Score\tSNP\tPos\tPrimer Seq\tProd Size\tSeq Length\tIncluded Length\tPair any compl\tPair 3' compl\n";
        print HTML_TABLE_FILE "<th>Index</th><th>Seq ID</th><th>Count</th><th>Primer type</th><th>Orientation</th><th>Start</th><th>Length</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Q Score</th><th>SNP</th><th>Pos</th><th>Primer Seq</th><th>Prod Size</th><th>Seq Length</th><th>Included Length</th><th>Pair Any compl</th><th>Pair 3' compl</th>\n";
    }
    elsif ($primer_type == 10) {   #Sequencing primers
    	print TABLE_FILE "Index\tSeq ID\tCount\tPrimer_type\tOrientation\tStart\tLength\tTm\tGC%\tAny compl\t3' compl\tQ Score\tPrimer Seq\n";
        print HTML_TABLE_FILE "<th>Index</th><th>Seq ID</th><th>Count</th><th>Primer type</th><th>Orientation</th><th>Start</th><th>Length</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Q Score</th><th>Primer Seq</th>\n";
    }
    
    else {  # 1 or 4
   	print TABLE_FILE "Index\tSeq ID\tCount\tPrimer_type\tOrientation\tStart\tLen\ttm\tGC%\tAny compl\t3' compl\tSeq\tProd Size\tSeq Length\tIncluded Length\tPair any compl\tPair 3' compl\n";
        print HTML_TABLE_FILE "<th>Index</th><th>Seq ID</th><th>Count</th><th>Primer_type</th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Seq</th><th>Prod Size</th><th>Seq Length</th><th>Included Length</th><th>Pair any compl</th><th>Pair 3' compl</th>\n";
    }

    print HTML_TABLE_FILE "</tr></thead>\n";
    if ($primer_type eq '2') {    
        print SSR_HTML "</tr></thead>\n";
    }
}
    
# This function is to run primer3 core program and return the primer design results.
sub run_primer3{
    my ($cmd, $input) = @_;

    my $primer3_pid;
    my %results = ();
    
    my ($childin, $childout, $childerr) = (FileHandle->new, FileHandle->new, FileHandle->new);
    local $^W = 0;
    $primer3_pid = open3($childin, $childout, $childerr, $cmd);
    if (!$primer3_pid) {
        print "Cannot excecure $cmd:<br>$!";
        print "$wrapup\n";
        exit;
    }

    print $childin @$input;
    $childin->close;
  
    my @results = $childout->getlines;
  
    waitpid $primer3_pid, 0;
    return @results;
}

# process user input and pick primers for all types of primers
sub process_input_and_pick_primers {
    my $start_time = time();
    my ($query) = @_;
    my ($v, $v1);

    &print_output_header();
    

    &check_server_side_configuration($query);

    my @names = $query->param;
    my $cmd = "$PRIMER_BIN -format_output -strict_tags";
    my $line;
    my $fasta_id;


    #save the all parameters in a file
    &save_primer_parameters($query, \@names);

        
    # SSR parameters
    if ($primer_type eq '2') {    
        if ($query->param('DINUCLEOTIDE_SSR')) {
        	$selected[0] = 1;
        }
        if ($query->param('TRINUCLEOTIDE_SSR')) {
        	$selected[1] = 1;
    	}
        if ($query->param('TETRANUCLEOTIDE_SSR')) {
        	$selected[2] = 1;
        }
        if ($query->param('PENTANUCLEOTIDE_SSR')) {
        	$selected[3] = 1;
        }
        if ($query->param('HEXANUCLEOTIDE_SSR')) {
        	$selected[4] = 1;
        }
    
        if ($query->param('DINUCLEOTIDE_SSR') && $query->param('DINUCLEOTIDE_SSR_REPEATS')) {
        	$repeats[0]->[1] = $query->param('DINUCLEOTIDE_SSR_REPEATS');
        }
        if ($query->param('TRINUCLEOTIDE_SSR') && $query->param('TRINUCLEOTIDE_SSR_REPEATS')) {
        	$repeats[1]->[1] = $query->param('TRINUCLEOTIDE_SSR_REPEATS');
    	}
        if ($query->param('TETRANUCLEOTIDE_SSR') && $query->param('TETRANUCLEOTIDE_SSR_REPEATS')) {
        	$repeats[2]->[1] = $query->param('TETRANUCLEOTIDE_SSR_REPEATS');
        }
        if ($query->param('PENTANUCLEOTIDE_SSR') && $query->param('PENTANUCLEOTIDE_SSR_REPEATS')) {
        	$repeats[3]->[1] = $query->param('PENTANUCLEOTIDE_SSR_REPEATS');
        }
        if ($query->param('HEXANUCLEOTIDE_SSR') && $query->param('HEXANUCLEOTIDE_SSR_REPEATS')) {
        	$repeats[4]->[1] = $query->param('HEXANUCLEOTIDE_SSR_REPEATS');
        }
    }
    
    $SNP_PRIMER_MIN_GC  = $query->param('SNP_PRIMER_MIN_GC');
    $SNP_PRIMER_MAX_GC  = $query->param('SNP_PRIMER_MAX_GC');
    $SNP_PRIMER_MIN_TM  = $query->param('SNP_PRIMER_MIN_TM');
    $SNP_PRIMER_OPT_TM  = $query->param('SNP_PRIMER_OPT_TM');
    $SNP_PRIMER_MAX_TM  = $query->param('SNP_PRIMER_MAX_TM');
    $SNP_PRIMER_MIN_SIZE  = $query->param('SNP_PRIMER_MIN_SIZE');
    $SNP_PRIMER_OPT_SIZE  = $query->param('SNP_PRIMER_OPT_SIZE');
    $SNP_PRIMER_MAX_SIZE  = $query->param('SNP_PRIMER_MAX_SIZE');
    $SNP_PRIMER_MAX_N  = $query->param('SNP_PRIMER_MAX_N');
    $SNP_PRIMER_SALT_CONC = $query->param('SNP_PRIMER_SALT_CONC');
    $SNP_PRIMER_DNA_CONC = $query->param('SNP_PRIMER_DNA_CONC');
    $SNP_MAX_SELF_COMPLEMENTARITY  = $query->param('SNP_MAX_SELF_COMPLEMENTARITY');  
    $SNP_MAX_3_SELF_COMPLEMENTARITY= $query->param('SNP_MAX_3_SELF_COMPLEMENTARITY');

    $PRIMER_MAX_DIFF_TM = $query->param('PRIMER_MAX_DIFF_TM');
    
    $SNP_INNER_PRODUCT_MIN_SIZE = $query->param('SNP_INNER_PRODUCT_MIN_SIZE');
    $SNP_INNER_PRODUCT_OPT_SIZE = $query->param('SNP_INNER_PRODUCT_OPT_SIZE');
    $SNP_INNER_PRODUCT_MAX_SIZE = $query->param('SNP_INNER_PRODUCT_MAX_SIZE');
    $SNP_INNER_PRODUCT_MIN_DIFF = $query->param('SNP_INNER_PRODUCT_MIN_DIFF');
    $SNP_INNER_PRODUCT_MAX_DIFF = $query->param('SNP_INNER_PRODUCT_MAX_DIFF');

    if ($primer_type eq '7' || $primer_type eq '8') {
        if ($query->param('SECOND_MISMATCH')) {
            $SECOND_MISMATCH = 1;
            $SECOND_MISMATCH_POS = $query->param('SECOND_MISMATCH_POS');
        }
    }


    if ($primer_type eq '10') {
        $SEQUENCING_DIRECTION = $query->param('SEQUENCING_DIRECTION');
        $PRIMER_NUM_RETURN = $query->param('PRIMER_NUM_RETURN');
    }
    
    my $size_range='';
    # my $sequence_id = $query->param('PRIMER_SEQUENCE_ID');
    my $sequence_id = '';   # sequence_id is obtained from the input sequence
    my $min_prod_size = $PR_DEFAULT_PRODUCT_MIN_SIZE;
    my $max_prod_size = $PR_DEFAULT_PRODUCT_MAX_SIZE;

    if ($primer_type ne '3' && $primer_type ne '6' && $primer_type ne '8' &&  $primer_type ne '9' &&  $primer_type ne '10') {
        $min_prod_size = $query->param('MUST_XLATE_PRODUCT_MIN_SIZE');
        $max_prod_size = $query->param('MUST_XLATE_PRODUCT_MAX_SIZE');
        $min_prod_size = $PR_DEFAULT_PRODUCT_MIN_SIZE unless $min_prod_size =~ /\S/;
        $max_prod_size = $PR_DEFAULT_PRODUCT_MAX_SIZE unless $max_prod_size =~ /\S/;
    	$size_range = "$min_prod_size-$max_prod_size";
	}

    $first_base_index = $query->param('PRIMER_FIRST_BASE_INDEX');
    if (!$first_base_index || $first_base_index !~ /\d/) {
		$first_base_index = 1;
    }



    my $pick_left  = $query->param('MUST_XLATE_PICK_LEFT');
    my $pick_hyb   = $query->param('MUST_XLATE_PICK_HYB_PROBE');
    my $pick_right = $query->param('MUST_XLATE_PICK_RIGHT');

    $pick_left  = 1 if $query->param('PRIMER_LEFT_INPUT');
    $pick_right = 1 if $query->param('PRIMER_RIGHT_INPUT');
    $pick_hyb   = 1 if $query->param('PRIMER_INTERNAL_OLIGO_INPUT');

    
    my $task;
    if ($primer_type eq '1' || $primer_type eq '2' || $primer_type eq '3' ||$primer_type eq '4' ||
        $primer_type eq '5' || $primer_type eq '7') {
        if ($pick_hyb) {
            if ($pick_right || $pick_left) {
                $task='pick_pcr_primers_and_hyb_probe';
                print "WARNING: Assuming you want to pick a right primer because\n",
    		"         you are picking a left primer and internal oligo\n"
                if !$pick_right;
                    print "WARNING: Assuming you want to pick a left primer because\n",
                	"         you are picking a righ primer and internal oligo\n"
                if !$pick_left;
            } else {
                $task='pick_hyb_probe_only';
            }
        } else {
            if ($pick_right && $pick_left) {
                $task='pick_pcr_primers';
            } elsif ($pick_right) {
                $task='pick_right_only';
            } elsif ($pick_left) {
                $task='pick_left_only';
            } else {
                print "WARNING: assuming you want to pick PCR primers\n";
                $task='pick_pcr_primers';
            }
        }
    }

    my $print_input = $query->param('MUST_XLATE_PRINT_INPUT');
    
    my $inferred_sequence = '';
    if (!$query->param('SEQUENCE')) {
	if ($query->param('PRIMER_LEFT_INPUT')) {
            $inferred_sequence = $query->param('PRIMER_LEFT_INPUT');
	}
	if ($query->param('PRIMER_INTERNAL_OLIGO_INPUT')) {
            $inferred_sequence .= $query->param('PRIMER_INTERNAL_OLIGO_INPUT');
	}
	if ($query->param('PRIMER_RIGHT_INPUT')) {
            my $tmp_seq = $query->param('PRIMER_RIGHT_INPUT');
            $tmp_seq = scalar(reverse($tmp_seq));
            $tmp_seq =~ tr/acgtACGT/tgcaTGCA/;
            $inferred_sequence .= $tmp_seq;
	}
    }


    my @input;
    push @input, "PRIMER_EXPLAIN_FLAG=1\n";
    $DO_NOT_PICK = 0;

    for (@names) {
	next if /^EMAIL$/;
	next if /^Pick Primers$/;
	next if /^old_obj_fn$/;
	next if /^PRIMER_SEQUENCE_ID$/;
	next if /^PRIMER_FIRST_BASE_INDEX$/;
	next if /^TARGET$/;
	next if /^INCLUDED_REGION$/;
	next if /^EXCLUDED_REGION$/;
	next if /^SEQUENCEFILE$/;
	next if /^SEQUENCE$/;
	next if /^PRIMER_TYPE$/;
	next if /^DINUCLEOTIDE_SSR$/;
	next if /^TRINUCLEOTIDE_SSR$/;
	next if /^TETRANUCLEOTIDE_SSR$/;
	next if /^PENTANUCLEOTIDE_SSR$/;
	next if /^HEXANUCLEOTIDE_SSR$/;
	next if /^DINUCLEOTIDE_SSR_REPEATS$/;
	next if /^TRINUCLEOTIDE_SSR_REPEATS$/;
	next if /^TETRANUCLEOTIDE_SSR_REPEATS$/;
	next if /^PENTANUCLEOTIDE_SSR_REPEATS$/;
	next if /^HEXANUCLEOTIDE_SSR_REPEATS$/;
	next if /^PICK_SSR$/;
                
        # SNP (allele-specific) primers
        next if /^SNP_PRIMER_MAX_GC$/;
        next if /^SNP_PRIMER_MIN_GC$/;
        next if /^SNP_PRIMER_MIN_TM$/;
        next if /^SNP_PRIMER_OPT_TM$/;
        next if /^SNP_PRIMER_MAX_TM$/;
        next if /^SNP_PRIMER_MIN_SIZE$/;
        next if /^SNP_PRIMER_OPT_SIZE$/;
        next if /^SNP_PRIMER_MAX_SIZE$/;
        next if /^SNP_INNER_PRODUCT_MIN_SIZE$/;
        next if /^SNP_INNER_PRODUCT_OPT_SIZE$/;
        next if /^SNP_INNER_PRODUCT_MAX_SIZE$/;
        next if /^SNP_INNER_PRODUCT_MIN_DIFF$/;
        next if /^SNP_INNER_PRODUCT_MAX_DIFF$/;
        next if /^SNP_PRIMER_MAX_N$/;
        next if /^SNP_PRIMER_SALT_CONC$/;
        next if /^SNP_PRIMER_DNA_CONC$/;
        next if /^SNP_MAX_SELF_COMPLEMENTARITY$/;  
        next if /^SNP_MAX_3_SELF_COMPLEMENTARITY$/;
                
        next if /^SEQUENCING_DIRECTION$/;

        next if /^SECOND_MISMATCH$/;
        next if /^SECOND_MISMATCH_POS$/;
                
  
        $v = $query->param($_);
	next if $v =~ /^\s*$/;   # Is this still the right behavior?


	if (/^MUST_XLATE/) {
	    next if /^MUST_XLATE_PRODUCT_M(IN|AX)_SIZE$/;
	    next if /^MUST_XLATE_PICK_(LEFT|RIGHT|HYB_PROBE)$/;
	    next if /^MUST_XLATE_PRINT_INPUT$/;
	} elsif (/^PRIMER_(MISPRIMING|INTERNAL_OLIGO_MISHYB)_LIBRARY$/) {
            $v = $SEQ_LIBRARY{$v};
	} elsif (/^PRIMER_MIN_SIZE$/ && $v < 1) {
	    print "WARNING: Changed illegal Primer Size Min of $v to 1\n";
	    $v = 1;
	}
	$line = "$_=$v\n";
	push @input, $line;
    }
    
    if ($DO_NOT_PICK) {
	print "$wrapup\n";
        print HTML_TABLE_FILE "$wrapup\n";
        close(TABLE_FILE);
        close(HTML_TABLE_FILE);
        return;
    }
    if ($primer_type eq '1' || $primer_type eq '2' || $primer_type eq '3' ||$primer_type eq '4' ||
        $primer_type eq '5' || $primer_type eq '7') {
        push @input, "PRIMER_TASK=$task\n";
    }
    
    if ($primer_type ne '3') {
    	push @input, "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
    }
    
    push @input, "PRIMER_FIRST_BASE_INDEX=$first_base_index\n";
    push @input, "SEQUENCE=$inferred_sequence\n" if $inferred_sequence;

    # SNP codes are allowed
    if ($primer_type >= 4 && $primer_type <= 9){
        push @input, "PRIMER_LIBERAL_BASE=1\n";
    }
    
    if ($primer_type ne '6' && $primer_type ne '8' && $primer_type ne '10') {
        push @input, "PRIMER_PICK_ANYWAY=1\n";
    }
    
       
    # Processing multiple sequences. For a sequence file, sequence format must be a FASTA format.
    # For sequences in text area, if there is only one sequence, sequence can be or not
    # FASTA format. If there is sequence file, the sequences in text area will be ignored.

    my $s = $query->param('SEQUENCE');
    my $s0 = $query->param('SEQUENCEFILE');
    if ($s && $s !~ /^\s*$/ || $s0 && $s0 !~ /^\s*$/) {
	my $sequences = "";
	if ($s0 || $s0 =~ /^\s*$/) {  # if no sequence file
            $sequences = $s;
	}

        my $first_sequence = 1;
	my $next_header = '';
	my $last_seq = 0;
 	while(!$last_seq) {
            $v = '';
            if ($s0 && $s0 !~ /^\s*$/) {  # for sequence file
                my $count = 0;
                $v = $next_header;
                if ($v !~ /^\s*$/) {  #there is header already
                    $count = 1;
                }

		while(<$s0>) {
                    chomp($_);
                    my $index1 = index ($_, '>');  # search the first '>'
                    if ($index1 >= 0) {   # found the next record
			$count++;         # $count should be 2
                        my $index2 = &find_excluded_index ($v . $_);
                        my $index3 = rindex($v . $_, '>');
                        if ($index2 == $index3) {
                            $v .= "$_";
                            $count--;
                            next;
                        }
                        
			if ($count == 2) {
                            $next_header = substr($_, $index1, length($_)-1)."\n";
                            if ($index1 > 0) {
				$v .= substr($_, 0, $index1 -1);
				$v .= "\n";
                            }
                            last;
                        }
		        elsif ($count == 1) {
                            $v .= "$_\n";
                        }
		    }else {
                        $v .= "$_";
                    }
		}
      
                if ($count == 1) {  # last sequence)
                    $last_seq = 1;
                }

	    } elsif ($s && $s !~ /^\s*$/) {   # no sequence file and there is sequencess in the text area)
        	my $index1 = index ($sequences, '>', 0);  # search the first '>'
		my $index2 = &find_index ($sequences, '>', $index1 + 1);  # search the second '>'
 		if ($index2 >= 0) {
                    $v = substr($sequences, $index1, $index2 - $index1);
                    $sequences = substr($sequences, $index2, length($sequences)-1);
		} else {
                    $v = $sequences;
                    $last_seq = 1;
		}
                
                
	    } else {
		print "No sequence input!!";
    	  	print  "$wrapup\n";
                print HTML_TABLE_FILE "$wrapup\n";
                close(TABLE_FILE);
                close(HTML_TABLE_FILE);
                last;
            }

	    next if $v =~ /^\s*$/;
            $total_seqs++;
            if ($v =~ /^\s*>([^\n]*)/) {
		# Sequence is in Fasta format.
		$fasta_id = $1;
		$fasta_id =~ s/^\s*//;
		$fasta_id =~ s/\s*$//;
		if (!$sequence_id) {
                    $sequence_id = $fasta_id;
		} else {
                    $sequence_id = $fasta_id;
		}
		$v =~ s/^\s*>([^\n]*)//;
	    }

            ($sequence_id) = $sequence_id =~ /^(\S+)/;

	    if ($v =~ /\d/) {
		print "WARNING: Numbers in input sequence were deleted.\n";
		$v =~ s/\d//g;
	    }
            $v =~ s/\s//g;
            my ($m_target, $m_excluded_region, $m_included_region)
		= read_sequence_markup($v, (['[', ']'], ['<','>'], ['{','}']));
             $v =~ s/[\[\]\<\>\{\}]//g;
  
  
            my $target = undef;
            my $excluded_region = undef;
            my $included_region = undef;
              
            if (@$m_target) {
		if ($target) {
                    print "WARNING Targets specified both as sequence ",
                          "markups and in Other Per-Sequence Inputs\n";
		}
		$target = add_start_len_list($target, $m_target, $first_base_index);
            }
            
            if (@$m_excluded_region) {
                if ($excluded_region) {
                    print "WARNING Excluded Regions specified both as sequence ",
                        "markups and in Other Per-Sequence Inputs\n";
                }
                $excluded_region = add_start_len_list($excluded_region,
		    $m_excluded_region,
		    $first_base_index);
            }
	
            if (@$m_included_region) {
                if (scalar @$m_included_region > 1) {
                    print "ERROR: Too many included regions\n";
                    $DO_NOT_PICK = 1;
                } elsif ($included_region) {
                    print "ERROR: Included region specified both as sequence\n",
                        "       markup and in Other Per-Sequence Inputs\n";
                    $DO_NOT_PICK = 1;
                }
            
                $included_region = add_start_len_list($included_region,
                    $m_included_region,
                    $first_base_index);
            }

            for (my $i = @input-1; $i >= 0; $i--) {
                if ($input[$i] =~ m/EXCLUDED_REGION=/) {
                    splice(@input, $i, 1);
                    $i--;
                }
                if ($input[$i] =~ m/INCLUDED_REGION=/) {
                    splice(@input, $i, 1);
                    $i--;
                }
            }
           
            splice(@input, 0, 0, "EXCLUDED_REGION=$excluded_region\n") if $excluded_region;
            splice(@input, 0, 0, "INCLUDED_REGION=$included_region\n") if $included_region;
           

            if ($primer_type eq '1' || $primer_type eq '3') {
                if ($first_sequence == 1) {
                    push(@input, "PRIMER_SEQUENCE_ID=$sequence_id\n");
                    push(@input, "TARGET=$target\n") if $target;
                    push(@input, "SEQUENCE=$v\n");
                    push(@input, "=\n");
                } elsif ($first_sequence == -1) {
                    for (my $i = @input-1; $i >= 0; $i--) {
                        if ($input[$i] =~ m/PRIMER_SEQUENCE_ID=/) {
                            $input[$i] = "PRIMER_SEQUENCE_ID=$sequence_id\n";
                        }
                        if ($input[$i] =~ m/TARGET=/) {
                            if ($target) {
                                $input[$i] = "TARGET=$target\n"
                            } else {
                                splice(@input, $i, 1);
                                $i--;
                            }
                        }
                        if ($input[$i] =~ m/SEQUENCE=/) {
                            $input[$i] = "SEQUENCE=$v\n";
                        }
                    }
                }
                if ($first_sequence == 1) {
                    $first_sequence = -1;
                }

		$unique_id++;
		
                my @results = &run_primer3($cmd, \@input);
                use Primer3Output;
                use PrimerPair;
                use Primer;
                my $primer3output = new Primer3Output();
                $primer3output->set_results(\@results);
                        
                &print_primer_results($included_region, $excluded_region, $target, $v,
                    $primer_type, $primer3output, $print_input, $sequence_id, \@names, \@input);

            } elsif ($primer_type eq '2') {  # SSR primers
            	# return references of arrays
            	my ($targets, $prod_size_ranges, $motifs, $ssrs) = find_ssrs($v, $included_region, $excluded_region);
                
                my $sucess = 0;
                if (@$targets) {
                    # print SSR screening results
                    &print_SSR_screening_results($sequence_id, $motifs, $ssrs, $targets);
                    
		    for (my $i = 0; $i < scalar(@$targets); $i++) {
                       if ($first_sequence == 1) {
                            push(@input, "PRIMER_SEQUENCE_ID=$sequence_id\n");
                            push(@input, "TARGET=@$targets[$i]\n") if @$targets[$i];
                            push(@input, "SEQUENCE=$v\n");
                            push(@input, "=\n");
                        } elsif ($first_sequence == -1) {
                            for (my $j = @input-1; $j > @input -5; $j--) {
                                if ($input[$j] =~ m/PRIMER_SEQUENCE_ID=/) {
                                    $input[$j] = "PRIMER_SEQUENCE_ID=$sequence_id\n";
                                }
                            if ($input[$j] =~ m/TARGET=/) {
                                if (@$targets[$i]) {
                                    $input[$j] = "TARGET=@$targets[$i]\n"
                                } else {
                                    splice(@input, $j, 1);
                                }
                            }
                            if ($input[$j] =~ m/SEQUENCE=/) {
                                    $input[$j] = "SEQUENCE=$v\n";
                                }
                            }
                        }

                        # count the total numbers of detected SSRs and number of different patterns
                        $total_ssrs++;
                        $detected_ssrs[length(@$motifs[$i])-2]++;
            
                        $target = @$targets[$i];
                
                        # create product size range for SSR primer
                        my $len = length(@$ssrs[$i]);
                        my $size_range;
                        my $min_size = $min_prod_size;
                        my $max_size = $max_prod_size;
                        if ($min_prod_size > $len) {
                           $min_size = $len - 10;
                        }

                        if ($max_prod_size < $len) {
                        	$max_size = $len + 10;
                        }
                        if ($min_prod_size > $len || $max_prod_size < $len) {
                            $size_range = "$min_size-$max_size";
                            for (my $j = @input-1; $j > @input -5; $j--) {
                                if ($input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                                    $input[$j] = "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
                                    last;
                                }
                            }
                        }
                        if ($first_sequence == 1) {
                            $first_sequence = -1;
                        }
			
			$unique_id++;
		
			my @results = &run_primer3($cmd, \@input);
			use Primer3Output;
			use PrimerPair;
			use Primer;
			my $primer3output = new Primer3Output();
			$primer3output->set_results(\@results);
			&print_primer_results($included_region, $excluded_region, $target, $v,
                            $primer_type, $primer3output, $print_input, $sequence_id, \@names, \@input, @$motifs[$i], @$ssrs[$i]);
                        
                        $sucess = 1 if ($primer3output->get_primer_list_size() > 0);
                    
                    }
                        
                    if ($sucess == 1) {    
                        $sucess_seqs++;
                    } else {
                        $failed_seqs++;
                    }

                }
                else {
                    $failed_seqs++;
                }
            } elsif ($primer_type >= 4 && $primer_type <= 9) {  # SNP flanking, single base extension and alle-specific primers
                my $direction = " ";
                my @primer_results = ();
                my $snp = "";
                my $pos = 0;


                my $snp_header = " <table><thead><tr><th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Score</th><th> SNP</th><th>Pos</th><th>Primer Seq</th></tr></thead>\n";

            	# return references of arrays
            	my $targets = find_snp_targets($v, $first_base_index);
                if (@$targets) {
                    for (my $i = 0; $i < scalar(@$targets); $i++) {
                        $unique_id++;
                        if ($first_sequence == 1) {
                            push(@input, "PRIMER_SEQUENCE_ID=$sequence_id\n");
                            push(@input, "TARGET=@$targets[$i]\n") if @$targets[$i];
                            push(@input, "SEQUENCE=$v\n");
                            push(@input, "=\n");
                        } elsif ($first_sequence == -1) {
                            for (my $j = @input-1; $j > 0; $j--) {
                                if ($input[$j] =~ m/PRIMER_SEQUENCE_ID=/) {
                                    $input[$j] = "PRIMER_SEQUENCE_ID=$sequence_id\n";
                                }
                                if ($input[$j] =~ m/TARGET=/) {
                                    if (@$targets[$i]) {
                                        $input[$j] = "TARGET=@$targets[$i]\n"
                                    } else {
                                        splice(@input, $j, 1);
                                        $j--;
                                    }
                                }
                                if ($input[$j] =~ m/SEQUENCE=/) {
                                    $input[$j] = "SEQUENCE=$v\n";
                                }
                            }
                        }
		        if ($first_sequence == 1) {
                    	    $first_sequence = -1;
                	}
                        
                        $target = @$targets[$i];
                                                
                        my $std_primer_found = 0;
                        
                        if ($primer_type eq '4') {  # flanking primers
			    my @results = &run_primer3($cmd, \@input);
			    use Primer3Output;
			    use PrimerPair;
			    use Primer;
			    my $primer3output = new Primer3Output();
			    $primer3output->set_results(\@results);
        
     			    &print_primer_results($included_region, $excluded_region, $target, $v,
                                $primer_type, $primer3output, $print_input, $sequence_id, \@names, \@input);
                        }
                        
                        elsif ($primer_type eq '5') { # flanking primers and SBEs
 
                            # design SBE first
                            my $include_snp = 0;
                            my ($results_f, $results_r) = &find_SBE_primers ($included_region, $excluded_region, $v, @$targets[$i], $include_snp);
                            my @results_all_f = ($results_f);
                            my @results_all_r = ($results_r);
                            
                            # design flanking primers
                            my @results = &run_primer3($cmd, \@input);
			    use Primer3Output;
			    use PrimerPair;
			    use Primer;
			    my $primer3output = new Primer3Output();
			    $primer3output->set_results(\@results);
                            
                            # print all results together
                            &print_primer_results_for_flanking_and_SBE_alleles(
                                $included_region, $excluded_region, $target, $v,
                                $primer_type, $primer3output, $print_input, $sequence_id,
                                \@names, \@input, \@results_all_f, \@results_all_r);
 
                        }
                        
                        elsif ($primer_type eq '6') {  # SBE primers only                          
                            my $include_snp = 0;
                            my ($results_f, $results_r) = &find_SBE_primers ($included_region, $excluded_region, $v, @$targets[$i], $include_snp);
                            my @results_all_f = ($results_f);
                            my @results_all_r = ($results_r);
 
                            &print_primer_results_for_SBE_alleles(
                                $included_region, $excluded_region, $target, $v,
                                $primer_type, $sequence_id,
                                 \@results_all_f, \@results_all_r);
 
                        } elsif ($primer_type == '7') {  # allele-specifc primer and flanking primers
                            
                            # $results_f_all, $results_r_all are two-dimentional arrays
                            my $include_snp = 1;
                            my $second_mismatch = $SECOND_MISMATCH;
                            my $second_mismatch_pos = $SECOND_MISMATCH_POS;
                            my ($results_f_all, $results_r_all) = &find_allele_primers($included_region, $excluded_region, $v, @$targets[$i],
                                $include_snp, $second_mismatch, $second_mismatch_pos, $primer_type);
                            
                            # design flanking primers
            		    my @results = &run_primer3($cmd, \@input);
                            use Primer3Output;
                            use PrimerPair;
                            use Primer;
                            my $primer3output = new Primer3Output();
                            $primer3output->set_results(\@results);

                            &print_primer_results_for_flanking_and_SBE_alleles(
                                $included_region, $excluded_region, $target, $v,
                                $primer_type, $primer3output, $print_input, $sequence_id,
                                \@names, \@input, $results_f_all, $results_r_all);
                            
                        }elsif ($primer_type eq '8') { # allele-specific primers only 
                            ($snp, $pos) = get_snp($v, @$targets[$i]);
                            
                            # $results_f_all, $results_r_all are two-dimentional arrays
                            my $include_snp = 1;
                            my $second_mismatch = $SECOND_MISMATCH;
                            my $second_mismatch_pos = $SECOND_MISMATCH_POS;
                            my ($results_f_all, $results_r_all) = &find_allele_primers($included_region, $excluded_region, $v, @$targets[$i],
                                $include_snp, $second_mismatch, $second_mismatch_pos, $primer_type);
                            &print_primer_results_for_SBE_alleles(
                                $included_region, $excluded_region, $target, $v,
                                $primer_type, $sequence_id,
                                 $results_f_all, $results_r_all, $snp, $pos);
 
                        } elsif ($primer_type eq '9') {  # tetra-primer ARMS PCR primers
                            # Step 1: design four optimal allele-specific primers for forward and reverse orientations
                            my $include_snp = 1;
                            my $second_mismatch = 1;
                            my $second_mismatch_pos = -3;
                            # $results_f_all, $results_r_all are two-dimentional arrays
                            my ($results_f_all, $results_r_all) = &find_allele_primers($included_region, $excluded_region, $v, @$targets[$i],
                                $include_snp, $second_mismatch, $second_mismatch_pos, $primer_type);
                            
                            my ($primer3outputs, $results_f, $results_r) = 
                                &design_ARMS_outer_primers ($included_region, $excluded_region, $target, $v,
                                    $results_f_all, $results_r_all, \@input, $cmd, $sequence_id);

                            if ($primer3outputs) {
                                
                                my $found = 0;
                                for (my $m = 0; $m < @$primer3outputs; $m++) {
                                    next if (@$primer3outputs[$m]->get_primer_list_size() == 0);
                                    my @r_f = (@$results_f[$m]);
                                    my @r_r = (@$results_r[$m]);
                                    if ($found == 0) {
                                        $sucess_seqs++;
                                        $found = 1;
                                    }
                                    &print_primer_results_for_flanking_and_SBE_alleles(
                                        $included_region, $excluded_region, $target, $v,
                                        $primer_type, @$primer3outputs[$m], $print_input, $sequence_id,
                                        \@names, \@input, \@r_f, \@r_r);
                                    $unique_id++;
                                    $total_primer_pairs++;
                                }
                            } else {
                                $failed_seqs++;
                            }
                        }     
                       
                    }
                } else {
                    if ($primer_type >= 4) {
                        print "<p><b>Error:</b><font color=red> No SNP is masked in the input sequence. ".
                            "IUB/IUPAC nucleic acid codes are used to <br>represent SNPs or alleles. ".
                            "Please see the example sequences.</font><p/>";
                    }
                }
            }
            elsif ($primer_type == 10) {  # sequencing primer
                my $snp_header = " <table><thead><tr><th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Score</th><th>Primer Seq</th></tr></thead>\n";
            	$unique_id++;

                my @results = &find_sequencing_primers ($included_region, $excluded_region, $v);
                &print_primer_results_for_sequencing_primers($included_region, $excluded_region, $v, $primer_type, $sequence_id, \@results, $target);
            }
	} # end while
    } # end if
    
    
    &print_summary_and_footer($start_time);
    
    # write user access information to log file
    &log($ip,$total_seqs, $primer_type_descriptions[$primer_type-1]);
    
    # send email to the user if user's email is provided
    if (defined($email) && $email !~ /^\s*$/) {
        &send_email();
    }
    
    # create a zip file for the result directory
    &create_zip_file ($dir_name, $user_result_dir);
 
}

#send email using UNIX/Linux sendmail
sub send_email() {
    
    unless (Email::Valid->address($email)) {
        print "You supplied an invalid email address: $email.";
        return;
    }
    
    open(MAIL, "$SENDMAIL_CMD -t");

    my $to = $email;
    my $from = $MAIL_FROM;
    my $subject = "Results of batch primer design for $primer_type_descriptions[$primer_type-1]";
    my $body = "<p>Dear User:</p><p>Thank you for using BatchPrimer3 to design primers.\n\n"
               . "Please click the following link to access the results:<p>"
               . "<p><a href=$html_report_file_url>$html_report_file_url</a><p>"
               . "<p>From the result page, you can <br>"
               . " (1) save the primer list as tab-delimited text file or Excel file;<br>"
               . " (2) view the primer de detail for each sequence; <br>"
               . " (3) and download all results including primer design report, primer list, and primer design view files to local machine</p>"
               . "<p><p>"
               . "<p> The results will be retained in one week.</p>"
	       . "<p></p><p></p>"
               . "<p>BatchPrimer3 development team</p>";
 
    ## Mail Header
    print MAIL "To: $to\n";
    print MAIL "From: $from\n";
    print MAIL "Subject: $subject\n";
    print MAIL "Content-type: text/html\n\n";

    ## Mail Body
    print MAIL $body;

    close(MAIL);
    
 # altanative approach to send email    
 #   use Mail::Mailer;
 #   my $type = 'sendmail';
 #   my $mailprog = Mail::Mailer->new($type);

    # mail headers to use in the message
 #   my %headers = (
 #       'To' => $to,
 #       'From' => $from,
 #       'Subject' => $subject
 #   );

 #   $mailprog->open(\%headers);
 #   print $mailprog "$body\n";
 #   $mailprog->close;
    
}


# print primer design summary on the primer design report page
sub print_summary_and_footer() {
    my $start_time = shift;
    
    print HTML_REPORT_FILE "<hr>";
    print HTML_REPORT_FILE "<h3><a name=STATISTICS>Primer Report Statistics</a></h3>";
    print HTML_REPORT_FILE "Total sequences input:  <b>$total_seqs</b><p>";

    if ($primer_type eq '2') {
	print HTML_REPORT_FILE "Total SSRs detected: <b>$total_ssrs</b><p>";
        print HTML_REPORT_FILE "number of dinucleotide SSRs detected: <b>$detected_ssrs[0]</b><p>";
        print HTML_REPORT_FILE "number of trinucleotide SSRs detected: <b>$detected_ssrs[1]</b><p>";
        print HTML_REPORT_FILE "number of tetranucleotide SSRs detected: <b>$detected_ssrs[2]</b><p>";
        print HTML_REPORT_FILE "number of pentanucleotide SSRs detected: <b>$detected_ssrs[3]</b><p>";
        print HTML_REPORT_FILE "number of hexanucleotide SSRs detected: <b>$detected_ssrs[4]</b><p>";
    }

    if ($primer_type ne '6' && $primer_type ne '9' && $primer_type ne '10' ) {
	print HTML_REPORT_FILE "Number of sequences with sucessful primer pairs: <b>$sucess_seqs</b><p>";
        print HTML_REPORT_FILE "Number of sequences without primer pair picked: <b>$failed_seqs</b><p>";
        print HTML_REPORT_FILE "Total primer pairs picked: <b>$total_primer_pairs</b><p>";
    }
    
    if ($primer_type eq '9') {
	print HTML_REPORT_FILE "Number of sequences with sucessful primer sets: <b>$sucess_seqs</b><p>";
        print HTML_REPORT_FILE "Number of sequences without primer picked: <b>$failed_seqs</b><p>";
        print HTML_REPORT_FILE "Total primer sets picked: <b>$total_primer_pairs</b><p>";
    }

    if ($primer_type eq '10') {
	print HTML_REPORT_FILE "Number of sequences with sucessful primers: <b>$sucess_seqs</b><p>";
        print HTML_REPORT_FILE "Number of sequences without primer picked: <b>$failed_seqs</b><p>";
        print HTML_REPORT_FILE "Total primers picked: <b>$total_primer_pairs</b><p>";
    }


    if ($primer_type >= 5) {
        if ($primer_type eq '7') {
            print HTML_REPORT_FILE "Number of sequences with sucessful allele-specific primers: <b>$snp_sucess_seqs</b><p>";
            print HTML_REPORT_FILE "Number of sequences without allele-specific primers picked: <b>$snp_failed_seqs</b><p>";
            print HTML_REPORT_FILE "Total allele-specific primers picked: <b>$snp_total_primers</b><p>";
        }
        elsif ($primer_type eq '9') {
        }
        elsif ($primer_type eq '10') {
            print HTML_REPORT_FILE "Number of sequences with sucessful sequencing primers: <b>$snp_sucess_seqs</b><p>";
            print HTML_REPORT_FILE "Number of sequences without sequencing primers picked: <b>$snp_failed_seqs</b><p>";
            print HTML_REPORT_FILE "Total sequencing primers picked: <b>$snp_total_primers</b><p>";
        }
        else {
            print HTML_REPORT_FILE "Number of sequences with sucessful SNP primers: <b>$snp_sucess_seqs</b><p>";
            print HTML_REPORT_FILE "Number of sequences without SNP primers picked: <b>$snp_failed_seqs</b><p>";
            print HTML_REPORT_FILE "Total SNP primers picked: <b>$snp_total_primers</b><p>";
        }
    }


    my $used_time = (time()- $start_time);
 
    print HTML_REPORT_FILE "Used time: <b>$used_time</b> seconds.";
    print HTML_REPORT_FILE "</td></tr></tbody></table>";
    print HTML_REPORT_FILE "$wrapup\n";


    print HTML_TABLE_FILE "</tbody></table>\n";
    print HTML_TABLE_FILE "$wrapup\n";
    close(TABLE_FILE);
    close(HTML_TABLE_FILE);
    close(HTML_REPORT_FILE);
    
    if ($primer_type eq '2') {
        print SSR_HTML "</tbody></table>\n";
        close(SSR_TXT);
        close(SSR_HTML);
    }

}


# print generic primers or flanking primers for SNPs/allele/SSRs
sub print_primer_results {
    my ($included, $excluded, $target, $seq, $primer_type, $primer3output, $print_input, $sequence_id, $names, $input, $motif, $ssr) = @_;

    # output the result for each sequence
    my $filename;
    if ($primer_type eq '2'){ # SSR primer
     	$filename = $unique_id."_".$motif.".html";
    } else {
     	$filename = $unique_id.".html";
    }
    $filename = $user_id.$filename;

    my $outfilename = ">$user_result_dir/$filename";
    if(!open (NEW, $outfilename)) {
        print "<p>can\'t create the file $outfilename. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777 </p>";
    }

    print NEW $query->start_html(
        -title => "Primers of $sequence_id",
        -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
    );
    print NEW "<div class='mainpanel'>\n";
    print NEW "<h3 align=left>Sequence ID: $sequence_id</h3>\n";
    print NEW "<table class='standard'><tbody><tr><td>\n";
    
    my $repeats;

    if ($primer_type eq '2') {
        $repeats = length($ssr)/length($motif);
        print NEW "<table><tbody><tr><td><b>Motif:</b></td><td>(".uc($motif).")".$repeats."</td></tr>";
	print NEW "<tr><td><b>Motif length:</b></td><td>".length($motif)."</td></tr>\n";
	print NEW "<tr><td><b>SSR:</b></td><td>".uc($ssr)."</td></tr>";
	print NEW "<tr><td><b>SSR length:</b></td><td>".length($ssr)."</td></tr></tbody></table>\n";
    }

    print HTML_REPORT_FILE "<p><b>Sequence Index: $total_seqs</b><br>\n";
    print HTML_REPORT_FILE "<b>Sequence ID:    <a target=_blank href=\"$filename\">$sequence_id</a></b></p>\n";
    
    #print motif info
    if ($primer_type eq '2') {
        print HTML_REPORT_FILE "<table><tbody><tr><td><b>Motif:</b></td><td>(".uc($motif).")".$repeats."</td></tr>";
	print HTML_REPORT_FILE "<tr><td><b>Motif length:</b></td><td>".length($motif)."</td></tr>\n";
	print HTML_REPORT_FILE "<tr><td><b>SSR:</b></td><td>".uc($ssr)."</td></tr>";
	print HTML_REPORT_FILE "<tr><td><b>SSR length:</b></td><td>".length($ssr)."</td></tr></tbody></table>\n";
    }

    my $file_handle = \*NEW;    
    my $primer_list = $primer3output->get_primer_list();
    
    &print_generic_primers($file_handle, $primer_list, $primer_type, $filename, $sequence_id, $motif, $ssr);
    &print_primer_design_view ($file_handle, $included, $excluded, $target, $seq, $primer_list); # no snp/allele-specific primer list
    
    print NEW "</td></tr></tbody></table>\n";
    print NEW "$wrapup\n";
    close NEW;

}

# print flanking primers and SBE or allele-specific primers
sub print_primer_results_for_flanking_and_SBE_alleles {
    my ($included, $excluded, $target, $seq, $primer_type, $primer3output,
        $print_input, $sequence_id, $names, $input, $snp_f_all, $snp_r_all) = @_;

    my $filename = $unique_id.".html";
    
    $filename = $user_id.$filename;

    my $outfilename = ">$user_result_dir/$filename";
    if(!open (NEW, $outfilename)) {
        print "can\'t create the file $outfilename. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777 ";
    }

    print NEW $query->start_html(
        -title => "Primers of $sequence_id",
        -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
    );
    print NEW "<div class='mainpanel'>\n";
    print NEW "<h3 align=left>Sequence ID: $sequence_id</h3>\n";
    
    print NEW "<table class='standard'><tbody><tr><td>\n";
    

    print HTML_REPORT_FILE "<p><b>Sequence Index: $total_seqs</b><br>\n";
    print HTML_REPORT_FILE "<b>Sequence ID:    <a target=_blank href=\"$filename\">$sequence_id</a></b></p>\n";
    
    my $file_handle = \*NEW;
    
    my $primer_list = $primer3output->get_primer_list();
    &print_generic_primers($file_handle, $primer_list, $primer_type, $filename, $sequence_id);
    &print_snp_primers($file_handle, $primer_type, $filename, $sequence_id, $snp_f_all, $snp_r_all);
    
    &print_primer_design_view($file_handle, $included, $excluded, $target, $seq, $primer_list, $snp_f_all, $snp_r_all);
    
    print NEW "</td></tr></tbody></table>\n";
    print NEW "$wrapup\n";
    close NEW;
}

# print SBE or allele-specific primers only
sub print_primer_results_for_SBE_alleles {
    my ($included, $excluded, $target, $seq, $primer_type, $sequence_id, $snp_f_all, $snp_r_all, $snp, $pos) = @_;
    
    my $filename = $unique_id.".html";
    $filename = $user_id.$filename;

    my $outfilename = ">$user_result_dir/$filename";
    if(!open (NEW, $outfilename)) {
        print "can\'t create the file $outfilename. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777 ";
    }

    print NEW $query->start_html(
        -title => "Primers of $sequence_id",
        -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
    );
    print NEW "<div class='mainpanel'>\n";
    print NEW "<h3 align=left>Sequence ID: $sequence_id</h3>\n";
    
    print NEW "<table class='standard'><tbody><tr><td>\n";
    

    print HTML_REPORT_FILE "<p><b>Sequence Index: $total_seqs</b><br>\n";
    print HTML_REPORT_FILE "<b>Sequence ID:    <a target=_blank href=\"$filename\">$sequence_id</a></b></p>\n";
    
    my $file_handle = \*NEW;
    
    &print_snp_primers($file_handle, $primer_type, $filename, $sequence_id, $snp_f_all, $snp_r_all, $snp, $pos);


    my $primer_list = undef; # no flanking primer list
    &print_primer_design_view ($file_handle, $included, $excluded, $target, $seq, $primer_list, $snp_f_all, $snp_r_all);
    
    print NEW "</td></tr></tbody></table>\n";
    print NEW "$wrapup\n";
    close NEW;
}

sub print_snp_primers {
    my ($file_handle, $primer_type, $filename, $sequence_id, $results_f_all, $results_r_all) = @_;

    my $count = 0;
                           
    my $has_primers = 0;
    for (my $i = 0; $i < scalar(@$results_f_all); $i++) {
        if (@$results_f_all[$i]) {
            $has_primers = 1;
            last;
        }
    }
    for (my $i = 0; $i < scalar(@$results_r_all) && !$has_primers; $i++) {
        if (@$results_r_all[$i]) {
            $has_primers = 1;
            last;
        }
    }
                            
    if ($has_primers) {
        if ($primer_type eq '5' || $primer_type eq '6') {
            print HTML_REPORT_FILE "<p>SBE primers:</p>\n";
            print $file_handle "<p><b>SBE primers:</b></p>\n";
        } elsif ($primer_type eq '7' || $primer_type eq '8') {
            print HTML_REPORT_FILE "<p>Allele-specific primers:</p>\n";
            print $file_handle "<p><b>Allele-specific primers:</b></p>\n";
        } elsif ($primer_type eq '9') {
            print HTML_REPORT_FILE "<p>Inner primers for tetra-primer ARMS PCR:</p>\n";
            print $file_handle "<p><b>Inner primers for tetra-primer ARMS PCR:</b></p>\n";
        }
        $snp_sucess_seqs++;
    } else {
        $snp_failed_seqs++;
        if ($primer_type eq '5' || $primer_type eq '6') {
            print HTML_REPORT_FILE "<p><font color=red>No SBE primers found.</font></p>\n";
            print $file_handle "<p><font color=red>No SBE primers found.</font></p>\n";
        } elsif ($primer_type eq '7' || $primer_type eq '8') {
            print HTML_REPORT_FILE "<p><font color=red>No allele-specific primers found.</font></p>\n";
            print $file_handle "<p><font color=red>No allele-specific primers found.</font></p>\n";
        } elsif ($primer_type eq '7' || $primer_type eq '8') {
            print HTML_REPORT_FILE "<p><font color=red>No allele-specific primers found.</font></p>\n";
            print $file_handle "<p><font color=red>No allele-specific primers found.</font></p>\n";
        } elsif ($primer_type eq '9') {
            print HTML_REPORT_FILE "<p><font color=red>No inner primers for tetra-primer ARMS PCR found.</font></p>\n";
            print $file_handle "<p><font color=red>No inner primers for tetra-primer ARMS PCR found.</font></p>\n";
        }
        return;
    }
 
    my $primer_count = 0;
                         
    my $snp_header = "<table><thead><tr><th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Score</th><th> SNP</th><th>Pos</th><th>Primer Seq</th></tr></thead>\n";
    if ($primer_type eq '9') {
        $snp_header = "<table><thead><tr><th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th>" .
            "<th>GC%</th><th>Any compl</th><th>3' compl</th><th>Score</th><th> SNP</th><th>Pos</th><th>Primer Seq</th>" .
            "<th>Product Size</th><th>Seq Size</th><th>Included Size</th><th>Pair_any</th><th>Pair 3'</th>".
            "</tr></thead>\n";
    }
    print HTML_REPORT_FILE $snp_header;
    print HTML_REPORT_FILE "<tbody>\n";
    
    print $file_handle $snp_header;
    print $file_handle "<tbody>\n";
    
    
    if ($results_f_all){
    for (my $i = 0; $i < scalar(@$results_f_all); $i++) {
        if (@$results_f_all[$i] ) {
            my @primer_results = split(/\s+/, @$results_f_all[$i]);
            $primer_count++;
            $snp_total_primers++;
            my $primer_type_str = "Allele-specific";
            if ($primer_type == 5 || $primer_type == 6) {
                $primer_type_str = "SBE";
            }
            elsif ($primer_type == 9) {
                $primer_type_str = "Inner primer";
            }
            my $direction = "FORWARD";
            $snp_total_primers++;
            $table_row_count++;
            &print_single_snp_primer($file_handle, \@primer_results, $direction, $sequence_id, $primer_type_str, $primer_count,$table_row_count, $primer_type, $filename);
        }
    }
    }
   
    if ($results_r_all){
    for (my $i = 0; $i < scalar(@$results_r_all); $i++) {
        if (@$results_r_all[$i] ) {
            my @primer_results = split(/\s+/, @$results_r_all[$i]);
            $primer_count++;
            my $direction = "REVERSE";
            $snp_total_primers++;
            my $primer_type_str = "Allele-specific";
            if ($primer_type == 5 || $primer_type == 6) {
                $primer_type_str = "SBE";
            }
            elsif ($primer_type == 9) {
                $primer_type_str = "Inner primer";
            }
            $snp_total_primers++;
            $table_row_count++;
            &print_single_snp_primer($file_handle, \@primer_results, $direction, $sequence_id, $primer_type_str, $primer_count,$table_row_count,$primer_type, $filename);
        }
    }
    }
    if ($has_primers) {
        print HTML_REPORT_FILE "</tbody></table>\n";
        print $file_handle "</tbody></table>\n";
    }
    
}

sub print_SSR_screening_results {
    my ($sequence_id, $motifs, $ssrs, $targets) = @_;

    my $ssr_count = 0;
    
 
    for (my $i = 0; $i < @$motifs; $i++) {    
        $ssr_table_row_count++;
        $ssr_count++;
        
        my ($ssr_start, $ssr_len) = split(/,/, @$targets[$i]);
        my $ssr_end = $ssr_start + $ssr_len;
        
        if ($ssr_table_row_count % 2 == 0) {
            print SSR_HTML "<tr class='even'>\n";
        }
        else {
            print SSR_HTML "<tr class='odd'>\n";
        }

        print SSR_TXT $ssr_table_row_count."\t".$sequence_id."\t";
        print SSR_TXT "$ssr_count\t";

        print SSR_HTML "<td>$ssr_table_row_count</td><td>$sequence_id</td>";
        print SSR_HTML "<td>"."$ssr_count</td>";
    
        print SSR_TXT "$ssr_start\t$ssr_end\t";
        print SSR_TXT uc(@$motifs[$i])."\t";
        print SSR_TXT length(@$motifs[$i])."\t";
        print SSR_TXT @$ssrs[$i]."\t";
        print SSR_TXT length(@$ssrs[$i])."\n";
        
        print SSR_HTML "<td>$ssr_start</td><td>$ssr_end</td>";
        print SSR_HTML "<td>".uc(@$motifs[$i])."</td>";
        print SSR_HTML "<td>".length(@$motifs[$i])."</td>";
        print SSR_HTML "<td>".@$ssrs[$i]."</td>";
        print SSR_HTML "<td>".length(@$ssrs[$i])."</td>";
        print SSR_HTML "</tr>\n";
    }
}

sub print_generic_primers {
    my ($file_handle, $primer_list, $primer_type, $filename, $sequence_id, $motif, $ssr) = @_;

    my $count = 0;
    my $list_size = scalar(@$primer_list);
    if ($list_size == 0) {
        my $str = "";
        if ($primer_type eq '1') {
            $str = "No generic primers found!";
        }
        elsif ($primer_type eq '2') {
            $str = "No SSR-flanking primers found!";
        }
        elsif ($primer_type eq '3') {
            $str = "No oligo primer found!";
        }
        
        elsif ($primer_type eq '4' || $primer_type eq '5') {
            $str = "No SNP/allele flanking primers found!";
        }
        
        elsif ($primer_type eq '7') {
            $str = "No allele flanking primers found!";
        }
        
        elsif ($primer_type eq '9') {
            $str = "No outer primers found!";
        }
        print HTML_REPORT_FILE "<font color=red>$str</font><br>\n";
        print $file_handle "<font color=red>$str</font><br>\n";
        
	$failed_seqs++ if $primer_type ne '2' && $primer_type ne '9' ;
	return;
    } else {
        $sucess_seqs++ if $primer_type ne '2' && $primer_type ne '9';
    }
       
    my @headers = (
        "Generic primers",
        "SSR flanking primers",
        "Oligo primers",
        "SNP/Allele flanking primers",
        "SNP/Allele flanking primers",
        "",
        "Allele flanking primers",
        "",
        "Outer primers for tetra-primer ARMS PCR");
    print HTML_REPORT_FILE "<p>$headers[$primer_type-1]:</p>";
    print $file_handle "<p><b>$headers[$primer_type-1]:</b></p>";
    
    print HTML_REPORT_FILE "<table><thead><tr>\n";
    print $file_handle "<table><thead><tr>\n";
    if ($primer_type ne '3') {
	print  HTML_REPORT_FILE "<th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Primer Seq</th><th>Product Size</th><th>Seq Size</th><th>Included Size</th><th>Pair_any</th><th>Pair_3'</th>\n";
	print  $file_handle "<th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Primer Seq</th><th>Product Size</th><th>Seq Size</th><th>Included Size</th><th>Pair_any</th><th>Pair_3'</th>\n";
    } else {
	print  HTML_REPORT_FILE "<th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Primer Seq</th><th>Seq Size</th><th>Included Size</th>\n";
	print  $file_handle "<th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3'compl</th><th>Primer Seq</th><th>Seq Size</th><th>Included Size</th>\n";
    }
    print HTML_REPORT_FILE "</tr></thead><tbody>\n";
    print $file_handle "</tr></thead><tbody>\n";
    
    
    my $primer_type_str = "";
    if ($primer_type eq '1') { #Generic
        $primer_type_str = "Generic";
    }
    elsif ($primer_type eq '2') { #SSR
        $primer_type_str = "SSR";
    }
    elsif ($primer_type eq '3') { #Oligo
        $primer_type_str = "Oligo";
    }
    elsif ($primer_type eq '4' || $primer_type eq '5' || $primer_type eq '6') { # SBE flanking
        $primer_type_str = "SNP flanking";
    }
    elsif ($primer_type eq '7' || $primer_type eq '8') { # allele flanking
        $primer_type_str = "Allele flanking";
    }
    elsif ($primer_type eq '9') { # allele flanking
        $primer_type_str = "Outer primer";
    }
    
    
    foreach my $primer_pair (@$primer_list) {
        $count++;
        $total_primer_pairs++ if ($primer_type ne '9');
        $table_row_count++;
        
        print TABLE_FILE $table_row_count."\t".$sequence_id."\t";
        print TABLE_FILE "$count\t";
        
        if ($table_row_count % 2 == 0) {
            print HTML_TABLE_FILE "<tr class='even'>\n";
        }
        else {
            print HTML_TABLE_FILE "<tr class='odd'>\n";
        }

        print HTML_TABLE_FILE "<td>$table_row_count</td><td><a target=_blank href=\"$filename\">$sequence_id</a></td>";
        print HTML_TABLE_FILE "<td>"."$count</td>";

    
	
	if ($primer_type ne '3') {
	    #print left primer
	    my $left_primer_obj = $primer_pair->left_primer();
	    my @row_left = (
		   $left_primer_obj->direction(),
		   $left_primer_obj->start(),
		   $left_primer_obj->length(),
		   $left_primer_obj->tm(),
		   $left_primer_obj->gc(),
	           $left_primer_obj->any_complementarity(),
	           $left_primer_obj->end_complementarity(),
	           uc($left_primer_obj->sequence()),
		   $primer_pair->product_size(),
		   $primer_pair->seq_size(),
		   $primer_pair->included_size(),
		   $primer_pair->pair_any_complementarity(),
		   $primer_pair->pair_end_complementarity());
                    
	    print HTML_REPORT_FILE "<tr><td>$count</td>";
	    print $file_handle "<tr><td>$count</td>";
            print TABLE_FILE "$primer_type_str\t";
            print HTML_TABLE_FILE "<td>$primer_type_str</td>";
	    for (my $i = 0; $i < scalar(@row_left); $i++) {
                print HTML_REPORT_FILE "<td>$row_left[$i]</td>";
                print $file_handle "<td>$row_left[$i]</td>";
                if ($i == 0) {
                    print TABLE_FILE $row_left[$i]."\t";
		    print HTML_TABLE_FILE "<td>".$row_left[$i]."</td>";
		    
		} else {
		    if ($i <= 7) {
		        print TABLE_FILE $row_left[$i]."\t";
		        print HTML_TABLE_FILE "<td>".$row_left[$i]."</td>";
                        if ($i == 6 && ($primer_type eq '5' || $primer_type eq '7' || $primer_type eq '9')) {
                            print TABLE_FILE "\t\t\t";
                            print HTML_TABLE_FILE "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>";
                        }
		    }
		}
	    }

	    for (my $i = 8; $i < scalar(@row_left); $i++) {
	        print TABLE_FILE $row_left[$i]."\t";
	        print HTML_TABLE_FILE "<td>".$row_left[$i]."</td>";
	    }

	    if ($primer_type eq '2') { #SSR
   
		print TABLE_FILE uc($motif)."\t";
		print TABLE_FILE length($motif)."\t";
		print TABLE_FILE "$ssr\t";
		print TABLE_FILE length($ssr)."\t";

		print HTML_TABLE_FILE "<td>".uc($motif)."</td>";
		print HTML_TABLE_FILE "<td>".length($motif)."</td>";
		print HTML_TABLE_FILE "<td>"."$ssr</td>";
		print HTML_TABLE_FILE "<td>".length($ssr)."</td>";
	    }	

            print HTML_REPORT_FILE "</tr>\n";
	    print $file_handle "</tr>\n";
            print TABLE_FILE "\n";
	    print HTML_TABLE_FILE "</tr>\n";

	    #print right primer
            print HTML_REPORT_FILE "<tr><td>&nbsp;</td>";
            print $file_handle "<tr><td>&nbsp;</td>";

            $table_row_count++;
            print TABLE_FILE $table_row_count."\t".$sequence_id."\t";
            print TABLE_FILE "$count\t";
            
            if ($table_row_count % 2 == 0) {
                print HTML_TABLE_FILE "<tr class='even'>\n";
            }
            else {
                print HTML_TABLE_FILE "<tr class='odd'>\n";
            }
            print HTML_TABLE_FILE "<td>$table_row_count</td><td><a target=_blank href=\"$filename\">$sequence_id</a></td>";
            print HTML_TABLE_FILE "<td>"."$count</td>";
	
            my $right_primer_obj = $primer_pair->right_primer();
	    my @row_right = (
		   $right_primer_obj->direction(),
		   $right_primer_obj->start(),
		   $right_primer_obj->length(),
		   $right_primer_obj->tm(),
		   $right_primer_obj->gc(),
	           $right_primer_obj->any_complementarity(),
	           $right_primer_obj->end_complementarity(),
	           uc($right_primer_obj->sequence()));
        
  	    print TABLE_FILE "$primer_type_str\t";
            print HTML_TABLE_FILE "<td>$primer_type_str</td>";
	
            for (my $i = 0; $i < scalar(@row_right); $i++) {
                print HTML_REPORT_FILE "<td>$row_right[$i]</td>";
                print $file_handle "<td>$row_right[$i]</td>";
	        print TABLE_FILE $row_right[$i]."\t";
	        print HTML_TABLE_FILE "<td>".$row_right[$i]."</td>";
		if ($i == 6 && ($primer_type eq '5' || $primer_type eq '7' || $primer_type eq '9')) {
                    print TABLE_FILE "\t\t\t";
                    print HTML_TABLE_FILE "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>";
                }
            }
            if ($primer_type eq '1' || $primer_type eq '4' || $primer_type eq '5' || $primer_type eq '7' || $primer_type eq '9') {
                print TABLE_FILE "\t\t\t\t\t";
                print HTML_TABLE_FILE "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>";
            }
            elsif ($primer_type eq '2') {
                print TABLE_FILE "\t\t\t\t\t\t\t\t\t";
                print HTML_TABLE_FILE "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>";
            }
            
	    print HTML_REPORT_FILE "</tr>\n";
	    print $file_handle "</tr>\n";
            print TABLE_FILE "\n";
	    print HTML_TABLE_FILE "</tr>\n";


	} elsif ($primer_type eq '3') { #oligo
	    my $oligo_primer_obj = $primer_pair->oligo_primer();
	    my @row = (
		   $oligo_primer_obj->direction(),
		   $oligo_primer_obj->start(),
		   $oligo_primer_obj->length(),
		   $oligo_primer_obj->tm(),
		   $oligo_primer_obj->gc(),
	           $oligo_primer_obj->any_complementarity(),
	           $oligo_primer_obj->end_complementarity(),
	           uc($oligo_primer_obj->sequence()),
		   $primer_pair->seq_size(),
		   $primer_pair->included_size());
 
	    print HTML_REPORT_FILE "<tr><td>$count</td>";
	    print $file_handle "<tr><td>$count</td>";
	    for (my $i = 0; $i < scalar(@row); $i++) {
	        print HTML_REPORT_FILE "<td>$row[$i]</td>";
	        print $file_handle "<td>$row[$i]</td>";
	        if ($i <=9 ) {
		    print TABLE_FILE $row[$i]."\t";
		    print HTML_TABLE_FILE "<td>".$row[$i]."</td>";
		}
	    }
	    print HTML_REPORT_FILE "</tr>\n";
	    print $file_handle "</tr>\n";
	    print TABLE_FILE "\n";
	    print HTML_TABLE_FILE "</tr>";
	    
	}
	     
    }
    print HTML_REPORT_FILE "</tbody></table>\n";
    print $file_handle "</tbody></table>\n";
}            


sub print_single_snp_primer_for_summary {
    my ($file_handle, $primer_results, $direction, $primer_count) = @_;
    print HTML_REPORT_FILE "<tr>";
    print HTML_REPORT_FILE "<td>$primer_count</td>";
    print HTML_REPORT_FILE "<td>$direction</td>";
    print HTML_REPORT_FILE "<td>".(@$primer_results[3]+1)."</td>";
    print HTML_REPORT_FILE "<td>".@$primer_results[4]."</td>";
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[5]*100 + 0.5)/100)."</td>";
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[6]*100 + 0.5)/100)."</td>";
    print HTML_REPORT_FILE "<td>".@$primer_results[7]."</td>";                       # complementarity  
    print HTML_REPORT_FILE "<td>".@$primer_results[8]."</td>";                       # 3' complementarity
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[9]*100 + 0.5)/100)."</td>";;  # score
    print HTML_REPORT_FILE "<td>".@$primer_results[1]."</td>";;
    print HTML_REPORT_FILE "<td>".@$primer_results[2]."</td>";;
    print HTML_REPORT_FILE "<td>".@$primer_results[10]."</td>";
    if ($primer_type eq '9' && scalar(@$primer_results) >= 16 ) {
        print HTML_REPORT_FILE "<td>".@$primer_results[11]."</td>";
        print HTML_REPORT_FILE "<td>".@$primer_results[12]."</td>";
        print HTML_REPORT_FILE "<td>".@$primer_results[13]."</td>";
        print HTML_REPORT_FILE "<td>".@$primer_results[14]."</td>";
        print HTML_REPORT_FILE "<td>".@$primer_results[15]."</td>";
    }
    
    print HTML_REPORT_FILE "</tr>";

    print $file_handle "<tr>";
    print $file_handle "<td>$primer_count</td>";
    print $file_handle "<td>$direction</td>";
    print $file_handle "<td>".(@$primer_results[3]+1)."</td>";
    print $file_handle "<td>".@$primer_results[4]."</td>";
    print $file_handle "<td>".(int(@$primer_results[5]*100 + 0.5)/100)."</td>";
    print $file_handle "<td>".(int(@$primer_results[6]*100 + 0.5)/100)."</td>";
    print $file_handle "<td>".@$primer_results[7]."</td>";                       # complementarity  
    print $file_handle "<td>".@$primer_results[8]."</td>";                       # 3' complementarity
    print $file_handle "<td>".(int(@$primer_results[9]*100 + 0.5)/100)."</td>";;  # score
    print $file_handle "<td>".@$primer_results[1]."</td>";;
    print $file_handle "<td>".@$primer_results[2]."</td>";;
    print $file_handle "<td>".@$primer_results[10]."</td>";
    if ($primer_type eq '9' && scalar(@$primer_results) >= 16) {
        print $file_handle "<td>".@$primer_results[11]."</td>";
        print $file_handle "<td>".@$primer_results[12]."</td>";
        print $file_handle "<td>".@$primer_results[13]."</td>";
        print $file_handle "<td>".@$primer_results[14]."</td>";
        print $file_handle "<td>".@$primer_results[15]."</td>";
    }
    print $file_handle "</tr>";


}

sub print_single_sequencing_primer {
    my ($file_handle, $primer_results, $direction, $sequence_id, $primer_type_str, $count, $index, $primer_type, $filename, $seq_length) = @_;

    # print summary
    print HTML_REPORT_FILE "<tr>";
    print HTML_REPORT_FILE "<td>$count</td>";
    print HTML_REPORT_FILE "<td>$direction</td>";
    if ($SEQUENCING_DIRECTION eq "toward5") {  # reverse
        print HTML_REPORT_FILE "<td>".(@$primer_results[1] - $seq_length-2)."</td>";
    } else {
        print HTML_REPORT_FILE "<td>".(@$primer_results[1]+1)."</td>";   # index start from 1
    }
    print HTML_REPORT_FILE "<td>".@$primer_results[2]."</td>";
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[3]*100 + 0.5)/100)."</td>";
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[4]*100 + 0.5)/100)."</td>";
    print HTML_REPORT_FILE "<td>".@$primer_results[5]."</td>";                       # complementarity  
    print HTML_REPORT_FILE "<td>".@$primer_results[6]."</td>";                       # 3' complementarity
    print HTML_REPORT_FILE "<td>".(int(@$primer_results[7]*100 + 0.5)/100)."</td>";;  # score
    print HTML_REPORT_FILE "<td>".@$primer_results[8]."</td>";
    print HTML_REPORT_FILE "</tr>";

    print $file_handle "<tr>";
    print $file_handle "<td>$count</td>";
    print $file_handle "<td>$direction</td>";
    if ($SEQUENCING_DIRECTION eq "toward5") {  # reverse
        print $file_handle "<td>".(@$primer_results[1] - $seq_length-2)."</td>";                       # start
    } else {
        print $file_handle "<td>".(@$primer_results[1]+1)."</td>";                       # start
    }
    print $file_handle "<td>".@$primer_results[2]."</td>";                       # length
    print $file_handle "<td>".(int(@$primer_results[3]*100 + 0.5)/100)."</td>";  # Tm
    print $file_handle "<td>".(int(@$primer_results[4]*100 + 0.5)/100)."</td>";  # GC %
    print $file_handle "<td>".@$primer_results[5]."</td>";                       # complementarity  
    print $file_handle "<td>".@$primer_results[6]."</td>";                       # 3' complementarity
    print $file_handle "<td>".(int(@$primer_results[7]*100 + 0.5)/100)."</td>";;  # score
    print $file_handle "<td>".@$primer_results[8]."</td>";
    print $file_handle "</tr>";

    # print table
    print TABLE_FILE $index."\t";
    print TABLE_FILE $sequence_id."\t";
    print TABLE_FILE $count."\t";
    print TABLE_FILE $primer_type_str."\t";
    print TABLE_FILE $direction."\t";
    if ($SEQUENCING_DIRECTION eq "toward5") {  # reverse
        print TABLE_FILE (@$primer_results[1] - $seq_length-2)."\t"; # start
    } else {
        print TABLE_FILE (@$primer_results[1] + 1)."\t";                   # start
    }
    print TABLE_FILE @$primer_results[2]."\t";                       # length
    print TABLE_FILE (int(@$primer_results[3]*100 + 0.5)/100)."\t";  # Tm
    print TABLE_FILE (int(@$primer_results[4]*100 + 0.5)/100)."\t";  # GC %
    print TABLE_FILE @$primer_results[5]."\t";                       # complementarity  
    print TABLE_FILE @$primer_results[6]."\t";                       # 3' complementarity
    print TABLE_FILE (int(@$primer_results[7]*100 + 0.5)/100)."\t";  # score
    print TABLE_FILE @$primer_results[8]."\n";

    
    if ($table_row_count % 2 == 0) {
        print HTML_TABLE_FILE "<tr class='even'>\n";
    }
    else {
        print HTML_TABLE_FILE "<tr class='odd'>\n";
    }
    print HTML_TABLE_FILE "<td>$index</td>";
    print HTML_TABLE_FILE "<td><a target=_blank href=\"$filename\">$sequence_id</a></td>";
    print HTML_TABLE_FILE "<td>$count</td>";
    print HTML_TABLE_FILE "<td>$primer_type_str</td>";
    print HTML_TABLE_FILE "<td>$direction</td>";
    if ($SEQUENCING_DIRECTION eq "toward5") {  # reverse
        print HTML_TABLE_FILE "<td>".(@$primer_results[1] - $seq_length-2)."</td>"; # start
    } else {
        print HTML_TABLE_FILE "<td>".(@$primer_results[1]+1)."</td>";                   # start
    }
    print HTML_TABLE_FILE "<td>@$primer_results[2]</td>";
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[3]*100 + 0.5)/100)."</td>";
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[4]*100 + 0.5)/100)."</td>";
    print HTML_TABLE_FILE "<td>@$primer_results[5]</td>";                       # complementarity  
    print HTML_TABLE_FILE "<td>@$primer_results[6]</td>";                       # 3' complementarity
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[7]*100 + 0.5)/100)."</td>";  # score
    print HTML_TABLE_FILE "<td>@$primer_results[8]</td>";
    print HTML_TABLE_FILE "</tr>\n";
}



sub print_single_snp_primer {
    my ($file_handle, $primer_results, $direction, $sequence_id, $primer_type_str, $count, $index, $primer_type, $filename) = @_;

    &print_single_snp_primer_for_summary($file_handle, $primer_results, $direction, $count, $primer_type);
   
    print TABLE_FILE $index."\t";
    print TABLE_FILE $sequence_id."\t";
    print TABLE_FILE $count."\t";
    print TABLE_FILE $primer_type_str."\t";
    print TABLE_FILE $direction."\t";
    print TABLE_FILE (@$primer_results[3]+1)."\t";
    print TABLE_FILE @$primer_results[4]."\t";
    print TABLE_FILE (int(@$primer_results[5]*100 + 0.5)/100)."\t";
    print TABLE_FILE (int(@$primer_results[6]*100 + 0.5)/100)."\t";
    print TABLE_FILE @$primer_results[7]."\t";                       # complementarity  
    print TABLE_FILE @$primer_results[8]."\t";                       # 3' complementarity
    print TABLE_FILE (int(@$primer_results[9]*100 + 0.5)/100)."\t";  # score
    print TABLE_FILE @$primer_results[1]."\t";
    print TABLE_FILE @$primer_results[2]."\t";
    print TABLE_FILE @$primer_results[10]."\t";
    if ($primer_type eq '9' && scalar(@$primer_results) >= 16 ) {
        print TABLE_FILE @$primer_results[11]."\t";
        print TABLE_FILE @$primer_results[12]."\t";
        print TABLE_FILE @$primer_results[13]."\t";
        print TABLE_FILE @$primer_results[14]."\t";
        print TABLE_FILE @$primer_results[15];
    }
    print TABLE_FILE "\n";
    
    if ($table_row_count % 2 == 0) {
        print HTML_TABLE_FILE "<tr class='even'>\n";
    }
    else {
        print HTML_TABLE_FILE "<tr class='odd'>\n";
    }
    print HTML_TABLE_FILE "<td>$index</td>";
    print HTML_TABLE_FILE "<td><a target=_blank href=\"$filename\">$sequence_id</a></td>";
    print HTML_TABLE_FILE "<td>$count</td>";
    print HTML_TABLE_FILE "<td>$primer_type_str</td>";
    print HTML_TABLE_FILE "<td>$direction</td>";
    print HTML_TABLE_FILE "<td>".(@$primer_results[3]+1)."</td>";
    print HTML_TABLE_FILE "<td>@$primer_results[4]</td>";
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[5]*100 + 0.5)/100)."</td>";
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[6]*100 + 0.5)/100)."</td>";
    print HTML_TABLE_FILE "<td>@$primer_results[7]</td>";                       # complementarity  
    print HTML_TABLE_FILE "<td>@$primer_results[8]</td>";                       # 3' complementarity
    print HTML_TABLE_FILE "<td>".(int(@$primer_results[9]*100 + 0.5)/100)."</td>";  # score
    print HTML_TABLE_FILE "<td>@$primer_results[1]</td>";
    print HTML_TABLE_FILE "<td>@$primer_results[2]</td>";
    print HTML_TABLE_FILE "<td>@$primer_results[10]</td>";
    if ($primer_type eq '9' && scalar(@$primer_results) >= 16 ) {
        print HTML_TABLE_FILE "<td>".@$primer_results[11]."</td>";
        print HTML_TABLE_FILE "<td>".@$primer_results[12]."</td>";
        print HTML_TABLE_FILE "<td>".@$primer_results[13]."</td>";
        print HTML_TABLE_FILE "<td>".@$primer_results[14]."</td>";
        print HTML_TABLE_FILE "<td>".@$primer_results[15]."</td>";
    }
    print HTML_TABLE_FILE "</tr>\n";
}



sub print_primer_totablefile{
	my $str = shift;
    my @tokens = split(/\s+/, $str);
    my $start = 2;
    if ($primer_type eq '3') {
        $start = 1;
    }
    for (my $i = $start; $i < scalar(@tokens); $i++) {
    	print TABLE_FILE uc($tokens[$i])."\t";
    	print HTML_TABLE_FILE "<td>".uc($tokens[$i])."</td>";
    }
}

sub print_sequence_size_totablefile {
	my $str = shift;
    my @tokens = split(/\s*:\s*/, $str);
    chomp($tokens[1]);
    print TABLE_FILE $tokens[1]."\t";
    print HTML_TABLE_FILE "<td>".$tokens[1]."</td>";
}


sub print_include_region_size_totablefile{
    my $str = shift;
    my @tokens = split(/\s*:\s*/, $str);
    chomp($tokens[1]);
    print TABLE_FILE $tokens[1]."\t";
    print HTML_TABLE_FILE "<td>".$tokens[1]."</td>";
}

sub print_product_size_totablefile {
	my $str = shift;
    my @tokens = split(/\s*,\s*/, $str);
    for (my $i = 0; $i < scalar(@tokens); $i++) {
    	my @tokens1 = split(/\s*:\s*/, $tokens[$i]);
    	chomp($tokens1[1]);
    	print TABLE_FILE $tokens1[1]."\t";
    	print HTML_TABLE_FILE "<td>".$tokens1[1]."</td>";
    }
}


sub no_primers_found {
    return qq{
	</pre>
	<h2>No Acceptable Primers Were Found</h2>
	The statistics below should indicate why no acceptable
	primers were found.
	Try relaxing various parameters, including the
	self-complementarity parameters and max and min oligo melting
	temperatures.  For example, for very A-T-rich regions you might
	have to increase maximum primer size or decrease minimum melting
	temperature.

	<hr>

	<pre>
    }
}

sub add_start_len_list($$$) {
    my ($list_string, $list, $plus) = @_;
    my $sp = $list_string ? ' ' : '' ;
    for (@$list) {
	$list_string .= ($sp . ($_->[0] + $plus) . "," . $_->[1]);
	$sp = ' ';
    }
    return $list_string;
}

sub read_sequence_markup($@) {
    my ($s, @delims) = @_;
    
    # E.g. ['/','/'] would be ok in @delims, but
    # no two pairs in @delims may share a character.
    my @out = ();
    for (@delims) {
	push @out, read_sequence_markup_1_delim($s, $_, @delims);
    }
    @out;
}

sub read_sequence_markup_1_delim($$@) {
    my ($s,  $d, @delims) = @_;
    my ($d0, $d1) = @$d;
    my $other_delims = '';
    
    for (@delims) {
	next if $_->[0] eq $d0 and $_->[1] eq $d1;
	confess 'Programming error' if $_->[0] eq $d0;
	confess 'Programming error' if $_->[1] eq $d1;
	$other_delims .= '\\' . $_->[0] . '\\' . $_->[1];
    }
    if ($other_delims) {
	$s =~ s/[$other_delims]//g;
    }
    # $s now contains only the delimters of interest.
    my @s = split(//, $s);
    my ($c, $pos) = (0, 0);
    my @out = ();
    my $len;
    while (@s) {
	$c = shift(@s);
	next if ($c eq ' '); # Already used delimeters are set to ' '
	if ($c eq $d0) {
            $len = len_to_delim($d0, $d1, \@s);
            return undef if (!defined $len);
            push @out, [$pos, $len];
	} elsif ($c eq $d1) {
            # There is a closing delimiter with no opening
            # delimeter, an input error.
            $DO_NOT_PICK = 1;
            print "ERROR IN SEQUENCE: closing delimiter $d1 not preceded by $d0\n";
            return undef;
	} else {
            $pos++;
	}
    }
    
    return \@out;
}

sub len_to_delim($$$) {
    my ($d0, $d1, $s) = @_;
    my $i;
    my $len = 0;
    for $i (0..$#{$s}) {
		if ($s->[$i] eq $d0) {
			# ignore it;
		} elsif ($s->[$i] eq $d1) {
			$s->[$i] = ' ';
			return $len;
		} else { $len++ }
    }
    # There was no terminating delim;
    $DO_NOT_PICK = 1;
    print "ERROR IN SEQUENCE: closing delimiter $d1 did not follow $d0\n";
    return undef;
}


################  functions for SSR primer design  ###############################

# find SSRs in the included regions. SSRs will be ignored in the exlcluded regions
sub find_ssrs {
    my ($sequence, $included_region, $excluded_region) = @_;
    my @target_regions;
    my @prod_size_ranges;
    my @motifs;
    my @ssrs;

    my ($included_start, $included_end) = ();
    my @excluded_starts = ();
    my @excluded_ends = ();
    
    if ($included_region) {
       ($included_start, $included_end) = $included_region =~ /(\d+),(\d+)/;
       $included_end += $included_start;  
    }
    
    if ($excluded_region) {
        my @tmp = split(/\s+/, $excluded_region);
        for (my $i = 0; $i < @tmp; $i++) {
            ($excluded_starts[$i], $excluded_ends[$i]) = $tmp[$i] =~ /(\d+),(\d+)/;
            $excluded_ends[$i] += $excluded_starts[$i];
        }
    }

    #find the SSRs
    study($sequence);     # may be necessary?
    my $seqlength = length($sequence);
    my $ssr_number = 1;   #track multiple ssrs within a single sequence
    my %locations;        #track location of SSRs as detected
    my $count = 0;
    for(my $i=0; $i<scalar(@repeats); $i++){ #test each spec against sequence
        next if ($selected[$i] == 0);

        my $motiflength = $repeats[$i]->[0];
        my $minreps = $repeats[$i]->[1] - 1;
        my $regexp = "(([gatc]{$motiflength})\\2{$minreps,})";
        while ($sequence =~ /$regexp/ig){
            my $motif = lc($2);
            my $ssr = $1;
            next if &homopolymer($motif,$motiflength); #comment out this line to report monomers
            my $ssrlength = length($ssr);          #overall SSR length
     #       my $repeats = $ssrlength/$motiflength; #number of rep units
            my $end = pos($sequence);              #where SSR ends
            pos($sequence) = $end - $motiflength;  #see docs
            my $start = $end - $ssrlength + 1;     #where SSR starts
            
            if ($included_region && ($start < $included_start || $end >$included_end)) {
                next;
            }

            my $s = $start - 1;  # base index = 0
            my $e = $end - 1;    # base index = 0
            if ($excluded_region) {
                my $excluded = 0;
                for (my $j = 0; $j < @excluded_starts; $j++) {
                    if ($s <= $excluded_starts[$j]-1 && $e >= $excluded_ends[$j]-2  ||
                        $s >= $excluded_starts[$j]-1 && $e <= $excluded_ends[$j]-2 ||
                        $s >= $excluded_starts[$j]-1 && $s <= $excluded_ends[$j]-2 && $e >= $excluded_ends[$j]-2 ||
                        $s <= $excluded_starts[$j]-1 && $e >= $excluded_starts[$j]-1 && $e <= $excluded_ends[$j]-2) { 
                        $excluded = 1;
                        last;
                    }
                }
                next if $excluded == 1;
            }
            
       
            if (&novel($start, \%locations)) {   #count SSR only once
            	$target_regions[$count] = "$start,$ssrlength";
                my $s1 = $ssrlength+20;
                my $e1 = $ssrlength+40;
                $prod_size_ranges[$count] = "$s1-$e1";
		$motifs[$count] = $motif;
                $ssrs[$count] = $ssr;
                $count++;
            }

        }
    }
    return (\@target_regions, \@prod_size_ranges, \@motifs, \@ssrs);
}

sub homopolymer {
    #return true if motif is repeat of single nucleotide
    my ($motif,$motiflength) = @_;
    my ($reps) = $motiflength - 1;
    return 1 if ($motif =~ /([gatc])\1{$reps}/);
    return 0;
}

sub novel {
    my($position, $locationsref) = @_;
    if(defined $locationsref->{$position}) {
       return undef;
   } else {
       $locationsref->{$position} = 1;
       return 1;
   }
}

###########################################################################

sub find_snp_targets {
    my $sequence = shift;
    $sequence = uc($sequence);

    my $first_base_index = shift;

    my @target_regions = ();
    my @codes = ('S', 'W', 'R', 'Y', 'K', 'M', 'V', 'H', 'D', 'B');
    my $offset = 0;
    for (my $i = 0; $i < scalar(@codes); $i++) {
        $offset = 0;
    	my $index = index($sequence, $codes[$i], $offset);
        while ($index != -1) {
            push(@target_regions, ($index + $first_base_index).",1");
            $offset = $index + 1;
    	    $index = index($sequence, $codes[$i], $offset);
        }
    }
    return \@target_regions;
}

################  functions for sequencing primer design  ###############################

sub find_sequencing_primers {
    my ($included_region, $excluded_region, $seq) = @_;
    
    my $window_width = 20;
    my $selected_index = -1;

    my @candidate_strings = ();
    my @scores = ();
    my $direction = "Towards_3'";
    if ($SEQUENCING_DIRECTION eq "toward3") {
        $direction = "Towards_5'";
    }
    
    my ($included_start, $included_end) = (0, 0);
    
    if ($included_region) {
       ($included_start, $included_end) = $included_region =~ /(\d+),(\d+)/;
       $included_end += $included_start;
       $included_start--;  # to base index = 0
    }

    my ($candidates, $start_poses) = ();
    my $oligo_tm = new OligoTM();
    my $start = $included_start;
    $start = length($seq) - $included_end if ($SEQUENCING_DIRECTION eq "toward5" && $included_region);

 
    my $count = 0;
    for (my $i = $start; $i < length($seq); $i += $window_width) {
    
        # Step 1: pick all primer candidate from one direction user specified
        # reference array of candidate primers
        ($candidates, $start_poses) = ();
        if ($SEQUENCING_DIRECTION eq "toward3") {
            ($candidates, $start_poses) = &find_sequencing_candidate_primers($included_region, $excluded_region, $seq, $i, $window_width);
        } else { # toward5
            ($candidates, $start_poses) = &find_sequencing_candidate_primers($included_region, $excluded_region, &reverse_complement($seq), $i, $window_width);
        }
 
        # Step 2: calcuate Tm of all candidate and scores    
        for (my $i = 0; $i < scalar(@$candidates); $i++) { 
            my $tm1 = $oligo_tm->oligo_tm(@$candidates[$i], $SNP_PRIMER_DNA_CONC, $SNP_PRIMER_SALT_CONC);
            my $gc1 = &oligo_gc(@$candidates[$i]);
            my $repeats1 = &find_repeat_bases(@$candidates[$i]);
            my $ns1 = &find_ambiguity_ns(@$candidates[$i]);
            my ($self_comp1, $self_comp2) = &self_complementarity(@$candidates[$i]);
            my $score = calculate_primer_score(length(@$candidates[$i]),
                $tm1,$gc1, $repeats1, $ns1, $self_comp1, $self_comp2);
            if ($score > 60) {
                $scores[$count] = $score;
                $candidate_strings[$count] = 
                    "$direction\t".
                    "@$start_poses[$i]\t".
                    length(@$candidates[$i])."\t".
                    "$tm1\t$gc1\t".
                    "$self_comp1\t$self_comp2\t".
                    "$score\t@$candidates[$i]";
                $count++;
            }
        }

        last if (scalar(@scores) >= $PRIMER_NUM_RETURN);
    }

    # sorting the scores
    my $quick_sort = new QuickSort();
    my ($sorted_scores, $sort_index) = $quick_sort->sort(\@scores);
    
    
    my @results = ();
    for (my $i = 0; $i < $PRIMER_NUM_RETURN && scalar(@candidate_strings) > $i ; $i++) {
        $results[$i] = $candidate_strings[@$sort_index[scalar(@$sorted_scores)-$i -1]];
    }
  
    
    # return result arrays: direction, primer sequence, primer start position, length, Tm, score
    return (@results);
    
}


# find_sequencing_candidate_primers
# Find the candidate primers meeting the requirements
sub find_sequencing_candidate_primers {
    my ($included_region, $excluded_region, $seq, $start, $end) = @_;

    my @candidates = ();
    my @start_poses = ();
    my $count = 0;


    my ($included_start, $included_end) = ();
    my @excluded_starts = ();
    my @excluded_ends = ();
    
    if ($included_region) {
       ($included_start, $included_end) = $included_region =~ /(\d+),(\d+)/;
       $included_end += $included_start;  
    }
    
    if ($excluded_region) {
        my @tmp = split(/\s+/, $excluded_region);
        for (my $i = 0; $i < @tmp; $i++) {
            ($excluded_starts[$i], $excluded_ends[$i]) = $tmp[$i] =~ /(\d+),(\d+)/;
            $excluded_ends[$i] += $excluded_starts[$i];
        }
    }
       

    for (my $i = $start; $i < $end; $i++) {

        for (my $j = $SNP_PRIMER_MIN_SIZE; $j <= $SNP_PRIMER_MAX_SIZE; $j++) {
            last if ($i + $j + 1 > length($seq));
            
            my $s = $i;  # base index = 0
            my $e = $i+$j - 1;
            if ($SEQUENCING_DIRECTION eq "toward5") {
                my $tmp = $e;
                $e = length($seq)-$s-1;
                $s = length($seq)-$tmp-1;
            }

            if ($included_region) {
                next if ($s < $included_start-1 || $e >$included_end-2);
            }
            
            if ($excluded_region) {
                my $excluded = 0;
                for (my $k = 0; $k < @excluded_starts; $k++) {
                    # -1 or -2 to convert the index to the base index 0
                    if ($s <= $excluded_starts[$k]-1 &&  $e >= $excluded_ends[$k]-2  ||
                        $s >= $excluded_starts[$k]-1 && $e <= $excluded_ends[$k]-2 ||
                        $s >= $excluded_starts[$k]-1 && $s <= $excluded_ends[$k]-2 && $e >= $excluded_ends[$k]-2 ||
                        $s <= $excluded_starts[$k]-1 && $e >= $excluded_starts[$k]-1 && $e <= $excluded_ends[$k]-2) { 
                       $excluded = 1;
                       last;
                    }
                }
                next if $excluded == 1;
            }

            
            # sweep out some bad candidate
            my $candidate_sprimer = uc(substr($seq, $i, $j));  
            my $ns = find_ambiguity_ns($candidate_sprimer);
            next if ($ns > $SNP_PRIMER_MAX_N);
            
            $candidates[$count] = $candidate_sprimer;
            if ($SEQUENCING_DIRECTION eq "toward3") {
                $start_poses[$count] = $i;
            } else {
                $start_poses[$count] = length($seq) - $i + 1;
            }
            $count++;
        }
    }
    return (\@candidates, \@start_poses);
}

                          
sub print_primer_results_for_sequencing_primers {
    my ($included_region, $excluded_region, $seq, $primer_type, $sequence_id, $results, $target) = @_;
    
    my $filename = $unique_id.".html";
    $filename = $user_id.$filename;

    my $outfilename = ">$user_result_dir/$filename";
    if(!open (NEW, $outfilename)) {
        print "can\'t create the file $outfilename. Please check if your directory is correct. Probably you need to change the mode of the directory: chmod 777 ";
    }

    print NEW $query->start_html(
        -title => "Primers of $sequence_id",
        -style => { -src => "$HTDOC/$HTML_STYLE_FILE" },
    );
    print NEW "<div class='mainpanel'>\n";
    print NEW "<h3 align=left>Sequence ID: $sequence_id</h3>\n";
    
    print NEW "<table class='standard'><tbody><tr><td>\n";
    

    print "<p><b>Sequence Index: $total_seqs</b><br>\n";
    print "<b>Sequence ID:    <a target=_blank href=\"$filename\">$sequence_id</a></b></p>\n";
    
    my $file_handle = \*NEW;
    
    &print_sequencing_primers($file_handle, $primer_type, $filename, $sequence_id, $results, length($seq));


    my $primer_list = undef; #
    if ($SEQUENCING_DIRECTION eq "toward3") {
        &print_primer_design_view ($file_handle, $included_region, $excluded_region, $target, $seq, $primer_list, $results, undef);
    } else {
        &print_primer_design_view ($file_handle, $included_region, $excluded_region, $target, $seq, $primer_list, undef, $results);
    }
    
    print NEW "</td></tr></tbody></table>\n";
    print NEW "$wrapup\n";
    close NEW;
    
}

# print sequencing primers
sub print_sequencing_primers {
    my ($file_handle, $primer_type, $filename, $sequence_id, $results, $seq_length) = @_;

    my $count = 0;
                           
    my $has_primers = 0;
    for (my $i = 0; $i < scalar(@$results); $i++) {
        if (@$results[$i]) {
            $has_primers = 1;
            last;
        }
    }
                            
    if ($has_primers) {
        print HTML_REPORT_FILE "<p>Sequencing primers:</p>\n";
        print $file_handle "<p><b>Sequencing primers:</b></p>\n";
        $snp_sucess_seqs++;
    } else {
        $snp_failed_seqs++;
        print HTML_REPORT_FILE "<p><font color=red>No sequencing primer found.</font></p>\n";
        print $file_handle "<p><font color=red>No sequencing primer found.</font></p>\n";
        return;
    }
 
    my $primer_count = 0;
                         
    my $snp_header = "<table><thead><tr><th></th><th>Orientation</th><th>Start</th><th>Len</th><th>Tm</th><th>GC%</th><th>Any compl</th><th>3' compl</th><th>Score</th><th>Primer Seq</th></tr></thead>\n";
    print HTML_REPORT_FILE $snp_header;
    print "<tbody>\n";
    
    print $file_handle $snp_header;
    print $file_handle "<tbody>\n";
    
    
    if ($results){
       for (my $i = 0; $i < scalar(@$results); $i++) {
           if (@$results[$i] ) {
               my @primer_results = split(/\t/, @$results[$i]);
               $primer_count++;
               $snp_total_primers++;
               my $primer_type_str = "Sequencing";
               my $direction = "FORWARD";
               $direction = "REVERSE" if ($SEQUENCING_DIRECTION eq "toward5");
               $snp_total_primers++;
               $table_row_count++;
               &print_single_sequencing_primer($file_handle, \@primer_results, $direction, $sequence_id, $primer_type_str, $primer_count,$table_row_count, $primer_type, $filename, $seq_length);
           }
       }
    }
    if ($has_primers) {
        print HTML_REPORT_FILE "</tbody></table>\n";
        print $file_handle "</tbody></table>\n";
    }
    
}

################  functions for SBE primer design  ###############################

sub find_SBE_primers {
    my ($included_region, $excluded_region, $sequence, $target, $include_snp) = @_;
    my ($snp, $target_pos) = get_snp($sequence, $target);
    
 
    # Step 1: pick all primer candidate from forward and reverse orientation:
    # need input data:
    # $sequence, $target_pos, $min_primer_len,
    # $max_primer_len, $include_snp
 
    # reference array of candidate primers
    my ($candidates_f, $start_pos_f) = &find_candidate_primers($included_region, $excluded_region, $sequence, $target_pos,
        $SNP_PRIMER_MIN_SIZE, $SNP_PRIMER_MAX_SIZE, $include_snp, 'FORWARD');

    my $target_pos_rev = length($sequence) - $target_pos +1;  # if $first_base_index = 1
    if ($first_base_index == 0) {
        $target_pos_rev = length($sequence) - $target_pos;
    }
    my ($candidates_r, $start_pos_r) = &find_candidate_primers($included_region, $excluded_region, &reverse_complement($sequence),
    	$target_pos_rev, $SNP_PRIMER_MIN_SIZE, $SNP_PRIMER_MAX_SIZE, $include_snp, 'REVERSE');

    # Step 2: calcuate Tm of all candidate and scores
    my @tm_f = ();
    my @tm_r = ();
    my @oligo_dg_f = ();
    my @end_oligo_dg_f = ();
    my @oligo_dg_r = ();
    my @end_oligo_dg_r = ();
    
    my @gc_f = ();
    my @gc_r = ();
    my @repeats_f = ();
    my @repeats_r = ();
    my @ns_f = ();
    my @ns_r = ();
    my @scores_f = ();
    my @scores_r = ();
    my @self_comp1_f = ();
    my @self_comp2_f = ();
    my @self_comp1_r = ();
    my @self_comp2_r = (); # 3' complementarity
    
    my $oligo_tm = new OligoTM();
 
#    print "forward primers:<br>\n";
    for (my $i = 0; $i < scalar(@$candidates_f); $i++) { 
    	$tm_f[$i] = $oligo_tm->oligo_tm(@$candidates_f[$i], $SNP_PRIMER_DNA_CONC, $SNP_PRIMER_SALT_CONC);
    #	$oligo_dg_f[$i] = &oligo_dg(@$candidates_f[$i]);
    #	$end_oligo_dg_f[$i] = &end_oligo_dg(@$candidates_f[$i], 5);
    	$gc_f[$i] = &oligo_gc(@$candidates_f[$i]);
        $repeats_f[$i] = &find_repeat_bases(@$candidates_f[$i]);
        $ns_f[$i] = &find_ambiguity_ns(@$candidates_f[$i]);
        ($self_comp1_f[$i], $self_comp2_f[$i]) = &self_complementarity(@$candidates_f[$i]);
        $scores_f[$i] = calculate_primer_score(length(@$candidates_f[$i]),
            $tm_f[$i],$gc_f[$i], $repeats_f[$i], $ns_f[$i], $self_comp1_f[$i], $self_comp2_f[$i]);
    }

#    print "reverse primers:<br>\n";
    for (my $i = 0; $i < scalar(@$candidates_r); $i++) { 
    	$tm_r[$i] = $oligo_tm->oligo_tm(@$candidates_r[$i], $SNP_PRIMER_DNA_CONC, $SNP_PRIMER_SALT_CONC);
    #	$oligo_dg_r[$i] = oligo_dg(@$candidates_r[$i]);
    #	$end_oligo_dg_r[$i] = end_oligo_dg(@$candidates_r[$i], 5);
    	$gc_r[$i] = &oligo_gc(@$candidates_r[$i]);
        $repeats_r[$i] = &find_repeat_bases(@$candidates_r[$i]);
        $ns_r[$i] = &find_ambiguity_ns(@$candidates_r[$i]);
        ($self_comp1_r[$i], $self_comp2_r[$i]) = &self_complementarity(@$candidates_r[$i]);
        $scores_r[$i] = calculate_primer_score(length(@$candidates_r[$i]),
            $tm_r[$i],$gc_r[$i], $repeats_r[$i], $ns_r[$i],$self_comp1_r[$i], $self_comp2_r[$i]);
    }

    # pick one best primer from each direction
    my $results_f = '';
    my $results_r = '';
    my $max_score_f = 0;
    my $selected_index = -1;
    for (my $i = 0; $i < scalar(@$candidates_f); $i++) { 
        if ($scores_f[$i] > $max_score_f) {
            $max_score_f = $scores_f[$i];
            $selected_index = $i;
        }
    }
    if ($selected_index != -1 && $max_score_f > 0) {
        $results_f =
            "Left  $snp $target_pos ".
            "@$start_pos_f[$selected_index] ".
            length(@$candidates_f[$selected_index])." ".
            "$tm_f[$selected_index] $gc_f[$selected_index] ".
            "$self_comp1_f[$selected_index] $self_comp2_f[$selected_index] ".
            "$scores_f[$selected_index] @$candidates_f[$selected_index]";
    }
    
    my $max_score_r = 0;
    $selected_index = -1;
    for (my $i = 0; $i < scalar(@$candidates_r); $i++) { 
        if ($scores_r[$i] > $max_score_r) {
            $max_score_r = $scores_r[$i];
            $selected_index = $i;
        }
    }   
    if ($selected_index != -1 && $max_score_r > 0) {
        $results_r =
            "Right $snp $target_pos ".
            "@$start_pos_r[$selected_index] ".
            length(@$candidates_r[$selected_index])." ".
            "$tm_r[$selected_index] $gc_r[$selected_index] ".
            "$self_comp1_r[$selected_index] $self_comp2_r[$selected_index] ".
            "$scores_r[$selected_index] @$candidates_r[$selected_index]";
            
    }

    # return results: direction, primer sequence, primer start position, length, Tm, score
    return ($results_f, $results_r);
}


# reverse_complement
# A subroutine to compute the reverse complement of DNA sequence
sub reverse_complement {
    my $seq = shift;
    # Reverse the sequence
    my $rev_com = reverse($seq);

    # Complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $rev_com =~ tr/ACGTacgt/TGCAtgca/;

    return $rev_com;
}

sub complement {
    my $seq = shift;
    
    # Complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}



# find_candidate_primers
# Find the candidate primers meeting the requirements
sub find_candidate_primers {
    my ($included_region, $excluded_region, $seq, $target_pos, $min_primer_len,
    	$max_primer_len, $include_snp, $orientation) = @_;

    if ($first_base_index != 0) {
       $target_pos -= $first_base_index;
    }

    if (!$include_snp) { # do not include SNP
       $target_pos--;
    }

    # for included and excluded regions
    my ($included_start, $included_end) = ();
    my @excluded_starts = ();
    my @excluded_ends = ();
    
    if ($included_region) {
       ($included_start, $included_end) = $included_region =~ /(\d+),(\d+)/;
       $included_end += $included_start;  
    }
    
    if ($excluded_region) {
        my @tmp = split(/\s+/, $excluded_region);
        for (my $i = 0; $i < @tmp; $i++) {
            ($excluded_starts[$i], $excluded_ends[$i]) = $tmp[$i] =~ /(\d+),(\d+)/;
            $excluded_ends[$i] += $excluded_starts[$i];
        }
    }

    my @candidates = ();
    my @start_poses = ();
    my $count = 0;
    my $start = ($target_pos - $max_primer_len + 1);
    if ($start < 0) {
    	$start = 0;
    }
    my $end = ($target_pos - $min_primer_len + 1);

    for (my $i = $start; $i <= $end && $end >= 0; $i++) {
        # check if the candidate primer is in the included region or excluded region;
        my $s = $i;
        my $e = $target_pos; # to base index = 0
        if ($orientation eq 'REVERSE') {
            my $tmp = $e;
            $e = length($seq)-$s-1;
            $s = length($seq)-$tmp-1;
        }
        
        if ($included_region) {
            next if ($s < $included_start-1 || $e >$included_end-2);
        }
            
        if ($excluded_region) {
            my $excluded = 0;
            for (my $k = 0; $k < @excluded_starts; $k++) {
                
#                print $s, "   ", $e, "   ", $excluded_starts[$k]-1, "   ", $excluded_ends[$k]-2, "   ", $orientation, " <br>";
                # -1 or -2 to convert the index to the base index 0
                if ($s <= $excluded_starts[$k]-1 &&  $e >= $excluded_ends[$k]-2  ||
                    $s >= $excluded_starts[$k]-1 && $e <= $excluded_ends[$k]-2 ||
                    $s >= $excluded_starts[$k]-1 && $s <= $excluded_ends[$k]-2 && $e >= $excluded_ends[$k]-2 ||
                    $s <= $excluded_starts[$k]-1 && $e >= $excluded_starts[$k]-1 && $e <= $excluded_ends[$k]-2) { 
                     $excluded = 1;
                    last;
                }
            }
            next if $excluded == 1;
        }
            
        # find the candidate primer        
        my $candidate_sprimer = uc(substr($seq, $i, $target_pos - $i + 1)); # sweep out some bad candidate
        my $ns = find_ambiguity_ns($candidate_sprimer);
        next if ($ns > $SNP_PRIMER_MAX_N);

       $candidates[$count] = $candidate_sprimer;
       if ($orientation eq 'REVERSE') {
            $start_poses[$count] = length($seq) - $i + 1;
       } else {
            $start_poses[$count] = $i;
       }
       $count++;
    }
    return (\@candidates, \@start_poses);
}

# Caluclate GC content
sub oligo_gc {
    my $seq = shift;
    my $gc = 0;
    my $length = length($seq);
    for (my $i=0; $i < $length; $i++) {
        if (substr($seq, $i, 1) eq "G" || substr($seq, $i, 1) eq "C" ||
            substr($seq, $i, 1) eq "g" || substr($seq, $i, 1) eq "c") {
            $gc++;
        }
    }
    return $gc/$length*100;
}

sub max {
    my ($v1, $v2) = @_;
    if ($v1 > $v2) {
        return $v1;
    } else {
        return $v2;
    }
}


# Calculate primer score for a primer candidate
sub calculate_primer_score {
    my ($length, $tm, $gc, $repeats, $n, $self_any, $self_end) = @_;
    my $score = 0;

    # primer size score: 20
    
    $score = 20 - abs($length - $SNP_PRIMER_OPT_SIZE)*20/&max(abs($SNP_PRIMER_MAX_SIZE - $SNP_PRIMER_OPT_SIZE), abs($SNP_PRIMER_MIN_SIZE - $SNP_PRIMER_OPT_SIZE));

    # GC score: 10
    if ($gc >= $SNP_PRIMER_MIN_GC && $gc <= $SNP_PRIMER_MAX_GC) {
        $score += 10;
    } else {
        $score = 0;
        return $score;
    }
    
    # Tm score: 30
    if ($tm >= $SNP_PRIMER_MIN_TM && $tm <= $SNP_PRIMER_MAX_TM) {
        $score += 30 - abs($tm - $SNP_PRIMER_OPT_TM)*30/&max(abs($SNP_PRIMER_MAX_TM - $SNP_PRIMER_OPT_TM), abs($SNP_PRIMER_MIN_TM - $SNP_PRIMER_OPT_TM));
    } else {
        $score = 0;
        return $score;
    }

   
   # repeat score: 10
   if ($repeats /$length <= 0.5) {
        $score += 10 - $repeats * 10/$length;
   } else {
        $score = 0;
        return $score;
   }
 

   # ambiguity N score: 20
   if ( $n <= $SNP_PRIMER_MAX_N) {
        if ($SNP_PRIMER_MAX_N > 0) {
            $score += 20 - $n/$SNP_PRIMER_MAX_N * 20;
        } else {
            $score += 20;
        }
   } else {
        $score = 0;
        return $score;
   }
   

   # self_complementarity check 10
   if ($self_any > $SNP_MAX_SELF_COMPLEMENTARITY || $self_end > $SNP_MAX_3_SELF_COMPLEMENTARITY) {
        $score = 0;
        return $score;
   } else {
        $score += 10 - (($SNP_MAX_SELF_COMPLEMENTARITY - $self_any) * 5/$SNP_MAX_SELF_COMPLEMENTARITY +
                        ($SNP_MAX_3_SELF_COMPLEMENTARITY - $self_end) * 5/$SNP_MAX_3_SELF_COMPLEMENTARITY);  
   }
    
   
   return $score;
}

# detect the repeat patterns and calculate the repeat bases
sub find_repeat_bases {
    my $sequence = shift;
    
    my @repeats = (
        [1,6],
        [2,6],  #dinucl. with >= 6 repeats
        [3,4]   #trinucl. with >= 4 repeats
    );  
  
    my $total_repeats = 0;
    
    #find the simple repeats
    study($sequence);     # may be necessary?
    my $seqlength = length($sequence);
    my $ssr_number = 1;   #track multiple ssrs within a single sequence
    my %locations;        #track location of simple repeats as detected
    my $i;
    my $count = 0;
    for($i=0; $i<scalar(@repeats); $i++){ #test each spec against sequence
        my $motiflength = $repeats[$i]->[0];
        my $minreps = $repeats[$i]->[1] - 1;
        my $regexp = "(([gatc]{$motiflength})\\2{$minreps,})";
        while ($sequence =~ /$regexp/ig){
            my $motif = lc($2);
            my $ssr = $1;
            my $ssrlength = length($ssr);          #overall SSR length
  
            my $end = pos($sequence);              #where SSR ends
            pos($sequence) = $end - $motiflength;  #see docs
            my $start = $end - $ssrlength + 1;     #where SSR starts

            if (&novel($start, \%locations)) {   #count SSR only once
            	$total_repeats += $ssrlength;
            }
        }
    }
    return $total_repeats;
}



# find number of ambiguilty Ns
sub find_ambiguity_ns {
    my $seq = shift;
    my $length = length($seq);
    my $ns = 0;
    for (my $i=0; $i < $length; $i++) {
        if (substr($seq, $i, 1) eq "N" || substr($seq, $i, 1) eq "n") {
            $ns++;
        }
    }
    return $ns;
}

# Find allele-specific primers. For each direction, two primers associated with two
# alleles are picked. Total of four primers will be designed if possible.
# Users should pick the best two primers based on the score.
sub find_allele_primers {
    my ($included_region, $excluded_region, $sequence, $target, $include_snp, $second_mismatch, $second_mismatch_pos, $primer_type) = @_;
    my ($snp, $target_pos) = get_snp($sequence, $target);
  
    # Step 1: pick all primer candidate from forward and reverse orientation:

    # reference array of candidate primers
    my ($candidates_f, $start_pos_f) = &find_candidate_primers($included_region, $excluded_region, $sequence, $target_pos,
        $SNP_PRIMER_MIN_SIZE, $SNP_PRIMER_MAX_SIZE, $include_snp, 'FORWARD');

    my $target_pos_rev = length($sequence) - $target_pos +1;  # if $first_base_index = 1
    if ($first_base_index == 0) {
        $target_pos_rev = length($sequence) - $target_pos;
    }
    my ($candidates_r, $start_pos_r) = &find_candidate_primers($included_region, $excluded_region, &reverse_complement($sequence),
    	$target_pos_rev, $SNP_PRIMER_MIN_SIZE, $SNP_PRIMER_MAX_SIZE, $include_snp, 'REVERSE');


    # convert the SNP degeneration code to regular code
    my @alleles = ();
    my $snp_code = "";
        
    if ($candidates_f && scalar(@$candidates_f) > 0) {
        $snp_code = substr(@$candidates_f[0], length(@$candidates_f[0])-1);
    } elsif ($candidates_r  && scalar(@$candidates_r) > 0) {
        $snp_code = substr(@$candidates_r[0], length(@$candidates_r[0])-1);
    }
    @alleles = &convert_code($snp_code);
    
    # add the second mismatch base
    if ($second_mismatch == 1) {
        $candidates_f = &add_second_mismatch($candidates_f, $snp_code, $second_mismatch_pos);  # two dimentional array
        $candidates_r = &add_second_mismatch($candidates_r, $snp_code, $second_mismatch_pos);  # two dimentional array
    }
 
    my @candidates_f_all = &replace_last_code_with_allele($candidates_f, \@alleles, $FORWARD);  # two dimentional array
    my @candidates_r_all = &replace_last_code_with_allele($candidates_r, \@alleles, $REVERSE);  # two dimentional array

    
    # Step 2: calcuate Tm of all candidate and scores
    my @tm_f_all = ();
    my @tm_r_all = ();
    
    my @gc_f_all = ();
    my @gc_r_all = ();
   
    my @repeats_f_all = ();
    my @repeats_r_all = ();

    my @ns_f_all = ();
    my @ns_r_all = ();
 
    my @scores_f_all = ();
    my @scores_r_all = ();
    my @self_comp1_f_all = ();
    my @self_comp2_f_all = ();
    my @self_comp1_r_all = ();
    my @self_comp2_r_all = (); # 3' complementarity
    
    my $oligo_tm = new OligoTM();

    for (my $j = 0; $j < scalar(@alleles); $j++) {
        my @tm = ();
        my @gc = ();
        my @repeats = ();
        my @ns = ();
        my @comp1 = ();
        my @comp2 = ();
        my @scores = ();
        my $tmp = $candidates_f_all[$j];
        for (my $i = 0; $i < scalar(@$tmp); $i++) {
            my $tm0 = $oligo_tm->oligo_tm(@$tmp[$i], $SNP_PRIMER_DNA_CONC, $SNP_PRIMER_SALT_CONC);
            my $gc0 = &oligo_gc(@$tmp[$i]);
            my $repeat0 = &find_repeat_bases(@$tmp[$i]);
            my $ns0 = &find_ambiguity_ns(@$tmp[$i]);

            my ($comp10, $comp20) = &self_complementarity(@$tmp[$i]);
            my $score = calculate_primer_score(length(@$tmp[$i]),
            $tm0,$gc0, $repeat0, $ns0, $comp10, $comp20);
            $tm[$i] = $tm0;
            $gc[$i] = $gc0;
            $repeats[$i] = $repeat0;
            $ns[$i] = $ns0;
            $comp1[$i] = $comp10;
            $comp2[$i] = $comp20;
            $scores[$i] = $score;
        }
        $tm_f_all[$j] = \@tm;
        $gc_f_all[$j] = \@gc;
        $repeats_f_all[$j] = \@repeats;
        $ns_f_all[$j] = \@ns;
        $self_comp1_f_all[$j] = \@comp1;
        $self_comp2_f_all[$j] = \@comp2;
        $scores_f_all[$j] = \@scores; 
 
        #reverse
        my @tm1 = ();
        my  @gc1 = ();
        my  @repeats1 = ();
        my  @ns1 = ();
        my  @comp11 = ();
        my  @comp21 = ();
        my  @scores1 = ();
        my  $tmp1 = $candidates_r_all[$j];
        for (my $i = 0; $i < scalar(@$tmp1); $i++) { 
            my $tm0 = $oligo_tm->oligo_tm(@$tmp1[$i], $SNP_PRIMER_DNA_CONC, $SNP_PRIMER_SALT_CONC);
            my $gc0 = &oligo_gc(@$tmp1[$i]);
            my $repeat0 = &find_repeat_bases(@$tmp1[$i]);
            my $ns0 = &find_ambiguity_ns(@$tmp1[$i]);

            my ($comp10, $comp20) = &self_complementarity(@$tmp1[$i]);
            my $score = calculate_primer_score(length(@$tmp1[$i]),
                $tm0,$gc0, $repeat0, $ns0, $comp10, $comp20);
            $tm1[$i] = $tm0;
            $gc1[$i] = $gc0;
            $repeats1[$i] = $repeat0;
            $ns1[$i] = $ns0;
            $comp11[$i] = $comp10;
            $comp21[$i] = $comp20;
            $scores1[$i] = $score;
        }
        $tm_r_all[$j] = \@tm1;
        $gc_r_all[$j] = \@gc1;
        $repeats_r_all[$j] = \@repeats1;
        $ns_r_all[$j] = \@ns1;
        $self_comp1_r_all[$j] = \@comp11;
        $self_comp2_r_all[$j] = \@comp21;
        $scores_r_all[$j] = \@scores1; 
    }

    # pick one best primer from each direction
    my @results_f_all = ();
    my @results_r_all = ();
   
    my $alleles_str = $snp;
   
    if ($primer_type ne '9') {
        for (my $j = 0; $j < scalar(@alleles); $j++) {
            my ($max_score, $selected_index) = &find_max_score_primer($scores_f_all[$j]);
            if ($primer_type >= 7) {  # 7, 8, 9
                $snp = substr($alleles_str, $j*2, 1);
            }
            if ($selected_index != -1 && $max_score > 0) {
                $results_f_all[$j] =
                    "Left ".
                    "$snp $target_pos ".
                    "@$start_pos_f[$selected_index] ".
                        length($candidates_f_all[$j]->[$selected_index])." ".
                    "$tm_f_all[$j]->[$selected_index] $gc_f_all[$j]->[$selected_index] ".
                    "$self_comp1_f_all[$j]->[$selected_index] $self_comp2_f_all[$j]->[$selected_index] ".
                    "$scores_f_all[$j]->[$selected_index] ".
                    "$candidates_f_all[$j]->[$selected_index]" ;
            }
        }
    
        for (my $j = 0; $j < scalar(@alleles); $j++) {
            my ($max_score, $selected_index) = &find_max_score_primer($scores_r_all[$j]);
            if ($primer_type >= 7) {  # 7, 8, 9
                $snp = substr($alleles_str, $j*2, 1);
            }
            if ($selected_index != -1 && $max_score > 0) {
                $results_r_all[$j] =
                    "Right ".
                    "$snp $target_pos ".
                    "@$start_pos_r[$selected_index] ".
                    length($candidates_r_all[$j]->[$selected_index])." ".
                    "$tm_r_all[$j]->[$selected_index] $gc_r_all[$j]->[$selected_index] ".
                    "$self_comp1_r_all[$j]->[$selected_index] $self_comp2_r_all[$j]->[$selected_index] ".
                    "$scores_r_all[$j]->[$selected_index] ".
                    $candidates_r_all[$j]->[$selected_index] ;
            }
        }
    } else  {  # tetra-primer ARMS PCR
        
        for (my $j = 0; $j < scalar(@alleles); $j++) {   # only two alleles
            my ($selected_index_f, $selected_index_r) = (-1, -1);
            # find_best_primer_pair
            my $tm_diff = 100;
            my $score_sum = 0;
  
            my $scores_f = $scores_f_all[$j];
            my $scores_r = $scores_r_all[1-$j];
 
            my $tm_f = $tm_f_all[$j];
            my $tm_r = $tm_r_all[1-$j];
 
                 
            if (defined($scores_f) && defined($scores_r)) {
                for (my $k1 = 0; $k1 < scalar(@$scores_f); $k1++) {
                    next if (@$scores_f[$k1] < 60);
                    for (my $k2 = 0; $k2 < scalar(@$scores_r); $k2++) {
                        next if (@$scores_r[$k2] < 60);
                        my $diff = abs(@$tm_f[$k1] - @$tm_r[$k2]);
                        my $sum = @$scores_f[$k1] + @$scores_r[$k2];
                        if ($tm_diff > $diff + 1  && $diff < 5 ||
                            abs($tm_diff - $diff) < 1 && $score_sum < $sum && $diff < 5) {
                            $selected_index_f = $k1;
                            $selected_index_r = $k2;
                            $tm_diff = $diff;
                            $score_sum = $sum;
                        }
                    }
                }
            }
   
            if ($selected_index_f >= 0 && $selected_index_r >= 0) {
                $snp = substr($alleles_str, $j*2, 1);
                $results_f_all[$j] = 
                    "Left ".
                    "$snp $target_pos ".
                    "@$start_pos_f[$selected_index_f] ".
                        length($candidates_f_all[$j]->[$selected_index_f])." ".
                    "@$tm_f[$selected_index_f] $gc_f_all[$j]->[$selected_index_f] ".
                    "$self_comp1_f_all[$j]->[$selected_index_f] $self_comp2_f_all[$j]->[$selected_index_f] ".
                    "@$scores_f[$selected_index_f] ".
                    "$candidates_f_all[$j]->[$selected_index_f]" ;
           
                $snp = substr($alleles_str, (1-$j)*2, 1);
                $results_r_all[1-$j] =
                    "Right ".
                    "$snp $target_pos ".
                    "@$start_pos_r[$selected_index_r] ".
                    length($candidates_r_all[1-$j]->[$selected_index_r])." ".
                    "@$tm_r[$selected_index_r] $gc_r_all[1-$j]->[$selected_index_r] ".
                    "$self_comp1_r_all[1-$j]->[$selected_index_r] $self_comp2_r_all[1-$j]->[$selected_index_r] ".
                    "@$scores_r[$selected_index_r] ".
                    $candidates_r_all[1-$j]->[$selected_index_r] ;
            }
        }
  
        
        for (my $j = 0; $j < scalar(@results_f_all); $j++) {   # only two alleles
            if (!$results_f_all[$j]) {
                splice(@results_f_all, $j, 1);
                $j--;
            }
        }
        for (my $j = 0; $j < scalar(@results_r_all); $j++) {   # only two alleles
            if (!$results_r_all[$j]) {
                splice(@results_r_all, $j, 1); 
                $j--;
            }
        }
       
    }
    
    return (\@results_f_all, \@results_r_all);
}

sub find_max_score_primer {
    my $scores = shift;  # array reference 

    my $max_score = 0;
    my $selected_index = -1;
    for (my $i = 0; $i < scalar(@$scores); $i++) {
        if (@$scores[$i] > $max_score) {
            $max_score = @$scores[$i];
            $selected_index = $i;
        }
    }

    return ($max_score, $selected_index);
}


# convert the degeneation code to alles
sub convert_code{

    my $code = shift;
    my @alleles;

    if($code eq "S"){
        $alleles[0] = "G";
        $alleles[1] = "C";
    }elsif($code eq "W"){ 
        $alleles[0] = "A";
        $alleles[1] = "T";
    }elsif($code eq "R"){ 
        $alleles[0] = "G";
        $alleles[1] = "A";
    }elsif($code eq "Y"){ 
        $alleles[0] = "T";
        $alleles[1] = "C";
    }elsif($code eq "K"){ 
        $alleles[0] = "G";
        $alleles[1] = "T";
    }elsif($code eq "M"){ 
        $alleles[0] = "A";
        $alleles[1] = "C";
    }elsif($code eq "V"){ 
        $alleles[0] = "A";
        $alleles[1] = "C";
        $alleles[2] = "G";
    }elsif($code eq "H"){ 
        $alleles[0] = "A";
        $alleles[1] = "C";
        $alleles[2] = "T";
     }elsif($code eq "D"){ 
        $alleles[0] = "A";
        $alleles[1] = "G";
        $alleles[2] = "T";
     }elsif($code eq "B"){ 
        $alleles[0] = "C";
        $alleles[1] = "G";
        $alleles[2] = "T";
    }
    return @alleles;
}

# return alleles and the position
sub get_snp {
    my ($seq, $target) = @_;
    my @target_pos = split(/,/, $target);
    my $code = substr($seq, $target_pos[0]-$first_base_index,$target_pos[1]); 
    my @alleles = &convert_code($code);
    if (scalar(@alleles) == 2) {
        return ("$alleles[0]/$alleles[1]", $target_pos[0]);
    } else {  # == 3
        return ("$alleles[0]/$alleles[1]/$alleles[2]", $target_pos[0]);
    }        
}

sub replace_last_code_with_allele {
    #$candidates is an reference array, $alleles is also an reference array 
    my ($candidates, $alleles, $direction) = @_;
    my @candidates_all = ();  # tow-dimentional array

    for (my $j = 0; $j < scalar(@$alleles); $j++) {
        my @candidates_snp = ();
        for (my $i = 0; $i < scalar(@$candidates); $i++) {
            if ($direction == $FORWARD) {
                $candidates_snp[$i] = substr(@$candidates[$i], 0, length(@$candidates[$i])-1).@$alleles[$j];
            } else {
                $candidates_snp[$i] = substr(@$candidates[$i], 0, length(@$candidates[$i])-1).reverse_complement(@$alleles[$j]);
            }
        }
        $candidates_all[$j] = \@candidates_snp;
     }
    return @candidates_all;
}

sub add_second_mismatch {
    #$candidates is an reference array  
    my ($candidates, $snp_code, $second_mismatch_pos) = @_;
    my @candidates2 = ();  # one-dimentional array
    
    for (my $i = 0; $i < scalar(@$candidates); $i++) {
        $candidates2[$i] = replace_mismatch_base(@$candidates[$i], $snp_code, $second_mismatch_pos);
    }
    return \@candidates2;
}


sub replace_mismatch_base {
    my ($primer_seq, $code, $mismatch_pos) = @_;
    
    # rules for adding second mismatch:
    # Little, S (1997) ARMS analysis of point mutations. In Taylor, G.R., (ed) Laboratory methods for the Detection of Mutations
    # and Polymorphisms in DNA. CRC Press, Boca Raton, FL, pp 45-51
    # Destabilization Strength        Mismatch
    #  Strong                         G/A, C/T
    #  Medium                         C/C, A/A, G/G, T/T
    #  Weak                           C/A, G/T
    
    # (1) A strong mismatch at the 3' end of an allele-specific primer will likely require a week second mismatch and vice visa
    # (2) A "medium" mismatch at the 3' end of an allele-specific primer will likely require a medium second mismatch
    #                    Mismatch at 3' end 
    #    Code Alleles Forward Reverse Mismatch  Type   Second mismatch
    #    R    G/A     G       T       G/T       week   G/A, T/C
    #                 A       C       A/C       week   A/G, C/T
    #    Y    T/C     T       G       G/T       week   G/A, T/C
    #                 C       A       A/C       week   A/G, C/T
    #    S    G/C     G       G       G/G       medium
    #                 C       C       C/C       medium
    #    W    A/T     A       A       A/A       medium
    #                 T       T       T/T       medium
    #    K    G/T     G       A       G/A       strong G/T, A/C
    #                 T       C       T/C       strong T/G, C/A 
    #    M    A/C     A       G       G/A       strong G/T, A/C
    #                 C       T       T/C       strong T/G, C/A

  
    # the key is the nucleotide at the template DNA ate the mismatch position, i.e., the nucleotide after complementation
    # at the mismatch position in a primer sequence   
    my %RY_hash = (
        'G' => 'A',
        'T' => 'C',
        'A' => 'G',
        'C' => 'T'
        );
  
    my %SW_hash = (
        'G' => 'G',
        'T' => 'T',
        'A' => 'A',
        'C' => 'C'
        );
    my %KM_hash = (
        'G' => 'T',
        'T' => 'G',
        'A' => 'C',
        'C' => 'A'
        );

        
    $code = uc($code);
    $primer_seq = uc($primer_seq);

    my $base = substr($primer_seq, $mismatch_pos, 1);
    if ($code eq 'R' || $code eq 'Y') {
        substr($primer_seq, $mismatch_pos, 1, $RY_hash{complement($base)}); 
    }
    elsif ($code eq 'S' || $code eq 'W') {
        substr($primer_seq, $mismatch_pos, 1, $SW_hash{complement($base)}); 
    }
    elsif ($code eq 'K' || $code eq 'M') {
        substr($primer_seq, $mismatch_pos, 1, $KM_hash{complement($base)}); 
    }
    
    return $primer_seq;
}

###################### Primer complementarity check ###############################

# check primer for max complementarity and max 3' complementarity
# return if this primer passed self complementarity and 3' self-complementarity
sub self_complementarity {
    my $primer = shift;
    
    my $rows = length($primer);
    my $comp = reverse_complement($primer);
        
    #  Check the primer for self-priming
    my ($row_max, $row_max_pos) = &self_complementarity_alignment($comp,$primer,$FORWARD);
    
    my $self_end = @$row_max[$rows-1];
    my $self_any = 0;
    for(my $i = 0; $i < $rows; $i++) {
        if (@$row_max[$i] > $self_any) {
            $self_any = @$row_max[$i];
        }
    }
    
    return ($self_any, $self_end);
}

# do complementarity check: local alighment using Smith-Waterman algorithm
sub self_complementarity_alignment {
    my ($template, $new_seq, $direction) = @_;
    
    my @row_max = ();
    my @row_max_pos = ();
    
    my $cols = length($template);   
    my $rows = length($new_seq);
    my $current;
    my $a;
    my $b;
    my $c;
    my @temp;   # int[]
    my $b1; #char 
    my $b2; #char
        
    my @prev_row = ();
    my @curr_row = ();
    
    if($direction == $FORWARD) {
        $b1 = substr($template, 0, 1);
        for(my $y = 0; $y < $rows; $y++) {
            if (substr($new_seq, $y, 1) eq $b1) {
                $current = $MATCH;
            } else {
                $current = 0;
            }
            $row_max[$y] = $current;
            $row_max_pos[$y] = 0;
        }
        
        $b1 = substr($new_seq, 0, 1);
        for(my $x = 0; $x < $cols; $x++) {
            if (substr($template, $x, 1) eq $b1) {
                $current = $MATCH;
            } else {
                $current = 0;
            }
            $prev_row[$x] = $current;
            if($current > $row_max[0]) {
                $row_max[0] = $current;
                $row_max_pos[0] = $x;
            }
        }
            
        for(my $y = 1; $y < $rows; $y++) {
            $b2 = substr($new_seq, $y, 1);
            if (substr($template, 0, 1) eq $b2) {
                $curr_row[0] = $MATCH;
            } else {
                $curr_row[0] = 0;
            }
            for(my $x = 1; $x < $cols; $x++) {
                $b1 = substr($template, $x, 1);
                if($b1 eq 'N' || $b2 eq 'N') {
                    $a = $prev_row[$x-1] + $N_SCORE;
                } else {
                    $a = $prev_row[$x-1] + ($b1 eq $b2 ? $MATCH : $MISMATCH);
                }
                $b = $prev_row[$x] + $GAP;
                $c = $curr_row[$x-1] + $GAP;
                    
                $current = ($a>$b ? ($a>$c ? $a:$c) : ($b>$c ? $b:$c));
                if($current<0) {
                    $current = 0;
                }
                $curr_row[$x] = $current;
                if($current >= $row_max[$y]) {
                    $row_max[$y] = $current;
                    $row_max_pos[$y] = $x;
                }
            }
            @temp = @prev_row;
            @prev_row = @curr_row;
            @curr_row = @temp;
        }
    } else {  # reverse
        $b1 = substr($template, $cols-1, 1);
        for(my $y = $rows-1; $y >= 0; $y--) {
            $current = (substr($new_seq, $y, 1) eq $b1 ? $MATCH : 0);
            $row_max[$y] = $current;
            $row_max_pos[$y] = $cols-1;
        }
        $b1 = substr($new_seq, $rows-1, 1);
        for(my $x = $cols-1; $x >= 0; $x--) {
            $current = (substr($template, $x, 1) eq $b1 ? $MATCH : 0);
            $prev_row[$x] = $current;
            if($current > $row_max[$rows-1]) {
                $row_max[$rows-1] = $current;
                $row_max_pos[$rows-1] = $x;
            }
        }
            
        for(my $y = $rows-2; $y >= 0; $y--) {
            $b2 = substr($new_seq, $y, 1);
            $curr_row[$cols-1] = (substr($template, $cols-1, 1) eq $b2) ? $MATCH : 0;
            for(my $x = $cols-2; $x >= 0; $x--) {
                $b1 = substr($template, $x, 1);
                if ($b1 == 'N' || $b2 == 'N') {
                    $a = $prev_row[$x+1] + $N_SCORE;
                } else {
                    $a = $prev_row[$x+1] + ($b1 eq $b2 ? $MATCH : $MISMATCH);
                }
                $b = $curr_row[$x+1] + $GAP;
                $c = $prev_row[$x] + $GAP;
                    
                $current = ($a>$b ? ($a>$c ? $a:$c) : ($b>$c ? $b:$c));
                if($current<0) {
                    $current = 0;
                }
                $curr_row[$x] = $current;
                if($current >= $row_max[$y]) {
                    $row_max[$y] = $current;
                    $row_max_pos[$y] = $x;
                }
            }
            @temp = @prev_row;
            @prev_row = @curr_row;
            @curr_row = @temp;
        }
    }
    return (\@row_max, \@row_max_pos);
}


############# The following functions is to design tetra-primer ARMS Primers (Ye et al 2000) ##################
    
# design outer primers for ARMS PCR
sub design_ARMS_outer_primers {
    my ($included, $excluded, $target, $seq, $snp_f_list, $snp_r_list, $input, $cmd, $sequence_id) = @_;
    
    # clean inner primer list array because some element in the array is not defined
    for (my $k = 0; $snp_r_list && $k < @$snp_r_list ; $k++) {
        if (!@$snp_r_list[$k]) {
            splice(@$snp_r_list, $k, 1);
            $k--;
        }
    }
    for (my $k = 0; $snp_f_list && $k < @$snp_f_list; $k++) {
        if (!@$snp_f_list[$k]) {
            splice(@$snp_f_list, $k, 1);
            $k--;
        }
    }

    # step 1: check if there is allele-specific primer found
    if (!$snp_f_list || scalar(@$snp_f_list) == 0) {
        #print error message
        print HTML_REPORT_FILE "<font color=red>No forward inner primers found </font> in $sequence_id!<br>\n";
        return;
    } elsif (!$snp_r_list || scalar(@$snp_r_list) == 0 ) {
        #print error message
        print HTML_REPORT_FILE "<font color=red>No reverse inner primers found</font> in $sequence_id!<br>\n";
        return;
    }
    
    # Step 2: calculate the forward and reverse maximum product size
    my @results = split(/\s+/, @$snp_f_list[0]);
    my $target_pos = $results[2];
    my $len = length($seq);

    my $max_prod_size_f = $target_pos;
    my $max_prod_size_r = $len - $target_pos + 1;
    
    # Step 3
    if ($max_prod_size_f < $SNP_INNER_PRODUCT_MIN_SIZE ||
        $max_prod_size_r < $SNP_INNER_PRODUCT_MIN_SIZE) {
        # print error message: and adjust the minimum product size
        print HTML_REPORT_FILE "<font color=red>Inner product size > sequence length</font> in $sequence_id!<br>\n";
        return;
    }
    
    # Step 4:
    
    if ($max_prod_size_f > $max_prod_size_r &&
        $max_prod_size_f / $max_prod_size_r < $SNP_INNER_PRODUCT_MIN_DIFF) {
        #print message
        print HTML_REPORT_FILE "<font color=red>Two inner product sizes can not meet the minimum relative size ratio $SNP_INNER_PRODUCT_MIN_DIFF </font> in $sequence_id!<br>\n";
        return;
    }
    if ($max_prod_size_f < $max_prod_size_r &&
        $max_prod_size_r / $max_prod_size_f < $SNP_INNER_PRODUCT_MIN_DIFF) {
        #print message
        print HTML_REPORT_FILE "<font color=red>Two inner product sizes can not meet the minimum relative size ratio $SNP_INNER_PRODUCT_MIN_DIFF </font> in $sequence_id!<br>\n";
        return;
    }


    # check if there is paired forward and reverse inner primer 
    if (scalar(@$snp_f_list) +  scalar(@$snp_r_list) < 4) {
        # case 1: @$snp_f_list = 2 and @$snp_r_list = 1
        if (scalar(@$snp_f_list) == 2 && scalar(@$snp_r_list) == 1) {
            
            my @tmp1 = split(/\s+/, @$snp_r_list[0]);
            for (my $k = 0; $k < @$snp_f_list; $k++) {
                my @tmp2 = split(/\s+/, @$snp_f_list[$k]);
                
                if ($tmp1[1] eq $tmp2[1]) {
                    splice(@$snp_f_list, $k, 1);
                    $k--;
                }
            }
            return if (scalar(@$snp_f_list) == 0);
        }
        
        # case 2: @$snp_f_list = 1 and @$snp_r_list = 2
        elsif (scalar(@$snp_f_list) == 1 && scalar(@$snp_r_list) == 2) {
            my @tmp1 = split(/\s+/, @$snp_f_list[0]);
            for (my $k = 0; $k < @$snp_r_list; $k++) {
                my @tmp2 = split(/\s+/, @$snp_r_list[$k]);
                if ($tmp1[1] eq $tmp2[1]) {
                    splice(@$snp_r_list, $k, 1);
                    $k--;
                }
            }
            return if (scalar(@$snp_r_list) == 0);
        }
        
        # case 3: @$snp_f_list = 1 and @$snp_r_list = 1
        elsif (scalar(@$snp_f_list) == 1 && scalar(@$snp_r_list) == 1) {
            my @tmp1 = split(/\s+/, @$snp_f_list[0]);
            my @tmp2 = split(/\s+/, @$snp_r_list[0]);
            if ($tmp1[1] eq $tmp2[1]) {
                return;
            }
        }
    }
    
    # Step 5: set parameters for picking out primers and pick out primers

    my @primer3outputs = ();
    my @inner_primers_f = ();
    my @inner_primers_r = ();
    my $count = 0;
    

    for (my $j = 0; $j < @$input;  $j++) {
        # remove the target
        if (@$input[$j] =~ m/TARGET=/ ||
            @$input[$j] =~ m/SEQUENCE=/ ||
            @$input[$j] =~ m/PRIMER_OPT_TM=/ ||
            @$input[$j] =~ m/PRIMER_MIN_TM=/ ||
            @$input[$j] =~ m/PRIMER_MAX_TM=/ ||
            @$input[$j] =~ m/PRIMER_LEFT_INPUT=/ ||
            @$input[$j] =~ m/PRIMER_RIGHT_INPUT=/ ||
            @$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/ ||
            @$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/ ||
            
            @$input[$j] =~ m/EXCLUDED_REGION=/ ||
            @$input[$j] =~ m/INCLUDED_REGION=/ ||
            
            @$input[$j] =~ m/PRIMER_MIN_SIZE=/ ||
            @$input[$j] =~ m/PRIMER_OPT_SIZE=/ ||
            @$input[$j] =~ m/PRIMER_MAX_SIZE=/ ||
            @$input[$j] =~ m/PRIMER_MIN_GC=/ ||
            @$input[$j] =~ m/PRIMER_MAX_GC=/ ||
            @$input[$j] =~ m/PRIMER_SELF_ANY=/ ||
            @$input[$j] =~ m/PRIMER_SELF_END=/ ||
            @$input[$j] =~ m/PRIMER_SALT_CONC=/ ||
            @$input[$j] =~ m/PRIMER_DNA_CONC=/ ||
            @$input[$j] =~ m/PRIMER_MAX_DIFF_TM=/ ||
            @$input[$j] =~ m/PRIMER_LIBERAL_BASE=/ ||
            @$input[$j] =~ m/PRIMER_NUM_RETURN=/ ||
            @$input[$j] =~ m/PRIMER_NUM_NS_ACCEPTED=/ ) {

            splice(@$input, $j, 1);
            $j--;
        }
    }
    
    # insert common parameter records
    splice(@$input, 0, 0, "PRIMER_MIN_SIZE=$SNP_PRIMER_MIN_SIZE\n");
    splice(@$input, 0, 0, "PRIMER_OPT_SIZE=$SNP_PRIMER_OPT_SIZE\n");
    splice(@$input, 0, 0, "PRIMER_MAX_SIZE=$SNP_PRIMER_MAX_SIZE\n");
    splice(@$input, 0, 0, "PRIMER_MAX_DIFF_TM=$PRIMER_MAX_DIFF_TM\n");
    splice(@$input, 0, 0, "PRIMER_MIN_GC=$SNP_PRIMER_MIN_GC\n");
    splice(@$input, 0, 0, "PRIMER_MAX_GC=$SNP_PRIMER_MAX_GC\n");   
    splice(@$input, 0, 0, "PRIMER_SELF_ANY=$SNP_MAX_SELF_COMPLEMENTARITY\n");
    splice(@$input, 0, 0, "PRIMER_SELF_END=$SNP_MAX_3_SELF_COMPLEMENTARITY\n");
    splice(@$input, 0, 0, "PRIMER_SALT_CONC=$SNP_PRIMER_SALT_CONC\n");
    splice(@$input, 0, 0, "PRIMER_DNA_CONC=$SNP_PRIMER_DNA_CONC\n");
    splice(@$input, 0, 0, "PRIMER_NUM_NS_ACCEPTED=$SNP_PRIMER_MAX_N\n");
    splice(@$input, 0, 0, "PRIMER_LIBERAL_BASE=1\n");
    splice(@$input, 0, 0, "PRIMER_NUM_RETURN=1\n");
    splice(@$input, 0, 0, "EXCLUDED_REGION=$excluded\n") if ($excluded);
    splice(@$input, 0, 0, "INCLUDED_REGION=$included\n") if ($included);
   
    for (my $i = 0; $i < scalar(@$snp_f_list); $i++) {
        my @results_f = split(/\s+/, @$snp_f_list[$i]);
        my @results_r = undef;
        if (scalar(@$snp_r_list) == 1) {
            @results_r = split(/\s+/, @$snp_r_list[0]);
        } else {
            @results_r = split(/\s+/, @$snp_r_list[1-$i]);
        }
        
        my $seq_f = $results_f[10];
        my $seq_r = $results_r[10];
        my $tm_f = $results_f[5];
        my $tm_r = $results_r[5];
        my $tm = int(($tm_f +  $tm_r)/2);
     
        my $size_range = "";
        if ($max_prod_size_f > $max_prod_size_r) {

            # choose reverse outer primer first
            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                
                if (@$input[$j] =~ m/PRIMER_LEFT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }

                if (@$input[$j] =~ m/PRIMER_RIGHT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/SEQUENCE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
            }

            my $max_product_size = int($max_prod_size_f/($SNP_INNER_PRODUCT_MIN_DIFF + ($SNP_INNER_PRODUCT_MAX_DIFF - $SNP_INNER_PRODUCT_MIN_DIFF)/2));
            $size_range = int($SNP_INNER_PRODUCT_MIN_SIZE) . "-$max_product_size";
            
            my $last_char = substr($seq_f, length($seq_f)-1,1);
            #replace allele in sequence
            my $seq1 = $seq;
            substr($seq1, $target_pos-1, 1, $last_char);
            #replace the second mismatch nucleotide
            my $mismatch_char = substr($seq_f, length($seq_f)-3, 1);
            substr($seq1, $target_pos-3, 1, $mismatch_char);
 
            splice(@$input, 0, 0, "PRIMER_OPT_TM=$tm\n");
            splice(@$input, 0, 0, "PRIMER_MIN_TM=$SNP_PRIMER_MIN_TM\n");
            splice(@$input, 0, 0, "PRIMER_MAX_TM=$SNP_PRIMER_MAX_TM\n");
            splice(@$input, 0, 0, "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n");
            splice(@$input, 0, 0, "PRIMER_PRODUCT_OPT_SIZE=$SNP_INNER_PRODUCT_OPT_SIZE\n");
            splice(@$input, 0, 0, "PRIMER_LEFT_INPUT=$seq_f\n");
            splice(@$input, 0, 0, "SEQUENCE=$seq1\n");

	    my @results = &run_primer3($cmd, $input);
            my $primer3output = new Primer3Output();
            $primer3output->set_results(\@results);
            
            my $primer_list = $primer3output->get_primer_list();
            if (!$primer_list || scalar(@$primer_list) == 0) {
                #print error message: no reverse outer primer found.
                print HTML_REPORT_FILE "<font color=red>No reverse outter primer found </font> in $sequence_id!<br>\n";
                next;
            }
            
            my $primer_pair_1 = @$primer_list[0];
            my $inner_primer_f_obj = $primer_pair_1->left_primer();   # forward inner primer, allele-specific primer
            my $outer_primer_r_obj = $primer_pair_1->right_primer(); # reverse outer primer

            $max_product_size = int($primer_pair_1->product_size() * $SNP_INNER_PRODUCT_MAX_DIFF);
            $size_range = int($primer_pair_1->product_size() * $SNP_INNER_PRODUCT_MIN_DIFF) . "-" . $max_product_size;
            my $opt_product_size = int($primer_pair_1->product_size() * ($SNP_INNER_PRODUCT_MAX_DIFF + $SNP_INNER_PRODUCT_MIN_DIFF)/2);
            
            if (!$outer_primer_r_obj) {
                # print message            
                next;
            }

            # then choose forward outer primer
            $last_char = &complement(substr($seq_r, length($seq_r)-1,1));
            #replace allele in sequence
            my $seq2 = $seq;
            substr($seq2, $target_pos-1, 1, $last_char);
         
            #replace the second mismatch nucleotide
            $mismatch_char = &complement(substr($seq_r, length($seq_r)-3, 1));
            substr($seq2, $target_pos+1, 1, $mismatch_char);

            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_OPT_SIZE=$opt_product_size\n";
                }
                
                if (@$input[$j] =~ m/PRIMER_LEFT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }

                if (@$input[$j] =~ m/SEQUENCE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
            }
 
            splice(@$input, 0, 0, "SEQUENCE=$seq2\n");
            splice(@$input, 0, 0, "PRIMER_RIGHT_INPUT=$seq_r\n");


	    my @results_2 = &run_primer3($cmd, $input);
            my $primer3output2 = new Primer3Output();
            $primer3output2->set_results(\@results_2);
            
            $primer_list = $primer3output2->get_primer_list();
            if (!$primer_list || scalar(@$primer_list) == 0) {
                #print error message: no forward outer primer found.
                next;
            }
                
            my $primer_pair_2 = @$primer_list[0];
            my $inner_primer_r_obj = $primer_pair_2->right_primer();   # reverse inner primer, allele-specific primer
            my $outer_primer_f_obj = $primer_pair_2->left_primer(); # forward outer primer

            if (!$outer_primer_f_obj) {
                # print message            
                next;
            }
            
            
            # last calculate the complementarity of outer primers
                        
            $size_range = int($SNP_INNER_PRODUCT_MIN_SIZE) . "-" . int($primer_pair_1->product_size() + $primer_pair_2->product_size() + 100);
            
            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                
                if (@$input[$j] =~ m/PRIMER_LEFT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/PRIMER_RIGHT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
            }
            
            splice(@$input, 0, 0, "PRIMER_RIGHT_INPUT=" . $outer_primer_r_obj->sequence . "\n");
            splice(@$input, 0, 0, "PRIMER_LEFT_INPUT=" . $outer_primer_f_obj->sequence . "\n");

            
	    my @results_3 = &run_primer3($cmd, $input);
            my $primer3output3 = new Primer3Output();
            $primer3output3->set_results(\@results_3);
            
           
            my $primer_pair_3 = undef;
            $primer_list = $primer3output3->get_primer_list();
            if ($primer_list && scalar(@$primer_list) > 0) {
                $primer_pair_3 = @$primer_list[0];
            } else {
                #print error message: no forward outer primer found.
            }
            
            
            # primer_pair_2: reverse inner primer and forward outer primer
            # primer_pair_1: forward inner primer and reverse outer primer
            if ($primer3output3) {
                $primer3outputs[$count] = $primer3output3;
            
                my $snp_f = $results_f[1];
                my $target_pos_f = $results_f[2];
                my $score_f = $results_f[9];
    
                my $snp_r = $results_r[1];
                my $target_pos_r = $results_r[2];
                my $score_r = $results_r[9];
            
                $inner_primers_f[$count] = 
                    "Left ".
                    "$snp_f ".
                    "$target_pos_f ".
                    ($inner_primer_f_obj->start() - 1)." ".
                    $inner_primer_f_obj->length()." ".
                    $inner_primer_f_obj->tm()." ".
                    $inner_primer_f_obj->gc()." ".
                    $inner_primer_f_obj->any_complementarity()." ".
                    $inner_primer_f_obj->end_complementarity()." ".
                    $score_f." ".
                    $inner_primer_f_obj->sequence(). " ".
                    $primer_pair_1->product_size(). " ".
                    $primer_pair_1->seq_size(). " ".
                    $primer_pair_1->included_size. " ".
                    $primer_pair_1->pair_any_complementarity(). " ".
                    $primer_pair_1->pair_end_complementarity();
                
                $inner_primers_r[$count] = 
                    "Right ".
                    "$snp_r ".
                    "$target_pos_r ".
                    ($inner_primer_r_obj->start() + 1)." ".
                    $inner_primer_r_obj->length()." ".
                    $inner_primer_r_obj->tm()." ".
                    $inner_primer_r_obj->gc()." ".
                    $inner_primer_r_obj->any_complementarity()." ".
                    $inner_primer_r_obj->end_complementarity()." ".
                    $score_r." ".
                    $inner_primer_r_obj->sequence(). " ".
                    $primer_pair_2->product_size(). " ".
                    $primer_pair_2->seq_size(). " ".
                    $primer_pair_2->included_size. " ".
                    $primer_pair_2->pair_any_complementarity(). " ".
                    $primer_pair_2->pair_end_complementarity();
                
                $count++;
            }            
        }
        else {
            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                
                if (@$input[$j] =~ m/PRIMER_LEFT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }

                if (@$input[$j] =~ m/PRIMER_RIGHT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/SEQUENCE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
            }

            
            #***************************************************          
            # choose forward outer primer first
            my $max_product_size = int($max_prod_size_r/($SNP_INNER_PRODUCT_MIN_DIFF + ($SNP_INNER_PRODUCT_MIN_DIFF + $SNP_INNER_PRODUCT_MAX_DIFF)/2));
            if ($max_product_size < $SNP_INNER_PRODUCT_MIN_SIZE) {
                $max_product_size =  $SNP_INNER_PRODUCT_MAX_SIZE;
            }
            $size_range = int($SNP_INNER_PRODUCT_MIN_SIZE) . "-$max_product_size";

            my $last_char = &complement(substr($seq_r, length($seq_r)-1,1));
            #replace allele in sequence
            my $seq1 = $seq;
            substr($seq1, $target_pos-1, 1, $last_char);
         
            #replace the second mismatch nucleotide
            my $mismatch_char = &complement(substr($seq_r, length($seq_r)-3, 1));
            substr($seq1, $target_pos+1, 1, $mismatch_char);
            
            splice(@$input, 0, 0, "PRIMER_OPT_TM=$tm\n");
            splice(@$input, 0, 0, "PRIMER_MIN_TM=$SNP_PRIMER_MIN_TM\n");
            splice(@$input, 0, 0, "PRIMER_MAX_TM=$SNP_PRIMER_MAX_TM\n");
            splice(@$input, 0, 0, "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n");
            splice(@$input, 0, 0, "PRIMER_PRODUCT_OPT_SIZE=$SNP_INNER_PRODUCT_OPT_SIZE\n");
            splice(@$input, 0, 0, "PRIMER_RIGHT_INPUT=$seq_r\n");
            splice(@$input, 0, 0, "SEQUENCE=$seq1\n");


	    my @results = &run_primer3($cmd, $input);
            my $primer3output = new Primer3Output();
            $primer3output->set_results(\@results);
       
            my $primer_list = $primer3output->get_primer_list();
            if (!$primer_list || scalar(@$primer_list) == 0) {
                #print error message: no forward outer primer found.
                print HTML_REPORT_FILE "<font color=red>No forward outer primers found </font> in $sequence_id!<br>\n";
                next;
            }
            
            my $primer_pair_1 = @$primer_list[0];
            my $inner_primer_r_obj = $primer_pair_1->right_primer();   # reverse inner primer, allele-specific primer
            my $outer_primer_f_obj = $primer_pair_1->left_primer();    # forward outer primer

            if (!$outer_primer_f_obj) {
                # print message            
                print HTML_REPORT_FILE "<font color=red>No forward outer primers found </font> in $sequence_id!<br>\n";
                next;
            }
 
            #***************************************************          
            # Then choose reverse outer primer
 
            $max_product_size = int($primer_pair_1->product_size() * $SNP_INNER_PRODUCT_MAX_DIFF);
            
            $size_range = int($primer_pair_1->product_size() * $SNP_INNER_PRODUCT_MIN_DIFF) . "-" . $max_product_size;
            my $opt_product_size = int($primer_pair_1->product_size() * ($SNP_INNER_PRODUCT_MAX_DIFF + $SNP_INNER_PRODUCT_MIN_DIFF)/2);
            
            $last_char = substr($seq_f, length($seq_f)-1,1);
            #replace allele in sequence
            my $seq2 = $seq;
            substr($seq2, $target_pos-1, 1, $last_char);
            #replace the second mismatch nucleotide
            $mismatch_char = substr($seq_f, length($seq_f)-3, 1);
            substr($seq2, $target_pos-3, 1, $mismatch_char);
 
  
            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_OPT_SIZE=$opt_product_size\n";
                }
                if (@$input[$j] =~ m/PRIMER_RIGHT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                
                if (@$input[$j] =~ m/SEQUENCE=/) {
                    @$input[$j] = "SEQUENCE=$seq2\n";
                }
            }
            
            splice(@$input, 0, 0, "PRIMER_LEFT_INPUT=$seq_f\n");

   
	    my @results_2 = &run_primer3($cmd, $input);
            my $primer3output2 = new Primer3Output();
            $primer3output2->set_results(\@results_2);
            
            $primer_list = $primer3output2->get_primer_list();
            if (!$primer_list || scalar(@$primer_list) == 0) {
                #print error message: no reverse outer primer found.
                print HTML_REPORT_FILE "<font color=red>No reverse outer primers found </font> in $sequence_id!<br>\n";
                next;
            }

            #***************************************************          
            # Last, calculate the complementarity of outer primers

                
            my $primer_pair_2 = @$primer_list[0];
            my $inner_primer_f_obj = $primer_pair_2->left_primer();   # forward inner primer, allele-specific primer
            my $outer_primer_r_obj = $primer_pair_2->right_primer();  # reverse outer primer

            if (!$outer_primer_r_obj) {
                # print message            
                print HTML_REPORT_FILE "<font color=red>No reverse outer primers found </font> in $sequence_id!<br>\n";
                next;
            }
                          
            $size_range = int($SNP_INNER_PRODUCT_MIN_SIZE). "-" . int($primer_pair_1->product_size() + $primer_pair_2->product_size() + 100);
            
            for (my $j = 0; $j < @$input; $j++) {
                if (@$input[$j] =~ m/PRIMER_PRODUCT_SIZE_RANGE=/) {
                    @$input[$j] = "PRIMER_PRODUCT_SIZE_RANGE=$size_range\n";
                }
                if (@$input[$j] =~ m/PRIMER_PRODUCT_OPT_SIZE=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                
                if (@$input[$j] =~ m/PRIMER_LEFT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
                if (@$input[$j] =~ m/PRIMER_RIGHT_INPUT=/) {
                    splice(@$input, $j, 1);
                    $j--;
                }
            }
            splice(@$input, 0, 0, "PRIMER_RIGHT_INPUT=" . $outer_primer_r_obj->sequence . "\n");
            splice(@$input, 0, 0, "PRIMER_LEFT_INPUT=" . $outer_primer_f_obj->sequence . "\n");
            
 	    my @results_3 = &run_primer3($cmd, $input);
            my $primer3output3 = new Primer3Output();
            $primer3output3->set_results(\@results_3);
            
            my $primer_pair_3 = undef;
            $primer_list = $primer3output->get_primer_list();
            if ($primer_list && scalar(@$primer_list) > 0) {
                $primer_pair_3 = @$primer_list[0];
            } else {
                #print error message: no forward outer primer found.
                next;
            }
            
            
            # primer_pair_1: reverse inner primer and forward outer primer
            # primer_pair_2: forward inner primer and reverse outer primer
            if ($primer3output3) {
                $primer3outputs[$count] = $primer3output3;
            
                my $snp_f = $results_f[1];
                my $target_pos_f = $results_f[2];
                my $score_f = $results_f[9];
    
                my $snp_r = $results_r[1];
                my $target_pos_r = $results_r[2];
                my $score_r = $results_r[9];
            
                $inner_primers_f[$count] = 
                    "Left ".
                    "$snp_f ".
                    "$target_pos_f ".
                    ($inner_primer_f_obj->start()-1)." ".
                    $inner_primer_f_obj->length()." ".
                    $inner_primer_f_obj->tm()." ".
                    $inner_primer_f_obj->gc()." ".
                    $inner_primer_f_obj->any_complementarity()." ".
                    $inner_primer_f_obj->end_complementarity()." ".
                    $score_f." ".
                    $inner_primer_f_obj->sequence(). " ".
                    $primer_pair_2->product_size(). " ".
                    $primer_pair_2->seq_size(). " ".
                    $primer_pair_2->included_size. " ".
                    $primer_pair_2->pair_any_complementarity(). " ".
                    $primer_pair_2->pair_end_complementarity();
                
                $inner_primers_r[$count] = 
                    "Right ".
                    "$snp_r ".
                    "$target_pos_r ".
                    ($inner_primer_r_obj->start()+1)." ".
                    $inner_primer_r_obj->length()." ".
                    $inner_primer_r_obj->tm()." ".
                    $inner_primer_r_obj->gc()." ".
                    $inner_primer_r_obj->any_complementarity()." ".
                    $inner_primer_r_obj->end_complementarity()." ".
                    $score_r." ".
                    $inner_primer_r_obj->sequence(). " ".
                    $primer_pair_1->product_size(). " ".
                    $primer_pair_1->seq_size(). " ".
                    $primer_pair_1->included_size. " ".
                    $primer_pair_1->pair_any_complementarity(). " ".
                    $primer_pair_1->pair_end_complementarity();
                    $count++;
            }
        }
    }
    
    return (\@primer3outputs, \@inner_primers_f, \@inner_primers_r);
 
}


########################### this function is to create primer view page for each sequence #################################

# print a view to show primer design results on the input sequence for all types of primers thus this is a
# complicated function.  
sub print_primer_design_view {
    my ($file_handle, $included, $excluded, $target, $seq, $primer_list, $snp_f_list, $snp_r_list) = @_;
    
    my $line_width = 60;
    
    #split the sequence into arrays in line_width
    my @seq_rows = &split_sequence($seq, $line_width);
    
    # hilight the sequence using colors for included, excluded and target regions
    @seq_rows = highlight_seq_rows(\@seq_rows, $included, $excluded, $target, $line_width);
    
    # print primers
    # find the maximum rows for primers
    my $max_rows = 0;
    
    $max_rows += scalar(@$primer_list) if ($primer_list);
    $max_rows += scalar(@$snp_f_list) if ($snp_f_list);
    $max_rows += scalar(@$snp_r_list) if ($snp_r_list);
    if (($primer_type eq '7' || $primer_type eq '8') && $SECOND_MISMATCH == 1) {
        if (!$snp_r_list   || @$snp_r_list == 0) {
            $max_rows++;
        }
    }
    
    
    # create an array with two dimentions and init with empty string
    my @snp_rows = ();
    for (my $i = 0; $i < scalar(@seq_rows); $i++) {
        my @rows = ();
        for (my $j = 0; $j < $max_rows; $j++) {
            $rows[$j] = "-" x $line_width;
        }
        $snp_rows[$i] = \@rows;
    }
    
    #add generic /franking primers
    my $count_f = 0;
    my $count_r = 0;
    my $count_o = 0;
    foreach my $primer_pair (@$primer_list) {
        #print left primer
	my $left_primer_obj = $primer_pair->left_primer();
	if ($left_primer_obj) {
            my ($start, $len) = ($left_primer_obj->start(), $left_primer_obj->length());
            my $end_row_index = int(($start+$len-1) / $line_width);
            my $r1 = ($start+$len-1) % $line_width;
            $end_row_index-- if ($r1 == 0);
            my $end_loc = $r1-1;
            $end_loc = $line_width if ($r1 == 0);
        
            my $start_row_index = int($start / $line_width);
            my $r2 = $start % $line_width;
            $start_row_index-- if ($r2 == 0);
            my $start_loc = $r2-1;
            $start_loc = $line_width - 1 if ($r2 == 0);

   
            if ($start_row_index == $end_row_index) {
                my $rows = $snp_rows[$start_row_index];
                substr(@$rows[$count_f], $start_loc, $len, "<font color=#000fff>" . (">" x $len) . "</font>");
            } else {
                my $rows = $snp_rows[$start_row_index];
                my $rows1 = $snp_rows[$start_row_index+1];
                my $len1 = $line_width-$start_loc;
                my $len2 = $len-$len1;
                
                substr(@$rows[$count_f], $start_loc, $len1, "<font color=#000fff>" . (">" x $len1) . "</font>");
                substr(@$rows1[$count_f], 0, $len2, "<font color=#000fff>" . (">" x $len2) . "</font>");
            }
            $count_f++;
        } 
	
        #print right primer
	my $right_primer_obj = $primer_pair->right_primer();
	if ($right_primer_obj) {
            my ($start, $len) = ($right_primer_obj->start(), $right_primer_obj->length());
            $start -= $len;
            
            my $end_row_index = int(($start+$len-1) / $line_width);
            my $r1 = ($start+$len-1) % $line_width;
            $end_row_index-- if ($r1 == 0);
            my $end_loc = $r1;
            $end_loc = $line_width if ($r1 == 0);
        
            my $start_row_index = int($start / $line_width);
            my $r2 = $start % $line_width;
            $start_row_index-- if ($r2 == 0);
            my $start_loc = $r2;
            $start_loc = $line_width - 1 if ($r2 == 0);
            
            if ($start_row_index == $end_row_index) {
                my $rows = $snp_rows[$start_row_index];
                substr(@$rows[$count_r], $start_loc, $len, "<font color=#993300>" . ("<" x $len) . "</font>");
            } else {
                my $rows = $snp_rows[$start_row_index];
                my $rows1 = $snp_rows[$start_row_index +1];
                my $len1 = $line_width-$start_loc;
                my $len2 = $len-$len1;
                substr(@$rows[$count_r], $start_loc, $len1, "<font color=#993300>" . ("<" x $len1). "</font>");
                substr(@$rows1[$count_r], 0, $len2, "<font color=#993300>" . ("<" x $len2). "</font>");
            }
            $count_r++;
        } 
	
        my $oligo_primer_obj = $primer_pair->oligo_primer();
	if ($oligo_primer_obj) {
            my ($start, $len) = ($oligo_primer_obj->start(), $oligo_primer_obj->length());
            my $end_row_index = int(($start+$len-1) / $line_width);
            my $r1 = ($start+$len-1) % $line_width;
            $end_row_index-- if ($r1 == 0);
            my $end_loc = $r1-1;
            $end_loc = $line_width if ($r1 == 0);
        
            my $start_row_index = int($start / $line_width);
            my $r2 = $start % $line_width;
            $start_row_index-- if ($r2 == 0);
            my $start_loc = $r2-1;
            $start_loc = $line_width - 1 if ($r2 == 0);
            
            if ($start_row_index == $end_row_index) {
                substr($snp_rows[$start_row_index]->[$count_o], $start_loc, $len, "<font color=#993300>" . ("^" x $len). "</font>");
            } else {
                my $len1 = $line_width-$start_loc;
                my $len2 = $len-$len1;
                
                substr($snp_rows[$start_row_index]->[$count_o], $start_loc, $len1, "<font color=#993300>" . ("^" x $len1). "</font>");
                substr($snp_rows[$start_row_index+1]->[$count_o], 0, $len2, "<font color=#993300>" . ("^" x $len2). "</font>");
            }
            $count_o++;
        }
    }
 
    #add SBE or allele-specific primers or sequencing primers
    my $count = 0;
    if ($snp_f_list){
        for (my $i = 0; $i < scalar(@$snp_f_list); $i++) {
            if (@$snp_f_list[$i] ) {
                my @primer_results = split(/\s+/, @$snp_f_list[$i]);
                my ($start, $len, $seq) = ();
                if ($primer_type == 10) {
                    ($start, $len, $seq) = ($primer_results[1],$primer_results[2], $primer_results[8]);
                } else {
                    ($start, $len, $seq) = ($primer_results[3],$primer_results[4], $primer_results[10]);
                }
                
                my $end_row_index = int(($start+$len-1) / $line_width);
                my $r1 = ($start+$len) % $line_width;
        
                $end_row_index-- if ($r1 == 0);
                my $end_loc = $r1;
                $end_loc = $line_width if ($r1 == 0);
        
                my $start_row_index = int($start / $line_width);
                my $r2 = $start % $line_width;
                $start_row_index-- if ($r2 == 0 && $start > 0);
                my $start_loc = $r2;
                $start_loc = $line_width - 1 if ($r2 == 0 && $start > 0);
                
                $count++ if (index($snp_rows[$start_row_index]->[$count], '<') != -1);  # chekc if there is already contents
                
                if ($start_row_index == $end_row_index) {
                    substr($snp_rows[$start_row_index]->[$count], $start_loc, $len, "<font color=#000fff>" . $seq . "</font>");
                } else {
                    my $len1 = $line_width-$start_loc;
                    my $len2 = $len-$len1;
                    substr($snp_rows[$start_row_index]->[$count], $start_loc, $len1, "<font color=#000fff>" . substr($seq, 0, $len1). "</font>");
                    substr($snp_rows[$start_row_index+1]->[$count], 0, $len2, "<font color=#000fff>" . substr($seq, length($seq)-$len2, $len2). "</font>");
                }                
                $count++;
            }
        }
    }
   
    if ($snp_r_list){
        for (my $i = 0; $i < scalar(@$snp_r_list); $i++) {
            if (@$snp_r_list[$i] ) {
                my @primer_results = split(/\s+/, @$snp_r_list[$i]);
                my ($start, $len, $seq) = ();
                if ($primer_type == 10) {
                    ($start, $len, $seq) = ($primer_results[1],$primer_results[2], $primer_results[8]);
                } else {
                    ($start, $len, $seq) = ($primer_results[3],$primer_results[4], $primer_results[10]);
                }
                $start -= $len;
                $seq = reverse($seq);
         
                my $end_row_index = int(($start+$len-1) / $line_width);
                my $r1 = ($start+$len-1) % $line_width;
                $end_row_index-- if ($r1 == 0);
                my $end_loc = $r1;
                $end_loc = $line_width if ($r1 == 0);
        
                my $start_row_index = int($start / $line_width);
                my $r2 = $start % $line_width;
                $start_row_index-- if ($r2 == 0);
                my $start_loc = $r2 - 1;
                $start_loc = $line_width - 1 if ($r2 == 0);
                       
            
                if ($start_row_index == $end_row_index) {
                    substr($snp_rows[$start_row_index]->[$count], $start_loc, $len, "<font color=#993300>" . $seq. "</font>");
                } else {
                    my $len1 = $line_width-$start_loc;
                    my $len2 = $len-$len1;
                    substr($snp_rows[$start_row_index]->[$count], $start_loc, $len1, "<font color=#993300>" . substr($seq, 0, $len1). "</font>");
                    substr($snp_rows[$start_row_index+1]->[$count], 0, $len2, "<font color=#993300>" . substr($seq, length($seq)-$len2, $len2). "</font>");
                }                
            }
            $count++;
        }
    }
  
    #add a mismatch symbol * to allele-specific primers in tetra-primer ARMS PCR or allele-specific primers
    # only have one inner forward primer and one inner reverse primer
    if ($primer_type eq '9') {
        my @target_regions = split(/\s+/,$target);  # multiple target regions are allowed.
        my ($start, $len) = split(/,/,$target_regions[0]);  # only one target
        my $start_row_index = int($start / $line_width);
        my $r2 = $start % $line_width;
        $start_row_index-- if ($r2 == 0);
        my $start_loc = $r2-1;
        $start_loc = $line_width-1 if ($r2 == 0);
        
        # for forward
        my $f_r = 1;
        my $pos = $start_loc - 2;
        my $start_row = $start_row_index;
        $start_row-- if ($start_loc == 0);

        if ($snp_f_list){
            $f_r++ if ($snp_rows[$start_row]->[0] =~ m/<font \S+>[<|>]+<\/font>/g || 
                   $snp_rows[$start_row]->[0] =~ m/^[\-]+$/);
       
            substr($snp_rows[$start_row]->[$f_r], $pos, 1, "<font color=black>*</font>");
        }
    
        # for reverse
        my $r_r = 0;
        if ($snp_r_list){

            $r_r++ if ($snp_rows[$start_row_index]->[0] =~ m/<font \S+>[<|>]+<\/font>/g ||
                       $snp_rows[$start_row_index]->[0] =~ m/^[\-]+$/); 
            $start_loc += &count_existing_offsets($snp_rows[$start_row_index]->[$r_r], "<font color=#000fff>", "</font>");
            $pos =  $start_loc + 2;

            my $temp_pos = length($snp_rows[$start_row_index]->[$r_r]);
            if ($pos > $temp_pos-1) {
                $pos = $pos - ($temp_pos-1) -1;
                $start_row_index++;
                $pos += &count_existing_offsets($snp_rows[$start_row_index]->[$r_r], "<font color=#000fff>", "</font>");
            
            }
            substr($snp_rows[$start_row_index]->[$r_r], $pos, 0, "<font color=black>*</font>");
        }
    }
    elsif (($primer_type eq '7' || $primer_type eq '8') && $SECOND_MISMATCH == 1) {
        my @target_regions = split(/\s+/,$target);  # multiple target regions are allowed.
        my ($start, $len) = split(/,/,$target_regions[0]);  # only one target
        my $start_row_index = int($start / $line_width);
        my $r2 = $start % $line_width;
        $start_row_index-- if ($r2 == 0);
        my $start_loc = $r2-1;
        $start_loc = $line_width-1 if ($r2 == 0);
        
        my $start_loc_bk = $start_loc;
        
        # for forward
        my $f_r = 1;
        my $pos = $start_loc - 2;
        my $start_row = $start_row_index;
        $start_row-- if ($start_loc == 0);

        if ($snp_f_list && @$snp_f_list > 0){
            $f_r++ if ($snp_rows[$start_row]->[0] =~ m/<font \S+>[<|>]+<\/font>/g || 
                   $snp_rows[$start_row]->[0] =~ m/^[\-]+$/);
            $f_r += @$snp_f_list-1; 
            substr($snp_rows[$start_row]->[$f_r], $pos, 1, "<font color=black>*</font>");
        }
    
        # for reverse
        my $r_r = 0;
        if ($snp_r_list && @$snp_r_list > 0){

            $r_r++ if ($snp_rows[$start_row_index]->[0] =~ m/<font \S+>[<|>]+<\/font>/g ||
                       $snp_rows[$start_row_index]->[0] =~ m/^[\-]+$/); 
            $start_loc += &count_existing_offsets($snp_rows[$start_row_index]->[$r_r], "<font color=#000fff>", "</font>");
            $pos =  $start_loc + 2;

            my $temp_pos = length($snp_rows[$start_row_index]->[$r_r]);
            if ($pos > $temp_pos-1) {
                $pos = $pos - ($temp_pos-1) -1;
                $start_row_index++;
                $pos += &count_existing_offsets($snp_rows[$start_row_index]->[$r_r], "<font color=#000fff>", "</font>");
            
            }

            if (!$snp_f_list || @$snp_f_list == 0) {
                # insert a new row into the first row
                my $new_row = "-" x $line_width;
   
                $pos = $start_loc_bk+2; 
                $temp_pos = length($new_row);
                if ($pos > $temp_pos-1) {
                   $pos = $pos - ($temp_pos-1) -1;
                   $start_row_index++;
                }
   
                my $tmp_ref = $snp_rows[$start_row_index];
                my @tmp_arr = ();
                $tmp_arr[0] = $new_row;
                for (my $i = 0; $i < $max_rows; $i++) {
                    $tmp_arr[$i + 1] = $tmp_ref->[$i];
                }
                splice(@tmp_arr, 0, 0, $new_row);
         
                substr($tmp_arr[$r_r], $pos, 0, "<font color=black>*</font>");
                $snp_rows[$start_row_index] = \@tmp_arr;
            } else {
                substr($snp_rows[$start_row_index]->[$r_r], $pos, 0, "<font color=black>*</font>");
            }
        }

    }

    #print a table
    print $file_handle "<p><b>Primer view:</b></p>\n";
   
    print $file_handle "<table><tr><td>\n";  # outer table
    
    print $file_handle "<table class=primerview><tbody>\n";
    for (my $i = 0; $i < scalar(@seq_rows); $i++) {
        print $file_handle "<tr class='even'><td>\n";
        print $file_handle $i * $line_width + 1;
        print $file_handle "</td><td>";
        print $file_handle $seq_rows[$i];
        print $file_handle "</td></tr>\n";
        
        my $rows = $snp_rows[$i];
        for (my $j = 0; $j < scalar(@$rows); $j++) {
            my $line = @$rows[$j];
            next if (!$line);
            next if ($line =~ /^\-+$/);  # if all spaces
            $line =~ s/-/&nbsp;/g;
            print $file_handle "<tr class='odd'><td>\n";
            print $file_handle "&nbsp;";
            print $file_handle "</td><td>";
            print $file_handle $line;
            print $file_handle "</td></tr>\n";
        }
    }
    print $file_handle "</tbody></table>\n";

    print $file_handle "</td><td>\n";
    print $file_handle "<table class=small><tbody>";
    if ($primer_type eq '1' || $primer_type eq '2' || $primer_type eq '5' || $primer_type eq '7' || $primer_type eq '9') {
        print $file_handle "<tr><td>Outer forward primer:</td><td><font family=Courier color=#000fff>>>>>>>>>>></font></td></tr>";
        print $file_handle "<tr><td>Outer reverse primer:</td><td><font family=Courier color=#993300><<<<<<<<<<</font></td></tr>";
    }
    if ($primer_type eq '3') {
        print $file_handle "<tr><td>Oligo primer:</td><td><font family=Courier color=#993300>^^^^^^^^^^</font></td></tr>";
    }
    if ($primer_type >= 5 && $primer_type <= 8) {
        print $file_handle "<tr><td>Forward SBE/allele-specific primer:</td><td><font color=#000fff>Blue</font></td></tr>";
        print $file_handle "<tr><td>Reverse SBE/allele-specific primer:</td><td><font color=#993300>Brown</font></td></tr>";
    }
    if ($primer_type == 9) {
        print $file_handle "<tr><td>Inner forward allele-specific primer:</td><td><font color=#000fff>Blue</font></td></tr>";
        print $file_handle "<tr><td>Inner reverse allele-specific primer:</td><td><font color=#993300>Brown</font></td></tr>";
    }
    
    if ($primer_type eq '10') {
        if ($SEQUENCING_DIRECTION eq "toward3") {
            print $file_handle "<tr><td>Forward sequencing primer:</td><td><font color=#000fff>Blue</font></td></tr>";
        } else {
            print $file_handle "<tr><td>Reverse sequencing primer:</td><td><font color=#993300>Brown</font></td></tr>";
        }
    }
    
    if ($target) {
        if ($primer_type >= 5 && $primer_type <= 10) { 
            print $file_handle "<tr><td>Target:</td><td><font color=red>Red</font></td></tr>";
        }
            
    }
    if ($included) {
        print $file_handle "<tr><td>Included region:</td><td><font color=green>Green</font></td></tr>";
    }
    if ($excluded) {
        print $file_handle "<tr><td>excluded region:</td><td><font color=orange>Orange</font></td></tr>";
    }
    
    if ($primer_type == 9 || ($primer_type == 7 || $primer_type == 8) && $SECOND_MISMATCH == 1) {
        print $file_handle "<tr><td>Mismatch in allele-specific primers:</td><td><font color=black>*</font></td></tr>";
    }
    
    print $file_handle "</tbody></table>\n";
    print $file_handle "</td><tr></table>\n";
}

sub split_sequence {
    my ($seq, $line_width) = @_;
    my $len = length($seq);
    my $rows = int($len/$line_width);
    $rows++ if ($len % $line_width != 0);
    
    my @seq_rows = ();
    for (my $i = 0; $i < $rows; $i++) {
        my $offset = $i * $line_width;
        my $width = $line_width;
        $width = ($len - $offset) if ($i == $rows-1); 
        $seq_rows[$i] = uc(substr($seq, $offset, $line_width));
    }
    
    return @seq_rows;
}
    
# hilight the sequence using colors for included, excluded and target regions
sub highlight_seq_rows {
    my ($seq_rows, $included, $excluded, $target, $line_width) = @_;
    my @included_rows = ();
    my @excluded_rows = ();
    my @target_rows = ();
    for (my $i = 0; $i < scalar(@$seq_rows); $i++) {
        $included_rows[$i] = @$seq_rows[$i];
        $excluded_rows[$i] = @$seq_rows[$i];
        $target_rows[$i] = @$seq_rows[$i];
    }
    
    # default $first_base_index = 1;
    
    if ($included) {  # only one included region is allowed
        my ($start, $len) = split(/\s*,\s*/,$included);
        my $end_index = int(($start+$len-1) / $line_width);
        my $r1 = ($start+$len-1) % $line_width;
        $end_index-- if ($r1 == 0);
        my $end_loc = $r1;
        $end_loc = $line_width if ($r1 == 0);
        substr($included_rows[$end_index], $end_loc, 0, "</font>");
        
        my $start_index = int($start / $line_width);
        my $r2 = $start % $line_width;
        $start_index-- if ($r2 == 0);
        my $start_loc = $r2 - 1;
        $start_loc = $line_width - 1 if ($r2 == 0);
        substr($included_rows[$start_index], $start_loc, 0, "<font color=green>");
        
         if ($start_index != $end_index) {
            for (my $k = $start_index+1; $k <= $end_index; $k++) {
                $included_rows[$k] = "<font color=green>".$included_rows[$k];
            }
            for (my $k = $start_index; $k < $end_index; $k++) {
                $included_rows[$k] .= "</font>";
            }
        }
    }

  
    if ($excluded) {
        my @ex_regions = split(/\s+/,$excluded);   # multiple excluded regions are allowed
        for (my $i = 0; $i < @ex_regions; $i++) {
            my ($start, $len) = split(/,/, $ex_regions[$i]);

            my $start_index = int($start / $line_width);
            my $r2 = $start % $line_width;
            $start_index-- if ($r2 == 0);
            my $start_loc = $r2 - 1;
            $start_loc = $line_width - 1 if ($r2 == 0);
            
            $start_loc += &count_existing_offsets($excluded_rows[$start_index], "<font color=orange>", "</font>");
            substr($excluded_rows[$start_index], $start_loc, 0, "<font color=orange>");

            my $end_index = int(($start+$len-1) / $line_width);
            my $r1 = ($start+$len-1) % $line_width;
            $end_index-- if ($r1 == 0);
            my $end_loc = $r1;
            $end_loc = $line_width if ($r1 == 0);
            
            $end_loc += &count_existing_offsets($excluded_rows[$end_index], "<font color=orange>", "</font>");
            substr($excluded_rows[$end_index], $end_loc, 0, "</font>");
        
        
            if ($start_index != $end_index) {
                for (my $k = $start_index+1; $k <= $end_index; $k++) {
                    $excluded_rows[$k] = "<font color=orange>".$excluded_rows[$k];
                }
                for (my $k = $start_index; $k < $end_index; $k++) {
                    $excluded_rows[$k] .= "</font>";
                }
            }
        }
    }

    if ($target) {
        my @target_regions = split(/\s+/,$target);  # multiple target regions are allowed.
        
        for (my $i = 0; $i < @target_regions; $i++) {
            my ($start, $len) = split(/,/,$target_regions[$i]);

            my $start_index = int($start / $line_width);
            my $r2 = $start % $line_width;
            $start_index-- if ($r2 == 0);
            my $start_loc = $r2-1;
            $start_loc = $line_width-1 if ($r2 == 0);

            $start_loc += &count_existing_offsets($target_rows[$start_index], "<font color=red>", "</font>");
            substr($target_rows[$start_index], $start_loc, 0, "<font color=red>");

            my $end_index = int(($start+$len-1) / $line_width);
            my $r1 = ($start+$len-1) % $line_width;
            $end_index-- if ($r1 == 0);
            my $end_loc = $r1;
            $end_loc = $line_width if ($r1 == 0);

            $end_loc += &count_existing_offsets($target_rows[$end_index], "<font color=red>", "</font>");
            substr($target_rows[$end_index], $end_loc, 0, "</font>");
        
    
            if ($start_index != $end_index) {
                for (my $k = $start_index+1; $k <= $end_index; $k++) {
                    $target_rows[$k] = "<font color=red>".$target_rows[$k];
                }
                for (my $k = $start_index; $k < $end_index; $k++) {
                    $target_rows[$k] .= "</font>";
                }
            }
        }
    }
    
    # merge the highlights together like an alignment
   
    my $len = length(@$seq_rows[0]);
    my @seqs = ();
    for (my $j = 0; $j < scalar(@$seq_rows); $j++) { 
        $seqs[$j] = '';
        my $p1 = 0;
        my $p2 = 0;
        my $p3 = 0;
        for (my $i = 0; $i < length(@$seq_rows[$j]); $i++) {
            if (substr($included_rows[$j], $p1, 1) eq substr($excluded_rows[$j], $p2, 1) &&
                substr($included_rows[$j], $p1, 1) eq substr($target_rows[$j], $p3, 1)) {
                $seqs[$j] .= substr($included_rows[$j], $p1, 1);
                $p1++;
                $p2++;
                $p3++;
            } else {
                if (substr($included_rows[$j], $p1, 1) eq '<') {
                    my $index = index($included_rows[$j], '>', $p1);
                    $seqs[$j] .= substr($included_rows[$j], $p1, $index-$p1+1);
                    $p1 = $index + 1;
                    $i--;
                }
                elsif (substr($excluded_rows[$j], $p2, 1) eq '<') {
                    my $index = index($excluded_rows[$j], '>', $p2);
                    $seqs[$j] .= substr($excluded_rows[$j], $p2, $index-$p2+1);
                    $p2 = $index + 1;
                    $i--;
                }
                
                elsif (substr($target_rows[$j], $p3, 1) eq '<') {
                    my $index = index($target_rows[$j], '>', $p3);
                    $seqs[$j] .= substr($target_rows[$j], $p3, $index-$p3+1);
                    $p3 = $index + 1;
                    $i--;
                }
                else {
                }
            }
         }
        
    }
 
    return @seqs;    
}
    
    
sub count_existing_offsets {
    my ($seq_row, $pattern1, $pattern2) = @_;
    my $len1 = 0;
    my $len2 = 0;

    while ($seq_row =~ m/$pattern1/g) {
        $len1 += length($pattern1);
    }
    
    while ($seq_row =~ m/$pattern2/g) {
        $len2 += length($pattern2);
    }
    
    return $len1+$len2;

}

###############################################################################################################
    

# find the start index of next record of sequence when there are multiple
# sequences in FASTA format
sub find_index {
    my ($seq, $sub_str, $start) = @_;
    my $ind = index($seq, $sub_str, $start);
    while ($seq =~ m/<[^<>]*>/g) {
        my $pos = pos($seq)-1;
        if ($ind != $pos || $ind < $pos) {
            last;
        }
        $ind = index($seq, $sub_str, $ind+1);
    }
    
    return $ind;
}
    
sub find_excluded_index {
    my $seq = shift;
    my $pos = -1;
    while ($seq =~ m/<[^<>]*>/g) {
        $pos = pos($seq)-1;
    }
    return $pos;
}


############# The following functions for pre-analysis of user's input sequences  ##################
    
#pre-analysis of input sequences
sub pre_analysis_of_input_sequences {
    # arrays to store the calculated data
    my @seq_lens = ();
    my @gcs = ();
    my @ns = ();
    my @polymorphisms = (0,0,0,0,0,0); # only for SNP genotyping primers
        
    my $s = $query->param('SEQUENCE');
    my $s0 = $query->param('SEQUENCEFILE');
    if ($s && $s !~ /^\s*$/ || $s0 && $s0 !~ /^\s*$/) {
	my $sequences = "";
	if ($s0 || $s0 =~ /^\s*$/) {  # if no sequence file
            $sequences = $s;
	}

        my $first_sequence = 1;
	my $next_header = '';
	my $last_seq = 0;
 	while(!$last_seq) {
            my $v = '';
            if ($s0 && $s0 !~ /^\s*$/) {  # for sequence file
                my $count = 0;
                $v = $next_header;
                if ($v !~ /^\s*$/) {  #there is header already
                    $count = 1;
                }

		while(<$s0>) {
                    chomp($_);
                    my $index1 = index ($_, '>');  # search the first '>'
                    if ($index1 >= 0) {   # found the next record
			$count++;         # $count should be 2
                        my $index2 = &find_excluded_index ($v . $_);
                        my $index3 = rindex($v . $_, '>');
                        if ($index2 == $index3) {
                            $v .= "$_";
                            $count--;
                            next;
                        }
                        
			if ($count == 2) {
                            $next_header = substr($_, $index1, length($_)-1)."\n";
                            if ($index1 > 0) {
				$v .= substr($_, 0, $index1 -1);
				$v .= "\n";
                            }
                            last;
                        }
		        elsif ($count == 1) {
                            $v .= "$_\n";
                        }
		    }else {
                        $v .= "$_";
                    }
                    
		}
      
                if ($count == 1) {  # last sequence)
                    $last_seq = 1;
                }

	    } elsif ($s && $s !~ /^\s*$/) {   # no sequence file and there is sequencess in the text area)
        	my $index1 = index ($sequences, '>', 0);  # search the first '>'
		my $index2 = &find_index ($sequences, '>', $index1 + 1);  # search the second '>'
 		if ($index2 >= 0) {
                    $v = substr($sequences, $index1, $index2 - $index1);
                    $sequences = substr($sequences, $index2, length($sequences)-1);
		} else {
                    $v = $sequences;
                    $last_seq = 1;
		}
                
            }
 
            #clean sequences
            # remove header in FASTA format
            if ($v =~ /^\s*>([^\n]*)/) {
		# Sequence is in Fasta format.
		$v =~ s/^\s*>([^\n]*)//;
	    }
            # remove preprocessed maskt symbols
            $v =~ s/[\[\]\<\>\{\}]//g;
            
            # remove spaces
            $v =~ s/\s*//g;
 
            &calculate_seq_properties($v, \@seq_lens, \@gcs, \@ns, \@polymorphisms);
        }
    }
    
    if (@seq_lens == 0) {
        print "<p><font color=red>No sequences are available!</font></p>";
        print "<p>Please go back to the homepage.</p>";
        return;
    }
    
    #create a HTML page
    
    print "<h1>Pre-analysis of Input Sequences</h1>\n";    
    print qq{
         <table class="standard">
         <thead>
         <h2>Table 1. Basic statistics of properties of input sequences</h2>
         <tr>\n
         <th>Item</th><th># of Input Sequences</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th><th>Coe. of Variation (%)</th>\n
         </tr></thead>\n
         <tbody>
    };
 
    my $len_bstats = new BStats(\@seq_lens);   
    print "<tr><td>Sequence length</td>";
    print "<td>" . $len_bstats->size() . "</td>";
    print "<td>" . $len_bstats->average() . "</td>";
    print "<td>" . $len_bstats->sd() . "</td>";
    print "<td>" . $len_bstats->min() . "</td>";
    print "<td>" . $len_bstats->max() . "</td>";
    print "<td>" . $len_bstats->cv() . "</td>";
    print "</tr>\n";
    
    my $gc_bstats = new BStats(\@gcs);
    print "<tr><td>GC contents (%)</td>";
    print "<td>" . $gc_bstats->size() . "</td>";
    print "<td>" . $gc_bstats->average() . "</td>";
    print "<td>" . $gc_bstats->sd() . "</td>";
    print "<td>" . $gc_bstats->min() . "</td>";
    print "<td>" . $gc_bstats->max() . "</td>";
    print "<td>" . $gc_bstats->cv() . "</td>";
    print "</tr>\n";
    
    my $n_bstats = new BStats(\@ns);
    print "<tr><td>Ambiguity code (N)</td>";
    print "<td>" . $n_bstats->size() . "</td>";
    print "<td>" . $n_bstats->average() . "</td>";
    print "<td>" . $n_bstats->sd() . "</td>";
    print "<td>" . $n_bstats->min() . "</td>";
    print "<td>" . $n_bstats->max() . "</td>";
    print "<td>" . $n_bstats->cv() . "</td>";
    print "</tr>\n";
    print "</tbody></table>";
    
    print "<table><tr><td>";
    print &len_graph(@seq_lens);
    print "</td></tr></table>";
    print "<p><table><tr><td>";
    print &gc_graph(@gcs);
    print "</td></tr></table>";

 
    if ($primer_type >= 4 && $primer_type <= 9) { 
        print qq{
             <table class="standard">fsen
             <thead>
             <h2>Table 2. Polymorphisms of input sequences</h2>
             <tr>\n
             <th>SNP type</th><th># of Occurences</th><th>%</th></th>\n
             </tr></thead>\n
             <tbody>
        };
 
        my @poly_types = ('S (G/C)', 'W (A/T)', 'R (G/A)', 'Y (T/C)', 'K (G/T)', 'M (A/C)');
        for (my $i = 0; $i < @polymorphisms; $i++) {
            print "<tr><td>$poly_types[$i]</td><td>$polymorphisms[$i] </td><td>";
            print int($polymorphisms[$i]/scalar(@seq_lens)*10000 + 0.5)/100.0 . "</td></tr>\n";
        }
        print "</tbody></table>\n";
    }
 
    print  "$wrapup\n";
  
}

# calculate properties of input sequences
sub calculate_seq_properties{
    my ($seq, $seq_lens, $gcs, $ns, $polymorphisms) = @_;
    $seq = uc($seq);
    my $gc = &oligo_gc($seq);
    $gc = int($gc * 100 + 0.5)/100.0;
    my $n = &find_ambiguity_ns($seq);
    my $len = length($seq);
    
    while ($seq =~ /S/ig){
        @$polymorphisms[0]++;
    }
    while ($seq =~ /W/ig){
        @$polymorphisms[1]++;
    }
    while ($seq =~ /R/ig){
        @$polymorphisms[2]++;
    }
    while ($seq =~ /Y/ig){
        @$polymorphisms[3]++;
    }
    while ($seq =~ /K/ig){
        @$polymorphisms[4]++;
    }
    while ($seq =~ /M/ig){
        @$polymorphisms[5]++;
    }

    push(@$seq_lens, $len);
    push(@$gcs, $gc);
    push(@$ns, $n);
}

# draw histogram of input sequence length
sub len_graph{
    my @seqs = @_;


    my $largest = 0;
    foreach (@seqs) {
        # find the largest sequence
        my $x = scalar($_);
        if($x > $largest){
            $largest = $x;
        }
    }

    # set the number of bars
    my $bar_num = floor($largest / 50);
    if(ceil($largest) > 0){
        $bar_num++;
    }

    # create the array of bars for the graph
    my @data;
    for(my $i = 0; $i <= $bar_num; $i++){
        my $begin = ($i * 50);
        my $bar = ($begin + 1) .'-'. ($begin + 50);
        $data[0][$i] = $bar;
        $data[1][$i] = 0;
    }

    foreach (@seqs) {
        # determine which group this seq belongs to
        my $pos = floor($_ / 50);
        if(ceil($_ / 50) > 0){
            $pos++;
        }
        # increase group by one
        $data[1][$pos]++;
    }

    my $mygraph = GD::Graph::bars->new(800, 200);
    $mygraph->set(
                  x_label     => 'Sequence Length (base pairs)',
                  y_label     => 'Number of Sequences',
                  title       => 'Distribution of Sequence Lengths',
	          long_ticks => 1,
	          #x_label_skip => 2,
	          x_label_position => 1/2,
	          textclr => 'black',
	          labelclr => 'black',
	          axislabelclr => 'black',
	          accent_treshold => 1110,
	          bar_spacing => 8,
	          shadow_depth => 4,
	          shadowclr => 'black',
	          l_margin => 0,
	          r_margin => 0,
	          x_labels_vertical =>1,
	          dclrs => [ qw(dblue pink blue cyan) ]
                  ) or warn $mygraph->error;

    my $myimage = $mygraph->plot(\@data) or die $mygraph->error;
    my $file = time;
    open(IMG, ">$user_result_dir/$file.gif") or die $!;
    binmode IMG;
    print IMG $myimage->gif;
    close IMG;
    #print "<img src=\"/tmp/file.png\">";
    my $img = "<img src=\"$relative_user_result_dir/$file.gif\">";
    return $img;

}

# draw histogram of input sequence GC content
sub gc_graph{
    my @seqs = @_;
    # create the array of bars for the graph
    my @data;
    for(my $i = 0; $i <= 100; $i++){;
        $data[0][$i] = $i;
        $data[1][$i] = 0;
    }

    my $largest = 0;
    foreach (@seqs) {
        # find the largest sequence
        my $perc = floor($_);
        if(ceil($perc) > 0){
            $perc++;
        }
	$data[1][$perc]++;
    }

    my $mygraph = GD::Graph::bars->new(600, 200);
    $mygraph->set(
                  x_label     => 'GC Content (%)',
                  y_label     => 'Number of Sequences',
                  title       => 'Distribution of GC Contents of Sequences',
	          long_ticks => 1,
	          x_label_skip => 5,
	          x_label_position => 1/2,
	          textclr => 'black',
	          labelclr => 'black',
	          axislabelclr => 'black',
	          #bar_spacing => 2,
	          #shadow_depth => 1,
	          #shadowclr => 'black',
	          dclrs => [ qw(dblue pink blue cyan) ]
                  ) or warn $mygraph->error;

    my $myimage = $mygraph->plot(\@data) or die $mygraph->error;
    my $file = (time + 1);
    open(IMG, ">$user_result_dir/$file.gif") or die $!;
    binmode IMG;
    print IMG $myimage->gif;
    close IMG;
    my $img = "<img src=\"$relative_user_result_dir/$file.gif\">";
    return $img;
}


# write log information for user access statistics
sub log {
    my ($ip,$seqs, $primer_type_str) = @_;
    my $date = strftime( "%A, %B %d, %Y", localtime(time()));
    
    my $log_file = "batchprimer3_log.txt";

    if (-e $log_file) {
	open (FILE, ">>$log_file") or die ("Can't open $!");
    } else {
	open (FILE, ">$log_file") or die ("Can't open $!");
    }
    
    flock(FILE, LOCK_EX) or die ("Can't get exclusive lock: $!");
	
    print (FILE "$date\t$ip\t$seqs\t$primer_type_str\n");
    flock (FILE, LOCK_UN) or die ("Can't  unlock file $log_file");
    close (FILE);
}

# create a zipped file for a job, which can be downloaded from a user
sub create_zip_file {
    my ($dir_name, $dir) = @_;  # $dir is a full path   
    use Archive::Zip;
    my $zip = Archive::Zip->new();
    # add all readable files and directories below . as xyz/*
    $zip->addTree( "$dir", "$dir_name");
    $zip->writeToFileNamed("$dir" . ".zip");
}

############# the following functions to clean and remove the user result directory ##################

# remove expired result directories. All the results are save in the directory $RESULT_DIR
sub remove_expired_results {
    return unless -e $RESULT_DIR;  # nothing to clean before the first result is ever saved
    opendir (DIR, $RESULT_DIR) or die ("Can't open the directory $RESULT_DIR! Please check the result directory and change mode to 777.<p>");

    my @files = readdir(DIR);
    foreach my $file (@files) {
        next if ($file eq '.');
        next if ($file eq '..');
        my $days = (-C "$RESULT_DIR/$file");
        if ($days > 7) {  # keep result in 7 days
            if (-d "$RESULT_DIR/$file") {
                # remove the directory
                &clean_directory("$RESULT_DIR/$file");
                rmdir("$RESULT_DIR/$file");
            } else {
                unlink("$RESULT_DIR/$file");
            }
        }
    }
}

# delete all the files in the directory and sub directory
sub clean_directory {
    my $dir = shift;
    opendir(DIR,$dir);
    my @files = readdir(DIR);
    foreach my $file (@files) {
        if (-d "$dir/$file") {
            &clean_directory("$dir/$file");
        } else {
            unlink("$dir/$file");
        }
    }
    closedir(DIR);
}

############# The following function is to save user's primer parameters  ##################

# save allparameters and options in a file for each user
# and keep the most recent parameters for a primer type
sub save_primer_parameters {
    my ($query, $names) = @_;
    if (! -e $PARAM_DIR) {
        if (!mkdir($PARAM_DIR)) {
            die ("<p>Can't create the parameter directory. Please check if your cgi-bin directory is correct. " .
                 "Probably you need to change the mode of the directory: chmod 777</p>");
        }
        chmod (0777, $PARAM_DIR);
    }
    my $file_name = "$PARAM_DIR/$ip" . "_" . $primer_type . ".txt";
    open (FILE, ">$file_name")
        or die ("<p>Can't not open the file $file_name!  Please check if your parameter directory is correct. " .
                "Probably you need to change the mode of the directory: chmod 777</p>");
                                       
    for (my $i = 0; $i < @$names; $i++) {
        my $key = @$names[$i];
        next if $key eq  'SEQUENCE';      # do not save the pasted sequences
        my $param = $query->param($key);
        print FILE "$key=$param\n";
    }
    close FILE;
}
                                       
