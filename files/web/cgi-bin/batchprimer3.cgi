#!/usr/bin/env perl

########################################################################################################
# Copyright Notice and Disclaimer for BatchPrimer3
#
# Copyright (c) 2007, Depatrtment of Plant Sciences, University of California,
# and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#
# Authors: Dr. Frank M. You.
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

use Carp;
use CGI qw(:standard);
use warnings;
use strict;


########################################################################################################
# You may wish to modify the following variables to suit your installation.

my $MAINTAINER
    ='<a href="mailto:frankyou@pw.usda.gov">frankyou@pw.usda.gov</a>';

my $BATCHPRIMER3_CGI_HOME = "__BP3_CGI__";
my $BATCHPRIMER3_HTDOC_HOME = "http://127.0.0.1:__BP3_PORT__";
my $PARAM_DIR = "$BATCHPRIMER3_CGI_HOME/parameters";

# This is the select box to which you should
# add options for mispriming libraries.  Modify the
# @REPEAT_LIBRARIES variable in batchprimer3_results.cgi as well

my @REPEAT_LIBRARIES = (
        "NONE",
        "Arabidopsis",  
        "Brachypodium",
        "Brassica",
        "Glycine",
        "Gramineae",
        "Hordeum",
        "Lotus",
        "Oryza",
        "Sorghum",
        "Solanum",
        "Triticum",
        "Zea",
        "Wheat");

my @primer_types = (
        "Generic primers",                       # 1
        "SSR screening and primers",             # 2
        "Hybridization oligos",                  # 3
        "SNP (allele) flanking primers",         # 4
        "Single base extension (SBE) primers and SNP (allele) flanking primers", # 5
        "Single base extension (SBE) primers",                                   # 6
        "Allele-specific primers and allele-flanking primers",                   # 7
        "Allele-specific primers",                                               # 8
        "Tetra-primer ARMS-PCR primers",                                         # 9
        "Sequencing primers");                                                   # 10
my @primer_type_defs = (
        "Design pairs of generic primers for any DNA sequences.",    #1
        "Screen SSR motifs and then design pairs of primers that flank a SSR motif. Users can define their own criteria to detect SSR motifs.",  #2
        "Design hybridization oligo primers for any DNA sequences.", #3
        "Design pairs of generic primers that flank SNP sites.",     #4
        "Design pairs of generic primers that flank SNP sites and SNP primers that neighbor SNP nucleotides.", #5
        "Design a SNP primer that is adjacent to a SNP nuecleotide.",                                          #6
        "Design a SNP primer that is adjacent to a SNP nuecleotide.",                                          #7
        "Design allele-specific primers that include two alleles at the 3' end.",                              #8
        "Design tetra primer ARMS PCR primers.",                                                               #9
        "Design DNA sequencing primers");                                                                      #10

# Documentation for this screen.
my $HTDOC = $BATCHPRIMER3_HTDOC_HOME;
my $DOC_URL = "$HTDOC/batchprimer_help.html";
my $BATCHPRIMER3_HOME_URL = "$HTDOC/index.html";

# URL that will process the data from this screen.
my $PROCESS_INPUT_URL = 'batchprimer3_results.cgi';
my $INPUT_FORM_URL = 'batchprimer3.cgi';


# If you make any substantial modifications give this code a new
# version designation.
my $CGI_VERSION = "BatchPrimer3 v1.0, October, 2007";


########################################################################################################

my $ip = $ENV{'REMOTE_ADDR'};
my $user_id = $ip;


my $CLEAR_SEQ = "clear_seq";
my $CLEAR_FORM = "yes";

BEGIN{
    # Ensure that errors go to the web browser.
#    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
}

my $query = new CGI;
my $primer_type = 1;
if ($query->param('PRIMER_TYPE')) {
    $primer_type = $query->param('PRIMER_TYPE');
}



my $clear_form = 0;
my $clear_seq = 0;
if ($query->param('CLEAR_FORM') && $query->param('CLEAR_FORM') eq $CLEAR_FORM) {
    $clear_form = 1;
} elsif ($query->param('CLEAR_FORM') && $query->param('CLEAR_FORM') eq $CLEAR_SEQ) {
    $clear_seq = 1;
}


#set html page title
print $query->header;

my $jscript = qq {
    function change_type(form) {
        document.location.href="$INPUT_FORM_URL?PRIMER_TYPE=" +
        document.form.PRIMER_TYPE.options[document.form.PRIMER_TYPE.selectedIndex].value +
        "&CLEAR_FORM=no";
    }
};

print $query->start_html(
    -title => "BatchPrimer3: a high-throughput web application for picking PCR and sequencing primers ($CGI_VERSION)",
    -style => { -src => "$HTDOC/style_batchprimer.css" },
    -script => $jscript
); 

my %param_hash = &read_user_primer_parameters();

# initial parameters
my $LIGHT_COLOR="#CCCCFF";
my $OLIGO_INPUT_SIZE=30;
my $SOURCE_SEQUENCE_WIDTH = 3 * ($OLIGO_INPUT_SIZE + 1) + 2;
my $SMALL_TEXT=15;
my $VSMALL_TEXT=10;

my $HTML_TABLE_WIDTH = 700;
my $PR_DEFAULT_PRODUCT_MIN_SIZE = 500;
my $PR_DEFAULT_PRODUCT_OPT_SIZE = 700;
my $PR_DEFAULT_PRODUCT_MAX_SIZE = 1000;

# repeat library selection
my $SELECT_SEQ_LIBRARY = "<select name=PRIMER_MISPRIMING_LIBRARY>\n";
for (my $i = 0; $i < @REPEAT_LIBRARIES; $i++) {
    if (exists($param_hash{'PRIMER_MISPRIMING_LIBRARY'}) && $param_hash{'PRIMER_MISPRIMING_LIBRARY'} eq $REPEAT_LIBRARIES[$i]) {
        $SELECT_SEQ_LIBRARY .= "<option selected>" . $REPEAT_LIBRARIES[$i] . "\n";
    } else {
        $SELECT_SEQ_LIBRARY .= "<option>" . $REPEAT_LIBRARIES[$i] . "\n";
    }
}
$SELECT_SEQ_LIBRARY .=  "</select><br>\n";


my $PRIMER_SEQUENCE_ID                       = "";
my $SEQUENCE                                 = exists($param_hash{'SEQUENCE'}) && $clear_form == 0 && $clear_seq == 0 ? $param_hash{'SEQUENCE'} : "";
my $INCLUDED_REGION                          = "";
my $TARGET                                   = "";
my $EXCLUDED_REGION                          = "";
my $PRIMER_SEQUENCE_QUALITY                  = "";
my $PRIMER_LEFT_INPUT                        = exists($param_hash{'PRIMER_LEFT_INPUT'}) && $clear_form == 0 ? $param_hash{'PRIMER_LEFT_INPUT'} : "";
my $PRIMER_RIGHT_INPUT                       = exists($param_hash{'PRIMER_RIGHT_INPUT'}) && $clear_form == 0 ? $param_hash{'PRIMER_RIGHT_INPUT'} : "";

my $SEQUENCEFILE                             = exists($param_hash{'SEQUENCEFILE'}) && $clear_form == 0 ? $param_hash{'SEQUENCEFILE'} : "";
my $EMAIL                                    = exists($param_hash{'EMAIL'}) && $clear_form == 0 ? $param_hash{'EMAIL'} : "";


my $PRIMER_MAX_MISPRIMING                    = exists($param_hash{'PRIMER_MAX_MISPRIMING'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_MISPRIMING'} : "12.00";
my $PRIMER_PAIR_MAX_MISPRIMING               = exists($param_hash{'PRIMER_PAIR_MAX_MISPRIMING'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_MAX_MISPRIMING'} : "24.00";
my $PRIMER_GC_CLAMP                          = exists($param_hash{'PRIMER_GC_CLAMP'}) && $clear_form == 0 ? $param_hash{'PRIMER_GC_CLAMP'} : "0";
my $PRIMER_OPT_SIZE                          = exists($param_hash{'PRIMER_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_SIZE'} : "20";
my $PRIMER_MIN_SIZE                          = exists($param_hash{'PRIMER_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_SIZE'} : "18";
my $PRIMER_MAX_SIZE                          = exists($param_hash{'PRIMER_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_SIZE'} : "27";
my $PRIMER_OPT_TM                            = exists($param_hash{'PRIMER_OPT_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_TM'} : "60.0";
my $PRIMER_MIN_TM                            = exists($param_hash{'PRIMER_MIN_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_TM'} : "57.0";
my $PRIMER_MAX_TM                            = exists($param_hash{'PRIMER_MAX_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_TM'} : "63.0";
my $PRIMER_MAX_DIFF_TM                       = exists($param_hash{'PRIMER_MAX_DIFF_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_DIFF_TM'} : "10.0";
my $PRIMER_MIN_GC                            = exists($param_hash{'PRIMER_MIN_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_GC'} : "20.0";
my $PRIMER_OPT_GC_PERCENT                    = exists($param_hash{'PRIMER_OPT_GC_PERCENT'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_GC_PERCENT'} : "";
my $PRIMER_MAX_GC                            = exists($param_hash{'PRIMER_MAX_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_GC'} : "80.0";
my $PRIMER_SALT_CONC                         = exists($param_hash{'PRIMER_SALT_CONC'}) && $clear_form == 0 ? $param_hash{'PRIMER_SALT_CONC'} : "50.0";
my $PRIMER_DNA_CONC                          = exists($param_hash{'PRIMER_DNA_CONC'}) && $clear_form == 0 ? $param_hash{'PRIMER_DNA_CONC'} : "50.0";
my $PRIMER_NUM_NS_ACCEPTED                   = exists($param_hash{'PRIMER_NUM_NS_ACCEPTED'}) && $clear_form == 0 ? $param_hash{'PRIMER_NUM_NS_ACCEPTED'} : "0";
my $PRIMER_SELF_ANY                          = exists($param_hash{'PRIMER_SELF_ANY'}) && $clear_form == 0 ? $param_hash{'PRIMER_SELF_ANY'} : "8.00";
my $PRIMER_SELF_END                          = exists($param_hash{'PRIMER_SELF_END'}) && $clear_form == 0 ? $param_hash{'PRIMER_SELF_END'} : "3.00";
my $PRIMER_MAX_POLY_X                        = exists($param_hash{'PRIMER_MAX_POLY_X'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_POLY_X'} : "5";
my $PRIMER_NUM_RETURN                        = exists($param_hash{'PRIMER_NUM_RETURN'}) && $clear_form == 0 ? $param_hash{'PRIMER_NUM_RETURN'} : "1";
my $PRIMER_FIRST_BASE_INDEX                  = exists($param_hash{'PRIMER_FIRST_BASE_INDEX'}) && $clear_form == 0 ? $param_hash{'PRIMER_FIRST_BASE_INDEX'} : "1";
my $PRIMER_MIN_QUALITY                       = exists($param_hash{'PRIMER_MIN_QUALITY'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_QUALITY'} : "0";
my $PRIMER_MIN_END_QUALITY                   = exists($param_hash{'PRIMER_MIN_END_QUALITY'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_END_QUALITY'} : "0";
my $PRIMER_QUALITY_RANGE_MIN                 = exists($param_hash{'PRIMER_QUALITY_RANGE_MIN'}) && $clear_form == 0 ? $param_hash{'PRIMER_QUALITY_RANGE_MIN'} : "0";
my $PRIMER_QUALITY_RANGE_MAX                 = exists($param_hash{'PRIMER_QUALITY_RANGE_MAX'}) && $clear_form == 0 ? $param_hash{'PRIMER_QUALITY_RANGE_MAX'} : "100";
my $PRIMER_INSIDE_PENALTY                    = exists($param_hash{'PRIMER_INSIDE_PENALTY'}) && $clear_form == 0 ? $param_hash{'PRIMER_INSIDE_PENALTY'} : "";
my $PRIMER_OUTSIDE_PENALTY                   = exists($param_hash{'PRIMER_OUTSIDE_PENALTY'}) && $clear_form == 0 ? $param_hash{'PRIMER_OUTSIDE_PENALTY'} : "0";
my $PR_DEFAULT_MAX_END_STABILITY             = exists($param_hash{'PR_DEFAULT_MAX_END_STABILITY'}) && $clear_form == 0 ? $param_hash{'PR_DEFAULT_MAX_END_STABILITY'} : "9.0";
my $PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION    = exists($param_hash{'PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION'} : "";
my $PRIMER_INTERNAL_OLIGO_INPUT              = exists($param_hash{'PRIMER_INTERNAL_OLIGO_INPUT'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_INPUT'} : "";
my $PRIMER_INTERNAL_OLIGO_OPT_SIZE           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_OPT_SIZE'} : "20";
my $PRIMER_INTERNAL_OLIGO_MIN_SIZE           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MIN_SIZE'} : "18";
my $PRIMER_INTERNAL_OLIGO_MAX_SIZE           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MAX_SIZE'} : "27";
my $PRIMER_INTERNAL_OLIGO_OPT_TM             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_OPT_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_OPT_TM'} : "60.0";
my $PRIMER_INTERNAL_OLIGO_MIN_TM             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MIN_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MIN_TM'} : "57.0";
my $PRIMER_INTERNAL_OLIGO_MAX_TM             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MAX_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MAX_TM'} : "63.0";
my $PRIMER_INTERNAL_OLIGO_MIN_GC             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MIN_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MIN_GC'} : "20.0";
my $PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT     = exists($param_hash{'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT'} : "";
my $PRIMER_INTERNAL_OLIGO_MAX_GC             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MAX_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MAX_GC'} : "80.0";
my $PRIMER_INTERNAL_OLIGO_SALT_CONC          = exists($param_hash{'PRIMER_INTERNAL_OLIGO_SALT_CONC'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_SALT_CONC'} : "50.0";
my $PRIMER_INTERNAL_OLIGO_DNA_CONC           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_DNA_CONC'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_DNA_CONC'} : "50.0";
my $PRIMER_INTERNAL_OLIGO_SELF_ANY           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_SELF_ANY'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_SELF_ANY'} : "12.00";
my $PRIMER_INTERNAL_OLIGO_MAX_POLY_X         = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MAX_POLY_X'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MAX_POLY_X'} : "5";
my $PRIMER_INTERNAL_OLIGO_SELF_END           = exists($param_hash{'PRIMER_INTERNAL_OLIGO_SELF_END'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_SELF_END'} : "12.00";
my $PRIMER_INTERNAL_OLIGO_MAX_MISHYB         = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MAX_MISHYB'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MAX_MISHYB'} : "12.00";
my $PRIMER_INTERNAL_OLIGO_MIN_QUALITY        = exists($param_hash{'PRIMER_INTERNAL_OLIGO_MIN_QUALITY'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_MIN_QUALITY'} : "0";
my $PRIMER_INTERNAL_OLIGO_NUM_NS             = exists($param_hash{'PRIMER_INTERNAL_OLIGO_NUM_NS'}) && $clear_form == 0 ? $param_hash{'PRIMER_INTERNAL_OLIGO_NUM_NS'} : "0";


my $MUST_XLATE_PRODUCT_MIN_SIZE              = exists($param_hash{'MUST_XLATE_PRODUCT_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'MUST_XLATE_PRODUCT_MIN_SIZE'} : $PR_DEFAULT_PRODUCT_MIN_SIZE;
my $MUST_XLATE_PRODUCT_MAX_SIZE              = exists($param_hash{'MUST_XLATE_PRODUCT_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'MUST_XLATE_PRODUCT_MAX_SIZE'} : $PR_DEFAULT_PRODUCT_MAX_SIZE;
my $PRIMER_PRODUCT_OPT_SIZE                  = exists($param_hash{'PRIMER_PRODUCT_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_PRODUCT_OPT_SIZE'} : $PR_DEFAULT_PRODUCT_OPT_SIZE;
#my $PRIMER_PRODUCT_MIN_TM                    = "-1000000.0";
#my $PRIMER_PRODUCT_OPT_TM                    = "0.0";
#$my $PRIMER_PRODUCT_MAX_TM                    = "1000000.0";

my $PRIMER_WT_TM_GT                          = exists($param_hash{'PRIMER_WT_TM_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_TM_GT'} : "1.0";
my $PRIMER_WT_TM_LT                          = exists($param_hash{'PRIMER_WT_TM_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_TM_LT'} : "1.0";
my $PRIMER_WT_SIZE_LT                        = exists($param_hash{'PRIMER_WT_SIZE_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_SIZE_LT'} : "1.0";
my $PRIMER_WT_SIZE_GT                        = exists($param_hash{'PRIMER_WT_SIZE_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_SIZE_GT'} : "1.0";
my $PRIMER_WT_GC_PERCENT_LT                  = exists($param_hash{'PRIMER_WT_GC_PERCENT_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_GC_PERCENT_LT'} : "0.0";
my $PRIMER_WT_GC_PERCENT_GT                  = exists($param_hash{'PRIMER_WT_GC_PERCENT_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_GC_PERCENT_GT'} : "0.0";
my $PRIMER_WT_COMPL_ANY                      = exists($param_hash{'PRIMER_WT_COMPL_ANY'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_COMPL_ANY'} : "0.0";
my $PRIMER_WT_COMPL_END                      = exists($param_hash{'PRIMER_WT_COMPL_END'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_COMPL_END'} : "0.0";
my $PRIMER_WT_NUM_NS                         = exists($param_hash{'PRIMER_WT_NUM_NS'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_NUM_NS'} : "0.0";
my $PRIMER_WT_REP_SIM                        = exists($param_hash{'PRIMER_WT_REP_SIM'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_REP_SIM'} : "0.0";
my $PRIMER_WT_SEQ_QUAL                       = exists($param_hash{'PRIMER_WT_SEQ_QUAL'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_SEQ_QUAL'} : "0.0";
my $PRIMER_WT_END_QUAL                       = exists($param_hash{'PRIMER_WT_END_QUAL'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_END_QUAL'} : "0.0";
my $PRIMER_WT_POS_PENALTY                    = exists($param_hash{'PRIMER_WT_POS_PENALTY'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_POS_PENALTY'} : "0.0";
my $PRIMER_WT_END_STABILITY                  = exists($param_hash{'PRIMER_WT_END_STABILITY'}) && $clear_form == 0 ? $param_hash{'PRIMER_WT_END_STABILITY'} : "0.0";

my $PRIMER_IO_WT_TM_GT                       = exists($param_hash{'PRIMER_IO_WT_TM_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_TM_GT'} : "1.0";
my $PRIMER_IO_WT_TM_LT                       = exists($param_hash{'PRIMER_IO_WT_TM_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_TM_LT'} : "1.0";
my $PRIMER_IO_WT_SIZE_LT                     = exists($param_hash{'PRIMER_IO_WT_SIZE_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_SIZE_LT'} : "1.0";
my $PRIMER_IO_WT_SIZE_GT                     = exists($param_hash{'PRIMER_IO_WT_SIZE_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_SIZE_GT'} : "1.0";
my $PRIMER_IO_WT_GC_PERCENT_LT               = exists($param_hash{'PRIMER_IO_WT_GC_PERCENT_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_GC_PERCENT_LT'} : "0.0";
my $PRIMER_IO_WT_GC_PERCENT_GT               = exists($param_hash{'PRIMER_IO_WT_GC_PERCENT_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_GC_PERCENT_GT'} : "0.0";
my $PRIMER_IO_WT_COMPL_ANY                   = exists($param_hash{'PRIMER_IO_WT_COMPL_ANY'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_COMPL_ANY'} : "0.0";
my $PRIMER_IO_WT_NUM_NS                      = exists($param_hash{'PRIMER_IO_WT_NUM_NS'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_NUM_NS'} : "0.0";
my $PRIMER_IO_WT_REP_SIM                     = exists($param_hash{'PRIMER_IO_WT_REP_SIM'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_REP_SIM'} : "0.0";
my $PRIMER_IO_WT_SEQ_QUAL                    = exists($param_hash{'PRIMER_IO_WT_SEQ_QUAL'}) && $clear_form == 0 ? $param_hash{'PRIMER_IO_WT_SEQ_QUAL'} : "0.0";

my $PRIMER_PAIR_WT_PR_PENALTY                = exists($param_hash{'PRIMER_PAIR_WT_PR_PENALTY'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_PR_PENALTY'} : "1.0";
my $PRIMER_PAIR_WT_IO_PENALTY                = exists($param_hash{'PRIMER_PAIR_WT_IO_PENALTY'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_IO_PENALTY'} : "0.0";

my $PRIMER_PAIR_WT_DIFF_TM                   = exists($param_hash{'PRIMER_PAIR_WT_DIFF_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_DIFF_TM'} : "0.0";
my $PRIMER_PAIR_WT_COMPL_ANY                 = exists($param_hash{'PRIMER_PAIR_WT_COMPL_ANY'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_COMPL_ANY'} : "0.0";
my $PRIMER_PAIR_WT_COMPL_END                 = exists($param_hash{'PRIMER_PAIR_WT_COMPL_END'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_COMPL_END'} : "0.0";
my $PRIMER_PAIR_WT_PRODUCT_TM_LT             = exists($param_hash{'PRIMER_PAIR_WT_PRODUCT_TM_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_PRODUCT_TM_LT'} : "0.0";
my $PRIMER_PAIR_WT_PRODUCT_TM_GT             = exists($param_hash{'PRIMER_PAIR_WT_PRODUCT_TM_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_PRODUCT_TM_GT'} : "0.0";

my $PRIMER_PAIR_WT_PRODUCT_SIZE_GT           = exists($param_hash{'PRIMER_PAIR_WT_PRODUCT_SIZE_GT'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_PRODUCT_SIZE_GT'} : "0.05";
my $PRIMER_PAIR_WT_PRODUCT_SIZE_LT           = exists($param_hash{'PRIMER_PAIR_WT_PRODUCT_SIZE_LT'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_PRODUCT_SIZE_LT'} : "0.05";

my $PRIMER_PAIR_WT_REP_SIM                   = exists($param_hash{'PRIMER_PAIR_WT_REP_SIM'}) && $clear_form == 0 ? $param_hash{'PRIMER_PAIR_WT_REP_SIM'} : "0.0";

# for SSR screening and primer picking
my $DINUCLEOTIDE_SSR_REPEATS_INPUT           = exists($param_hash{'DINUCLEOTIDE_SSR_REPEATS_INPUT'}) && $clear_form == 0 ? $param_hash{'DINUCLEOTIDE_SSR_REPEATS_INPUT'} : "6";
my $TRINUCLEOTIDE_SSR_REPEATS_INPUT          = exists($param_hash{'TRINUCLEOTIDE_SSR_REPEATS_INPUT'}) && $clear_form == 0 ? $param_hash{'TRINUCLEOTIDE_SSR_REPEATS_INPUT'} : "4";
my $TETRANUCLEOTIDE_SSR_REPEATS_INPUT        = exists($param_hash{'TETRANUCLEOTIDE_SSR_REPEATS_INPUT'}) && $clear_form == 0 ? $param_hash{'TETRANUCLEOTIDE_SSR_REPEATS_INPUT'} : "3";
my $PENTANUCLEOTIDE_SSR_REPEATS_INPUT        = exists($param_hash{'PENTANUCLEOTIDE_SSR_REPEATS_INPUT'}) && $clear_form == 0 ? $param_hash{'PENTANUCLEOTIDE_SSR_REPEATS_INPUT'} : "3";
my $HEXANUCLEOTIDE_SSR_REPEATS_INPUT         = exists($param_hash{'HEXANUCLEOTIDE_SSR_REPEATS_INPUT'}) && $clear_form == 0 ? $param_hash{'HEXANUCLEOTIDE_SSR_REPEATS_INPUT'} : "3";


# for SNP primer or allele-specific primers
my $SNP_PRIMER_MIN_GC                        = exists($param_hash{'SNP_PRIMER_MIN_GC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_GC'} : 20.0;
my $SNP_PRIMER_MAX_GC                        = exists($param_hash{'SNP_PRIMER_MAX_GC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_GC'} : 80.0;
my $SNP_PRIMER_MIN_TM                        = exists($param_hash{'SNP_PRIMER_MIN_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_TM'} : 57.0;
my $SNP_PRIMER_OPT_TM                        = exists($param_hash{'SNP_PRIMER_OPT_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_OPT_TM'} : 60.0;
my $SNP_PRIMER_MAX_TM                        = exists($param_hash{'SNP_PRIMER_MAX_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_TM'} : 63.0;
my $SNP_PRIMER_MIN_SIZE                      = exists($param_hash{'SNP_PRIMER_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_SIZE'} : 15;
my $SNP_PRIMER_OPT_SIZE                      = exists($param_hash{'SNP_PRIMER_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_OPT_SIZE'} : 20;
my $SNP_PRIMER_MAX_SIZE                      = exists($param_hash{'SNP_PRIMER_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_SIZE'} : 30;
my $SNP_INNER_PRODUCT_MIN_SIZE               = exists($param_hash{'SNP_INNER_PRODUCT_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MIN_SIZE'} : 100;
my $SNP_INNER_PRODUCT_OPT_SIZE               = exists($param_hash{'SNP_INNER_PRODUCT_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_OPT_SIZE'} : 200;
my $SNP_INNER_PRODUCT_MAX_SIZE               = exists($param_hash{'SNP_INNER_PRODUCT_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MAX_SIZE'} : 300;
my $SNP_INNER_PRODUCT_MIN_DIFF               = exists($param_hash{'SNP_INNER_PRODUCT_MIN_DIFF'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MIN_DIFF'} : 1.1;
my $SNP_INNER_PRODUCT_MAX_DIFF               = exists($param_hash{'SNP_INNER_PRODUCT_MAX_DIFF'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MAX_DIFF'} : 1.5;
my $SNP_PRIMER_MAX_N                         = exists($param_hash{'SNP_PRIMER_MAX_N'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_N'} : 0;
my $SNP_PRIMER_SALT_CONC                     = exists($param_hash{'SNP_PRIMER_SALT_CONC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_SALT_CONC'} : 50;
my $SNP_PRIMER_DNA_CONC                      = exists($param_hash{'SNP_PRIMER_DNA_CONC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_DNA_CONC'} : 50;
my $SNP_MAX_SELF_COMPLEMENTARITY             = exists($param_hash{'SNP_MAX_SELF_COMPLEMENTARITY'}) && $clear_form == 0 ? $param_hash{'SNP_MAX_SELF_COMPLEMENTARITY'} : 8;
my $SNP_MAX_3_SELF_COMPLEMENTARITY           = exists($param_hash{'SNP_MAX_3_SELF_COMPLEMENTARITY'}) && $clear_form == 0 ? $param_hash{'SNP_MAX_3_SELF_COMPLEMENTARITY'} : 3;

my $SECOND_MISMATCH_POS                      = exists($param_hash{'SECOND_MISMATCH_POS'}) && $clear_form == 0 ? $param_hash{'SECOND_MISMATCH_POS'} : -3;

&input_screen($query);

##################### end of the main program ######################




sub input_screen {
    my ($query) = @_;

 
    my $hr = qq{
       <hr height=3 aligh=left width=$HTML_TABLE_WIDTH> <p>
    };


   # $SEQUENCE = $param_hash{'SEQUENCE'} if exists($param_hash{'SEQUENCE'}); 


    if ($primer_type eq '2') { # SSR primer
        $PRIMER_OPT_SIZE = exists($param_hash{'PRIMER_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_SIZE'} : "21";
	$PRIMER_MIN_SIZE = exists($param_hash{'PRIMER_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_SIZE'} : "18";
	$PRIMER_MAX_SIZE = exists($param_hash{'PRIMER_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_SIZE'} : "23";
        $PRIMER_OPT_TM   = exists($param_hash{'PRIMER_OPT_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_TM'} : "55.0";
	$PRIMER_MIN_TM   = exists($param_hash{'PRIMER_MIN_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_TM'} : "50.0";
	$PRIMER_MAX_TM   = exists($param_hash{'PRIMER_MAX_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_TM'} : "70.0";
        $PRIMER_MAX_DIFF_TM = exists($param_hash{'PRIMER_MAX_DIFF_TM'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_DIFF_TM'} : "20.0";
 	$PRIMER_MIN_GC   = exists($param_hash{'PRIMER_MIN_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_MIN_GC'} : "30.0";
	$PRIMER_OPT_GC_PERCENT  = exists($param_hash{'PRIMER_OPT_GC_PERCENT'}) && $clear_form == 0 ? $param_hash{'PRIMER_OPT_GC_PERCENT'} : "50.0";
	$PRIMER_MAX_GC = exists($param_hash{'PRIMER_MAX_GC'}) && $clear_form == 0 ? $param_hash{'PRIMER_MAX_GC'} : "70.0";
	$MUST_XLATE_PRODUCT_MIN_SIZE = exists($param_hash{'MUST_XLATE_PRODUCT_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'MUST_XLATE_PRODUCT_MIN_SIZE'} : 100;
	$PRIMER_PRODUCT_OPT_SIZE = exists($param_hash{'PRIMER_PRODUCT_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'PRIMER_PRODUCT_OPT_SIZE'} : 150;
	$MUST_XLATE_PRODUCT_MAX_SIZE = exists($param_hash{'MUST_XLATE_PRODUCT_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'MUST_XLATE_PRODUCT_MAX_SIZE'} : 300;
    }
    
    if ($primer_type eq '9') { # tetra-primer ARMS PCR
        $SNP_PRIMER_MIN_GC   = exists($param_hash{'SNP_PRIMER_MIN_GC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_GC'} : 20.0;
        $SNP_PRIMER_MAX_GC   = exists($param_hash{'SNP_PRIMER_MAX_GC'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_GC'} : 80.0;
        $SNP_PRIMER_MIN_TM   = exists($param_hash{'SNP_PRIMER_MIN_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_TM'} : 50.0;
        $SNP_PRIMER_OPT_TM   = exists($param_hash{'SNP_PRIMER_OPT_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_OPT_TM'} : 65.0;
        $SNP_PRIMER_MAX_TM   = exists($param_hash{'SNP_PRIMER_MAX_TM'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_TM'} : 80.0;
        $SNP_PRIMER_MIN_SIZE = exists($param_hash{'SNP_PRIMER_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_SIZE'} : 26;
        $SNP_PRIMER_OPT_SIZE = exists($param_hash{'SNP_PRIMER_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_OPT_SIZE'} : 28;
        $SNP_PRIMER_MAX_SIZE = exists($param_hash{'SNP_PRIMER_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_SIZE'} : 30;
        $SNP_INNER_PRODUCT_MIN_SIZE = exists($param_hash{'SNP_INNER_PRODUCT_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MIN_SIZE'} : 100;
        $SNP_INNER_PRODUCT_OPT_SIZE = exists($param_hash{'SNP_INNER_PRODUCT_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_OPT_SIZE'} : 200;
        $SNP_INNER_PRODUCT_MAX_SIZE = exists($param_hash{'SNP_INNER_PRODUCT_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MAX_SIZE'} : 300;
        $SNP_INNER_PRODUCT_MIN_DIFF = exists($param_hash{'SNP_INNER_PRODUCT_MIN_DIFF'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MIN_DIFF'} : 1.1;
        $SNP_INNER_PRODUCT_MAX_DIFF = exists($param_hash{'SNP_INNER_PRODUCT_MAX_DIFF'}) && $clear_form == 0 ? $param_hash{'SNP_INNER_PRODUCT_MAX_DIFF'} : 1.5;

        $PRIMER_MAX_DIFF_TM = exists($param_hash{'PRIMER_RIGHT_INPUT'}) && $clear_form == 0 ? $param_hash{'PRIMER_RIGHT_INPUT'} : "5";

        $PRIMER_OPT_SIZE    = $SNP_PRIMER_OPT_SIZE;
        $PRIMER_MIN_SIZE    = $SNP_PRIMER_MIN_SIZE;
        $PRIMER_MAX_SIZE    = $SNP_PRIMER_MAX_SIZE;
        $PRIMER_OPT_TM      = $SNP_PRIMER_OPT_TM;
        $PRIMER_MIN_TM      = $SNP_PRIMER_MIN_TM;
        $PRIMER_MAX_TM      = $SNP_PRIMER_MAX_TM;
 	$PRIMER_MIN_GC      = $SNP_PRIMER_MIN_GC;
	$PRIMER_MAX_GC      = $SNP_PRIMER_MAX_GC;
   }
    
    if ($primer_type eq '10') { # sequencing primer
        $SNP_PRIMER_MIN_SIZE = exists($param_hash{'SNP_PRIMER_MIN_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MIN_SIZE'} : 18;
        $SNP_PRIMER_OPT_SIZE = exists($param_hash{'SNP_PRIMER_OPT_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_OPT_SIZE'} : 24;
        $SNP_PRIMER_MAX_SIZE = exists($param_hash{'SNP_PRIMER_MAX_SIZE'}) && $clear_form == 0 ? $param_hash{'SNP_PRIMER_MAX_SIZE'} : 30;
    }
 


    print qq{
       	<table width = $HTML_TABLE_WIDTH align =left cellpadding=0 cellspacing=10 class="standard">
 	    <tr><td>
           	<form name="form" action="$PROCESS_INPUT_URL" method="post" enctype="multipart/form-data">
	    	<table border=0 width="100%" cellpadding=0 cellspacing=10>
  	    	    <tr>
    		    	<td colspan=4 align="left"> <font size="+3"> BatchPrimer3</font><br><font size="+2">a high-throughput web tool for picking PCR and sequencing primers</font><p></td>
 		    </tr>

                    <tr>
                        <td colspan=4 align=left>
        		    <strong><a target=_blank href=$BATCHPRIMER3_HOME_URL><b>BatchPrimer3 Home</b></a></strong> |
        		    <strong><a target=_blank href=$DOC_URL>Help</a></strong> |
                            <strong><a target=_blank href="http://jura.wi.mit.edu/primer3/index.php?title=Primer3_Wiki&printable=$CLEAR_FORM">Primer3 Wiki</a></strong> |
                            <strong><a target=_blank href="$HTDOC/disclaimer.html">Copyright Notice and Disclaimer of Primer3</a></strong> |
                            <strong><a target=_blank href="$HTDOC/batchprimer_ack.html">Acknowledgements</a></strong>
                        </td>
                    </tr>
            	</table>

            <p><p>
            <p><p>
            <div class="mainpanel">
            <div class="toppanel">
            <table border=0>
            	<tr>
                    <td><b>Choose primer type:</b></td>
            	    <td>
                          <SELECT NAME="PRIMER_TYPE" onChange="change_type(this.form)">
    };

    
    # print primer type combo box and definition
    for (my $i = 1; $i <= @primer_types; $i++) {
        print "<OPTION value=$i ";
        print "SELECTED" if ($i == $primer_type);
        print ">" . $primer_types[$i-1] . "</OPTION>\n";
    }
    print "</SELECT><p>\n";
    print $primer_type_defs[$primer_type-1] . "\n";

   
    print qq {
            	    </p></td>
                    <td><p>
                    <input type="submit" name="Pick Primers" value="Pick Primers" class="button"> <p>
                    <a href=$INPUT_FORM_URL?PRIMER_TYPE=$primer_type&CLEAR_FORM=yes><b>Reset the entire form</b></a>
                    </td>

                </tr>
     		</table>
            </div>
    };
    
    # user email: the result page cab be sent to a user if a job will take long time.
    print qq {
        <table border=0>
            <tr>
                <td> <b><a href="$DOC_URL#EMAIL">User's e-mail</a> </b>(results will be sent to the user) </td>
                <td> <input type="text"  size=50
                        name="EMAIL" value=$EMAIL> (Optional)</td>
            <tr>
        </table>
    };
    
    # input sequence block
    print qq {
            <h3>Input Sequences:</h3>
            <table border=0>
            	<tr>
                    <td colspan=2 nowrap><a href="$DOC_URL#UPLOAD">Upload sequence file in FASTA format</a>:&nbsp;&nbsp;
                    <input type="file" name="SEQUENCEFILE" size="50" value=$SEQUENCEFILE enctype="multipart/form-data">
  					&nbsp;&nbsp;
    };


    print qq {
                    </td>
            	</tr>

            	<tr>
                	<td colspan=2> <b>OR</b> copy/paste <a href="$DOC_URL#SEQUENCE"> source sequences</a>
                                in FASTA format. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          
        };
        
    if ($primer_type eq '1') {
    	print "<a target=_blank href=$HTDOC/example1.txt>Example sequences</a>";
	} elsif ($primer_type eq '2') {
        	print "<a target=_blank href=$HTDOC/example2.txt>Example sequences</a>";
    } elsif ($primer_type eq '3') {
    	print "<a target=_blank href=$HTDOC/example1.txt>Example sequences</a>";
    } elsif ($primer_type >= 4 && $primer_type <= 9) {
    	print "<a target=_blank href=$HTDOC/example3.txt>Example sequences</a>";
    } elsif ($primer_type == 10) {
    	print "<a target=_blank href=$HTDOC/example4.txt>Example sequences</a>";
    }

    print "&nbsp;&nbsp;&nbsp;<input type=submit name=\"pre-analysis\" value=\"Pre-analysis of input sequences\"  class=\"button\"> ";
    print "&nbsp;&nbsp;<a href=$INPUT_FORM_URL?PRIMER_TYPE=$primer_type&CLEAR_FORM=$CLEAR_SEQ><b>Clear sequence only</b></a><br>";
 

    print qq {        
                    </td>
                </tr>
                <tr>
                    <td colspan=2>
           				<textarea name='SEQUENCE' rows=10 cols=$SOURCE_SEQUENCE_WIDTH>$SEQUENCE</textarea>
                    </td>
                    <p>
                </tr>

                <tr>
                	<td colspan=2>
                		<table border=0>
    };


    if ($primer_type eq '3') {
    	print qq {
            <tr>
                <td> <input type=checkbox name="MUST_XLATE_PICK_HYB_PROBE" checked>
                    Pick hybridization probe (internal oligo) or use the oligo</td>
                <td> <input type="text"  size=50
                    name=PRIMER_INTERNAL_OLIGO_INPUT value=$PRIMER_INTERNAL_OLIGO_INPUT></td>
            </tr>
        };
    }
    elsif ($primer_type ne  '6' && $primer_type ne '8' && $primer_type ne '9'  && $primer_type ne '10') {
    	print qq {
            <tr>
                <td> <input type=checkbox name="MUST_XLATE_PICK_LEFT" value=1 checked>
                    Pick left primer or use the left primer </td>
                <td> <input type="text" size=50
                    name=PRIMER_LEFT_INPUT value=$PRIMER_LEFT_INPUT> </td>
            </tr>
            <tr>
                <td> <input type=checkbox name="MUST_XLATE_PICK_RIGHT" value=1 checked>
                    Pick right primer or use the right primer</td>
                <td> <input type="text"  size=50
                    name=PRIMER_RIGHT_INPUT value=$PRIMER_RIGHT_INPUT></td>
            </tr>
 
            <tr> <td><a name=PRIMER_MISPRIMING_LIBRARY_INPUT href="$DOC_URL#PRIMER_MISPRIMING_LIBRARY">
                    <b>Mispriming/repeat library</b></a>:&nbsp;&nbsp;$SELECT_SEQ_LIBRARY</td>
            </tr>
       };
    }
    
    print qq{
                        </table>
                    </td>
                </tr>

             </table>
            </div>
    };

    #print parameter form for SSR screening
    if ($primer_type eq '2') {
    	print qq {
            <div class="mainpanel">
            <h3> <a name="SSR_SCREENING">SSR Screening</a> </h3>

            <table border=0 cellpadding=0 cellspacing=10>
                <tr>
                	<td colspan=10 align=left><b>Pattern types to be screened:</b></td>
                </tr>
                <tr>
        };
    
        print "<td colspan=2 nowrap>";
        if (exists($param_hash{'DINUCLEOTIDE_SSR'}) || !exists($param_hash{'Pick Primers'}) || $clear_form  == 1) {
            print "<input type=checkbox name=DINUCLEOTIDE_SSR checked> Di-nucleotide</td>";
        } else {
            print "<input type=checkbox name=DINUCLEOTIDE_SSR> Di-nucleotide</td>";
        }
        print "<td colspan=2 nowrap>";
    
        if (exists($param_hash{'TRINUCLEOTIDE_SSR'}) || !exists($param_hash{'Pick Primers'}) || $clear_form  == 1) {
            print "<input type=checkbox name=TRINUCLEOTIDE_SSR checked> Tri-nucleotide</td>";
        } else {
            print "<input type=checkbox name=TRINUCLEOTIDE_SSR> Tri-nucleotide</td>";
        }
        print"<td colspan=2 nowrap>";
    
        if (exists($param_hash{'TETRANUCLEOTIDE_SSR'}) || !exists($param_hash{'Pick Primers'})  || $clear_form  == 1) {
            print "<input type=checkbox name=TETRANUCLEOTIDE_SSR checked> Tetra-nucleotide</td>";
        } else {
            print "<input type=checkbox name=TETRANUCLEOTIDE_SSR> Tetra-nucleotide</td>";
        }
        print "<td colspan=2 nowrap>";
    
        if (exists($param_hash{'PENTANUCLEOTIDE_SSR'}) || !exists($param_hash{'Pick Primers'}) || $clear_form  == 1) {
            print "<input type=checkbox name=PENTANUCLEOTIDE_SSR checked> Penta-nucleotide</td>";
        } else {
            print "<input type=checkbox name=PENTANUCLEOTIDE_SSR> Penta-nucleotide</td>";
        }        
        print "<td colspan=2 nowrap>";
        
        if (exists($param_hash{'HEXANUCLEOTIDE_SSR'}) || !exists($param_hash{'Pick Primers'}) || $clear_form  == 1) {
            print "<input type=checkbox name=HEXANUCLEOTIDE_SSR checked> Hexa-nucleotide</td>";
        } else {
            print "<input type=checkbox name=HEXANUCLEOTIDE_SSR> Hexa-nucleotide</td>";
        }
        print "</tr>";
            
    
        print qq {
                <tr>
                	<td colspan=10 align=left><B>Minimum number of SSR pattern repeats:</B></td>
                </tr>
                <tr>
                	<td>Di-nucleotide</td>
                    <td> <input type="text" size=3 name=DINUCLEOTIDE_SSR_REPEATS value=$DINUCLEOTIDE_SSR_REPEATS_INPUT> </td>
                	<td>Tri-nucleotide</td>
                    <td> <input type="text" size=3 name=TRINUCLEOTIDE_SSR_REPEATS value=$TRINUCLEOTIDE_SSR_REPEATS_INPUT> </td>
                	<td>Tetra-nucleotide</td>
                    <td> <input type="text" size=3 name=TETRANUCLEOTIDE_SSR_REPEATS value=$TETRANUCLEOTIDE_SSR_REPEATS_INPUT> </td>
                	<td>Penta-nucleotide</td>
                    <td> <input type="text" size=3 name=PENTANUCLEOTIDE_SSR_REPEATS value=$PENTANUCLEOTIDE_SSR_REPEATS_INPUT> </td>
                	<td>Hexa-nucleotide</td>
                    <td> <input type="text" size=3 name=HEXANUCLEOTIDE_SSR_REPEATS value=$HEXANUCLEOTIDE_SSR_REPEATS_INPUT> </td>
				</tr>

            </table>

            </div>
        };
    }

    #print parameter form for SNP (allele-specific) primers setting
    if ($primer_type >= 5) {
        # for 5 (SNP flanking and SBE primers), 6 (SBE only) and 7 (allele-specific primers
        # and flanking primers) and 8 (allele-specific primersonly), and 9 (tetra-primer ARMS PCR primers
        if ($primer_type == 9){
            print qq {
                <div class="mainpanel">
                <h3> <a name="GENERAL_SETTINGS">Tetra Primer ARMS PCR Primer Settings</a> </h3>
                <table border=0 cellpadding=0 cellspacing=10>
            };
        }
        elsif ($primer_type == 10){
            print qq {
                <div class="mainpanel">
                <h3> <a name="GENERAL_SETTINGS">Sequencing Primer Settings</a> </h3>
                <table border=0 cellpadding=0 cellspacing=10>
            };
        }
        else {
            print qq {
                <div class="mainpanel">
                <h3> <a name="GENERAL_SETTINGS">SNP or Allele-specific Primer Settings</a> </h3>
                <table border=0 cellpadding=0 cellspacing=10>
            };
        }
        
        
        if ($primer_type == 10){
            print qq {
                <tr>
                    <td><a href="$DOC_URL#SEQUENCING_DIRECTION">Sequencing direction:</a> </td>
                    <td >
                    <SELECT NAME="SEQUENCING_DIRECTION">
                        <OPTION value="toward3" SELECTED>Towards -------> 3'</OPTION>\n";
                        <OPTION value="toward5">         Towards 5'<--------</OPTION>\n";
                    </SELECT>
                    </td>
                    <td><a name=PRIMER_NUM_RETURN_INPUT href="$DOC_URL#PRIMER_NUM_RETURN">Number To Return:</a> </td>
                    <td><input type="text" size=$SMALL_TEXT name=PRIMER_NUM_RETURN value=$PRIMER_NUM_RETURN> </td>
                </tr>
            };
        }
        
        if ($primer_type == 7 || $primer_type == 8){
            print qq {
                <tr>
                    <td colspan=3> <input type=checkbox name="SECOND_MISMATCH" value=1>
                    Additional mismatch is introduced at the specificed position counting from the 3' end</td>
                    <td><input type="text" size=$SMALL_TEXT name=SECOND_MISMATCH_POS value=$SECOND_MISMATCH_POS> </td>
                </tr>
            }
        }
        

        print qq {

            <tr>
                <td><a href="$DOC_URL#PRIMER_SIZE">Primer Size</a> </td>
                <td>Min: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MIN_SIZE value=$SNP_PRIMER_MIN_SIZE> </td>
                <td>Opt: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_OPT_SIZE value=$SNP_PRIMER_OPT_SIZE> </td>
                <td>Max: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MAX_SIZE value=$SNP_PRIMER_MAX_SIZE> </td>
            </tr>

            <tr>
                <td><a href="$DOC_URL#PRIMER_TM">Primer Tm</a> </td>
                <td>Min: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MIN_TM value=$SNP_PRIMER_MIN_TM> </td>
                <td>Opt: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_OPT_TM value=$SNP_PRIMER_OPT_TM> </td>
                <td>Max: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MAX_TM value=$SNP_PRIMER_MAX_TM> </td>
            </tr>
        };
        
        if ($primer_type eq '9') {
            print qq {
                <tr>
                   	<td><a name=PRIMER_MAX_DIFF_TM_INPUT href="$DOC_URL#PRIMER_MAX_DIFF_TM">Max Tm Difference:</a> </td>
                    <td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type="text" size=$SMALL_TEXT name=PRIMER_MAX_DIFF_TM value=$PRIMER_MAX_DIFF_TM> </td>
                    <td></td>
                    <td></td>
                </tr>
            };
        }

        print qq {
            <tr>
                <td><a name=PRIMER_GC_PERCENT_INPUT href="$DOC_URL#PRIMER_GC_PERCENT">Primer GC%</a> </td>
                <td>Min: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MIN_GC value=$SNP_PRIMER_MIN_GC></td>
                <td>Max: <input type="text" size=$SMALL_TEXT name=SNP_PRIMER_MAX_GC value=$SNP_PRIMER_MAX_GC></td>
                <td></td>
            </tr>
        };
        
        if ($primer_type eq '9') {
            print qq {
            <tr>
                <td><a name=PRIMER_PRODUCT_SIZE_INPUT href="$DOC_URL#PRIMER_PRODUCT_SIZE_INPUT">Inner product size</a> </td>
                <td>Min: <input type="text" size=$SMALL_TEXT name=SNP_INNER_PRODUCT_MIN_SIZE value=$SNP_INNER_PRODUCT_MIN_SIZE></td>
                <td>Opt: <input type="text" size=$SMALL_TEXT name=SNP_INNER_PRODUCT_OPT_SIZE value=$SNP_INNER_PRODUCT_OPT_SIZE></td>
                <td>Max: <input type="text" size=$SMALL_TEXT name=SNP_INNER_PRODUCT_MAX_SIZE value=$SNP_INNER_PRODUCT_MAX_SIZE></td>
            </tr>

            <tr>
                <td><a name=PRIMER_PRODUCT_SIZE_INPUT href="$DOC_URL#PRIMER_PRODUCT_SIZE_INPUT">Relative size difference <br>between inner product sizes</a> </td>
                <td>Min: <input type="text" size=$SMALL_TEXT name=SNP_INNER_PRODUCT_MIN_DIFF value=$SNP_INNER_PRODUCT_MIN_DIFF></td>
                <td>Max: <input type="text" size=$SMALL_TEXT name=SNP_INNER_PRODUCT_MAX_DIFF value=$SNP_INNER_PRODUCT_MAX_DIFF></td>
               	<td></td>
            </tr>
            }
        }
        
        print qq {
            <tr>
            	<td><a href="$DOC_URL#PRIMER_NUM_NS_ACCEPTED">Max \#N\'s:</a></td>
              	<td><input type="text" size=$VSMALL_TEXT name=SNP_PRIMER_MAX_N value=$SNP_PRIMER_MAX_N> </td>
               	<td></td>
                <td></td>
            </tr>

            <tr>
                <td><a href="$DOC_URL#PRIMER_SALT_CONC">Salt Concentration:</a></td>
                <td><input type="text" size=$VSMALL_TEXT name=SNP_PRIMER_SALT_CONC value=$SNP_PRIMER_SALT_CONC></td>
                <td><a href="$DOC_URL#PRIMER_DNA_CONC">DNA Concentration:</a></td>
                <td><input type="text" size=$VSMALL_TEXT name=SNP_PRIMER_DNA_CONC value=$SNP_PRIMER_DNA_CONC> </td>
            </tr>
            
            <tr>
            	<td><a href="$DOC_URL#PRIMER_SELF_ANY">Max Self Complementarity:</a></td>
               	<td><input type="text" size=$VSMALL_TEXT name=SNP_MAX_SELF_COMPLEMENTARITY value=$SNP_MAX_SELF_COMPLEMENTARITY></td>
               	<td><a href="$DOC_URL#PRIMER_SELF_END">Max 3\' Self Complementarity:</a></td>
               	<td><input type="text" size=$VSMALL_TEXT name=SNP_MAX_3_SELF_COMPLEMENTARITY value=$SNP_MAX_3_SELF_COMPLEMENTARITY></td>
            </tr>
            
            </table>
            </div>
        };
    }


    #print parameter form for general settings
    if ($primer_type ne '6' && $primer_type ne '8' && $primer_type ne '9' && $primer_type ne '10' ) {
        if ($primer_type eq '3') {
            print qq {
                <div class="mainpanel">
                <h3> <a name="GENERAL_SETTINGS">General Settings for Oligo Primers</a> </h3>
                <table border=0 cellpadding=0 cellspacing=10>
            };
        }
        
        if ($primer_type ne '3') {
            print qq {
                <div class="mainpanel">
                <h3> <a name="GENERAL_SETTINGS">General Settings for Generic Primers</a> </h3>
                <table border=0 cellpadding=0 cellspacing=10>
            };

            print qq {
                 <tr>
                 	<td nowrap><a name=PRIMER_PRODUCT_SIZE_INPUT href="$DOC_URL#PRIMER_PRODUCT_SIZE">Product Size</a>
                 		Min: </td>
                    <td>
                        <input type="text" size=$SMALL_TEXT name=MUST_XLATE_PRODUCT_MIN_SIZE value=$MUST_XLATE_PRODUCT_MIN_SIZE>
                    </td>
                 	<td>Opt:</td>
                    <td>
                     	<input type="text" size=$SMALL_TEXT name=PRIMER_PRODUCT_OPT_SIZE value=$PRIMER_PRODUCT_OPT_SIZE>
                    </td>
                    <td> Max: </td>
                    <td>
                    	<input type="text" size=$SMALL_TEXT name=MUST_XLATE_PRODUCT_MAX_SIZE value=$MUST_XLATE_PRODUCT_MAX_SIZE>
                 	</td>
                 </tr>
            };
        }
        
        print qq {
                 <tr>
                 	<td nowrap><a name=PRIMER_NUM_RETURN_INPUT href="$DOC_URL#PRIMER_NUM_RETURN">Number To Return:</a> </td>
                 	<td><input type="text" size=$SMALL_TEXT name=PRIMER_NUM_RETURN value=$PRIMER_NUM_RETURN> </td>
                 	<td><a name=PRIMER_MAX_END_STABILITY_INPUT href="$DOC_URL#PRIMER_MAX_END_STABILITY">Max 3\' Stability:</a> </td>
                 	<td><input type="text" size=$SMALL_TEXT name=PRIMER_MAX_END_STABILITY value=$PR_DEFAULT_MAX_END_STABILITY></td>
                    <td>&nbsp;</td>
                    <td>&nbsp;</td>
                </tr>
                <tr>
                 	<td nowrap> <a name=PRIMER_MAX_MISPRIMING_INPUT href="$DOC_URL#PRIMER_MAX_MISPRIMING">Max Mispriming:</a>
                 	<td> <input type="text" size=$SMALL_TEXT name=PRIMER_MAX_MISPRIMING value=$PRIMER_MAX_MISPRIMING>
                 	<td> <a name=PRIMER_PAIR_MAX_MISPRIMING_INPUT href="$DOC_URL#PRIMER_PAIR_MAX_MISPRIMING">Pair Max Mispriming:</a>
                 	<td> <input type="text" size=$SMALL_TEXT name=PRIMER_PAIR_MAX_MISPRIMING value=$PRIMER_PAIR_MAX_MISPRIMING>
                    <td>&nbsp;</td>
                    <td>&nbsp;</td>
   				</tr>
        };

   	if ($primer_type ne '3') {
            print qq {

                <tr>
                 	<td><a name=PRIMER_OPT_SIZE_INPUT href="$DOC_URL#PRIMER_SIZE">Primer Size</a> </td>
                 	<td>Min: <input type="text" size=$SMALL_TEXT name=PRIMER_MIN_SIZE value=$PRIMER_MIN_SIZE> </td>
                 	<td>Opt: <input type="text" size=$SMALL_TEXT name=PRIMER_OPT_SIZE value=$PRIMER_OPT_SIZE> </td>
                 	<td>Max: <input type="text" size=$SMALL_TEXT name=PRIMER_MAX_SIZE value=$PRIMER_MAX_SIZE> </td>
                    <td></td>
                    <td></td>
                </tr>

                <tr>
                 	<td><a name=PRIMER_OPT_TM_INPUT href="$DOC_URL#PRIMER_TM">Primer Tm</a> </td>
                 	<td>Min: <input type="text" size=$SMALL_TEXT name=PRIMER_MIN_TM value=$PRIMER_MIN_TM> </td>
                 	<td>Opt: <input type="text" size=$SMALL_TEXT name=PRIMER_OPT_TM value=$PRIMER_OPT_TM> </td>
                 	<td>Max: <input type="text" size=$SMALL_TEXT name=PRIMER_MAX_TM value=$PRIMER_MAX_TM> </td>
                 	<td><a name=PRIMER_MAX_DIFF_TM_INPUT href="$DOC_URL#PRIMER_MAX_DIFF_TM">Max Tm Difference:</a> </td>
                    <td><input type="text" size=$SMALL_TEXT name=PRIMER_MAX_DIFF_TM value=$PRIMER_MAX_DIFF_TM> </td>
                </tr>

                <tr>
                 	<td><a name=PRIMER_PRODUCT_TM_INPUT href="$DOC_URL#PRIMER_PRODUCT_TM">Product Tm</a> </td>
                 	<td>Min: <input type="text" size=$SMALL_TEXT name=PRIMER_PRODUCT_MIN_TM value=''></td>
                 	<td>Opt: <input type="text" size=$SMALL_TEXT name=PRIMER_PRODUCT_OPT_TM value=''></td>
                 	<td>Max: <input type="text" size=$SMALL_TEXT name=PRIMER_PRODUCT_MAX_TM value=''></td>
                    <td></td>
                    <td></td>
                </tr>

                <tr>
                 	<td><a name=PRIMER_GC_PERCENT_INPUT href="$DOC_URL#PRIMER_GC_PERCENT">Primer GC%</a> </td>
                 	<td>Min: <input type="text" size=$SMALL_TEXT name=PRIMER_MIN_GC value=$PRIMER_MIN_GC></td>
                 	<td>Opt: <input type="text" size=$SMALL_TEXT name=PRIMER_OPT_GC_PERCENT value=$PRIMER_OPT_GC_PERCENT></td>
                 	<td>Max: <input type="text" size=$SMALL_TEXT name=PRIMER_MAX_GC value=$PRIMER_MAX_GC></td>
                    <td></td>
                    <td></td>
                </tr>
            </table>

            <table border=0 cellpadding=0 cellspacing=10>

                <tr>
               		<td><a name=PRIMER_SELF_ANY_INPUT href="$DOC_URL#PRIMER_SELF_ANY">Max Self Complementarity:</a></td>
               		<td><input type="text" size=$VSMALL_TEXT name=PRIMER_SELF_ANY value=$PRIMER_SELF_ANY></td>
               		<td><a name=PRIMER_SELF_END_INPUT href="$DOC_URL#PRIMER_SELF_END">Max 3\' Self Complementarity:</a></td>
               		<td><input type="text" size=$VSMALL_TEXT name=PRIMER_SELF_END value=$PRIMER_SELF_END></td>
                </tr>

                <tr>
                	<td><a name=PRIMER_NUM_NS_ACCEPTED_INPUT href="$DOC_URL#PRIMER_NUM_NS_ACCEPTED">Max \#N\'s:</a></td>
                 	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_NUM_NS_ACCEPTED value=$PRIMER_NUM_NS_ACCEPTED> </td>

                 	<td><a name=PRIMER_MAX_POLY_X_INPUT href="$DOC_URL#PRIMER_MAX_POLY_X">Max Poly-X:</a></td>
                 	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_MAX_POLY_X value=$PRIMER_MAX_POLY_X></td>
                </tr>
            };
        }

        print qq {
                <tr>
                  	<td><a name=PRIMER_INSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_INSIDE_PENALTY">Inside Target Penalty:</a></td>
                  	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_INSIDE_PENALTY value=$PRIMER_INSIDE_PENALTY></td>
                  	<td><a name=PRIMER_OUTSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_OUTSIDE_PENALTY">Outside Target Penalty:</a></td>
                  	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_OUTSIDE_PENALTY value=$PRIMER_OUTSIDE_PENALTY></td>
                  	<td><a name=PRIMER_INSIDE_PENALTY_INPUT href="$DOC_URL#PRIMER_INSIDE_PENALTY">
                      Set Inside Target Penalty to allow primers inside a target.
                      </a></td>
                </tr>

                <tr>
               <!--
                 	<td><a name=PRIMER_FIRST_BASE_INDEX_INPUT href="$DOC_URL#PRIMER_FIRST_BASE_INDEX">First Base Index:</a> </td>
                 	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_FIRST_BASE_INDEX value=$PRIMER_FIRST_BASE_INDEX> </td>
                 -->
                 	<td><a name=PRIMER_GC_CLAMP_INPUT href="$DOC_URL#PRIMER_GC_CLAMP">CG Clamp:</a></td>
                 	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_GC_CLAMP value=$PRIMER_GC_CLAMP></td>
                 	<td></td>
                        <td></td>
                        <td></td>
                </tr>

                <tr>
                  	<td><a name=PRIMER_SALT_CONC_INPUT href="$DOC_URL#PRIMER_SALT_CONC">Salt Concentration:</a></td>
                  	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_SALT_CONC value=$PRIMER_SALT_CONC></td>

                  	<td><a name=PRIMER_DNA_CONC_INPUT href="$DOC_URL#PRIMER_DNA_CONC">Annealing Oligo Concentration:</a></td>
                  	<td><input type="text" size=$VSMALL_TEXT name=PRIMER_DNA_CONC value=$PRIMER_DNA_CONC> </td>

                  	<td><a name=PRIMER_DNA_CONC_INPUT href="$DOC_URL#PRIMER_DNA_CONC">
                      (Not the concentration of oligos in the reaction mix but of those
                       annealing to template.)</a></td>
                </tr>
            </table>
        };
        
        if ($primer_type < 4 && $primer_type > 9) { 
            print qq {
                <table border=0 cellpadding=0 cellspacing=10>
                    <tr>
                     	<td><input type=checkbox name=PRIMER_LIBERAL_BASE value=1 checked>
                      		<a name=PRIMER_LIBERAL_BASE_INPUT href="$DOC_URL#PRIMER_LIBERAL_BASE">Liberal Base</a>
                        </td>
                        <td>&nbsp;</td>
      <!--              
                 	<td><input type=checkbox name=MUST_XLATE_PRINT_INPUT value=1>
                  		<a name=SHOW_DEBUGGING_INPUT href="$DOC_URL#SHOW_DEBUGGING">Show Debuging Info</a>
                        </td>
      -->              
                    </tr>
                </table>
            };
        }

        print "    </div>\n";
    }


    #print oligo form
    if ($primer_type eq '3') {
        print qq{
            <div class="mainpanel">

            <h3><a name="INTERNAL_OLIGOS">Hybridization Oligo (Internal Oligo)</a></h3>

            <table border=0>
                <tr><td><a name=internal_oligo_generic_INPUT href="primer3_www_doc.html#internal_oligo_generic">Hyb Oligo Excluded Region:</a></td>
                    <td><input type="text" name=PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION value=$PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION></td>
                </tr>
            </table>

            <h3><a name="Internal+Oligo+Global+Parameters">Hybridization Oligo (Internal Oligo) General Conditions</a></h3>

            <table border=0>
                <tr>
                    <td><a name=PRIMER_INTERNAL_OLIGO_SIZE_INPUT href="$DOC_URL#PRIMER_SIZE">Hyb Oligo Size:</a></td>
                    <td>Min <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_SIZE value=$PRIMER_INTERNAL_OLIGO_MIN_SIZE></td>
                    <td>Opt <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_SIZE value=$PRIMER_INTERNAL_OLIGO_OPT_SIZE></td>
                    <td>Max <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_SIZE value=$PRIMER_INTERNAL_OLIGO_MAX_SIZE></td>
                </tr>

                <tr>
                    <td><a name=PRIMER_OPT_TM_INPUT href="$DOC_URL#PRIMER_TM">Hyb Oligo Tm:</a></td>
                    <td>Min <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_TM value=$PRIMER_INTERNAL_OLIGO_MIN_TM></td>
                    <td>Opt <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_TM value=$PRIMER_INTERNAL_OLIGO_OPT_TM></td>
                    <td>Max <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_TM value=$PRIMER_INTERNAL_OLIGO_MAX_TM></td>
                </tr>

                <tr>
                    <td><a name=PRIMER_INTERNAL_OLIGO_GC_INPUT href="$DOC_URL#PRIMER_GC">Hyb Oligo GC%</a></td>
                    <td>Min: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MIN_GC value=$PRIMER_INTERNAL_OLIGO_MIN_GC></td>
                    <td>Opt: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT value=$PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT></td>
                    <td>Max: <input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_GC value=$PRIMER_INTERNAL_OLIGO_MAX_GC></td>
                </tr>
            </table>

            <table border=0>

                <tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
                    Hyb Oligo Self Complementarity:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SELF_ANY
                               value=$PRIMER_INTERNAL_OLIGO_SELF_ANY></td>
                    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
                        Hyb Oligo Max 3\' Self Complementarity:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SELF_END
                               value=$PRIMER_INTERNAL_OLIGO_SELF_END></td>
                </tr>

                <tr>
                    <td><a name=PRIMER_INTERNAL_OLIGO_NUM_NS_INPUT
                      href="$DOC_URL#internal_oligo_generic">Max #Ns:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_NUM_NS
                            value=$PRIMER_INTERNAL_OLIGO_NUM_NS></td>
                    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Max Poly-X:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_POLY_X value=$PRIMER_INTERNAL_OLIGO_MAX_POLY_X></td>
                </tr>

                <tr>
                    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Mishyb Library:</a></td>
                    <td><select name=PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY>
                      <option selected>NONE</select></td>

                    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Max Mishyb:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_MAX_MISHYB
                             value=$PRIMER_INTERNAL_OLIGO_MAX_MISHYB></td>
                </tr>

                <tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">
                        Hyb Oligo Min Sequence Quality:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT
                         name=PRIMER_INTERNAL_OLIGO_MIN_QUALITY value=$PRIMER_INTERNAL_OLIGO_MIN_QUALITY></td>
                </tr>

                <tr><td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo Salt Concentration:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_SALT_CONC
                               value=$PRIMER_INTERNAL_OLIGO_SALT_CONC></td>
                    <td><a name=internal_oligo_generic_INPUT href="$DOC_URL#internal_oligo_generic">Hyb Oligo DNA Concentration:</a></td>
                    <td><input type="text" size=$VSMALL_TEXT name=PRIMER_INTERNAL_OLIGO_DNA_CONC
                               value=$PRIMER_INTERNAL_OLIGO_DNA_CONC></td>
                </tr>
            </table>

             <h3><a name=OBJECTIVE_PENALTY>Objective Function Penalty Weights for Hyb Oligos (Internal Oligos)</a></h3>

            <table border=0>
        };

        obj_fn_weight_lt_gt('IO_WT_TM',   $PRIMER_IO_WT_TM_LT,   $PRIMER_IO_WT_TM_GT,   'Hyb Oligo Tm');
        obj_fn_weight_lt_gt('IO_WT_SIZE', $PRIMER_IO_WT_SIZE_LT, $PRIMER_IO_WT_SIZE_GT, 'Hyb Oligo Size');
        obj_fn_weight_lt_gt('IO_WT_GC_PERCENT', $PRIMER_IO_WT_GC_PERCENT_LT, $PRIMER_IO_WT_GC_PERCENT_GT, 'Hyb Oligo GC%');

        print "</table>\n\n<table border=0>\n";


        obj_fn_weight('IO_WT_COMPL_ANY', $PRIMER_IO_WT_COMPL_ANY, 'Hyb Oligo Self Complementarity');
        obj_fn_weight('IO_WT_NUM_NS', $PRIMER_IO_WT_NUM_NS, 'Hyb Oligo #N\'s');
        obj_fn_weight('IO_WT_REP_SIM', $PRIMER_IO_WT_REP_SIM, 'Hyb Oligo Mishybing');
        obj_fn_weight('IO_WT_SEQ_QUAL', $PRIMER_IO_WT_SEQ_QUAL, 'Hyb Oligo Sequence Quality');

        print qq{
            </table>
            </div>
        };
    }

	# print form for only single sequence primer design
#    print qq {
#             <div class="mainpanel">
#             <h3> <a name="SINGLE_SEQUENCE">Settings only for single sequence</a> </h3>
#			 <table border=0 width="100%" cellpadding=0 cellspacing=10>
#
#                <tr>
#                    <td nowrap><a name=PRIMER_SEQUENCE_ID_INPUT href="$DOC_URL#PRIMER_SEQUENCE_ID">Sequence Id:</a>  </td>
#                    <td><input type="text" size=20 name=PRIMER_SEQUENCE_ID value=$PRIMER_SEQUENCE_ID> </td>
#                    <td> A string to identify your output. </td>
#                </tr>

#                <tr>
#                    <td nowrap><a name=TARGET_INPUT href="$DOC_URL#TARGET">Targets []:</a> </td>
#                    <td><input type="text" size=20 name=TARGET value=$TARGET></td>
#                    <td colspan=4>E.g. 50,2 requires primers to surround the 2 bases at positions 50 and 51.
#                       Or mark the <a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
#                       with [ and ]: e.g. ...ATCT[CCCC]TCAT.. means
#                       that primers must flank the central CCCC.
#                    </td>
#                </tr>

#                <tr>
#                    <td nowrap><a name=EXCLUDED_REGION_INPUT href="$DOC_URL#EXCLUDED_REGION">Excluded Regions {}:</a> </td>
#                    <td><input type="text" size=20 name=EXCLUDED_REGION value=$EXCLUDED_REGION></td>
#                    <td colspan=4>E.g. 401,7 68,3 forbids selection of primers in the 7 bases starting at 401 and the 3 bases at 68.
#                     Or mark the <a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
#                     with &lt; and &gt;: e.g. ...ATCT&lt;CCCC&gt;TCAT.. forbids
#                     primers in the central CCCC.</td>
#                </tr>

#            	<tr><td nowrap><a name=INCLUDED_REGION_INPUT href="$DOC_URL#INCLUDED_REGION">Included Region <>:</a></td>
#            		<td><input type="text"  size=20 name=INCLUDED_REGION value=$INCLUDED_REGION></td>
#            		<td colspan=4> E.g. 20,400: only pick primers in the 400 base region starting at position 20.
#                		Or use { and } in the
#                		<a href="#PRIMER_SEQUENCE_INPUT">source sequence</a>
#                		to mark the beginning and end of the included region: e.g.
#                		in ATC{TTC...TCT}AT the included region is TTC...TCT. </td>
#                </tr>

#            	<tr><td><a name=PRIMER_START_CODON_POSITION_INPUT href="$DOC_URL#PRIMER_START_CODON_POSITION">Start Codon Position:</a> </td>
#            		<td colspan=4><input type="text"  size=20 name=PRIMER_START_CODON_POSITION value=''> </td>
#                </tr>
#            </table>
#            </div>

# 	};

  #	print_seq_quality_input();

    if ($primer_type ne '3' && $primer_type ne '6' && $primer_type ne '8' && $primer_type ne '9' && $primer_type ne '10') {

        print qq{
            <div class="mainpanel">
            <table border=0  width="100%">
                <tr>
                    <td valign=top width="50%">
                        <h3><a name="PENALTY_WEIGHTS">Penalty Weights for Generic Primers</a></h3>
                        <table border=0>
        };

        obj_fn_weight_lt_gt('WT_TM', $PRIMER_WT_TM_LT, $PRIMER_WT_TM_GT, 'Tm');
        obj_fn_weight_lt_gt('WT_SIZE', $PRIMER_WT_SIZE_LT, $PRIMER_WT_SIZE_GT, 'Size');
        obj_fn_weight_lt_gt('WT_GC_PERCENT', $PRIMER_WT_GC_PERCENT_LT, $PRIMER_WT_GC_PERCENT_GT, 'GC%');

        print "</table><table border=0>";

        obj_fn_weight('WT_COMPL_ANY', $PRIMER_WT_COMPL_ANY, 'Self Complementarity');
        obj_fn_weight('WT_COMPL_END', $PRIMER_WT_COMPL_END, '3\' Self Complementarity');
        obj_fn_weight('WT_NUM_NS', $PRIMER_WT_NUM_NS, '#N\'s');
        obj_fn_weight('WT_REP_SIM', $PRIMER_WT_REP_SIM, 'Mispriming');
        obj_fn_weight('WT_SEQ_QUAL', $PRIMER_WT_SEQ_QUAL, 'Sequence Quality');
        obj_fn_weight('WT_END_QUAL',     $PRIMER_WT_END_QUAL,      'End Sequence Quality');
        obj_fn_weight('WT_POS_PENALTY',  $PRIMER_WT_POS_PENALTY,   'Position Penalty');
        obj_fn_weight('WT_END_STABILITY',$PRIMER_WT_END_STABILITY, 'End Stability');

        print qq{
                        </table>
                    </td>
                    <td valign=top width="50%">
                        <h3><a name="PENALTY_WEIGHTS">Penalty Weights for Generic Primer Pairs</a></h3>
                        <table border=0>
        };

        obj_fn_weight_lt_gt('PAIR_WT_PRODUCT_SIZE', $PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
            $PRIMER_PAIR_WT_PRODUCT_SIZE_GT, 'Product Size');
        obj_fn_weight_lt_gt('PAIR_WT_PRODUCT_TM', $PRIMER_PAIR_WT_PRODUCT_TM_LT,
            $PRIMER_PAIR_WT_PRODUCT_TM_GT, 'Product Tm');

        obj_fn_weight('PAIR_WT_DIFF_TM',   $PRIMER_PAIR_WT_DIFF_TM,   'Tm Difference');
        obj_fn_weight('PAIR_WT_COMPL_ANY', $PRIMER_PAIR_WT_COMPL_ANY, 'Any Complementarity');
        obj_fn_weight('PAIR_WT_COMPL_END', $PRIMER_PAIR_WT_COMPL_END, '3\' Complementarity');
        obj_fn_weight('PAIR_WT_REP_SIM',   $PRIMER_PAIR_WT_REP_SIM,   'Pair Mispriming');
        obj_fn_weight('PAIR_WT_PR_PENALTY',$PRIMER_PAIR_WT_PR_PENALTY, 'Primer Penalty Weight');
        obj_fn_weight('PAIR_WT_IO_PENALTY',$PRIMER_PAIR_WT_IO_PENALTY, 'Hyb Oligo Penalty Weight');

        print qq{
                        </table>
                    <td>
                </tr>
            </table>
            </div>
        };
    }


    # print the end of the html page
    print "    </form>\n";

    print qq {
		</td></tr>
	   	</table>
        </table>
	};

    print $query->end_html;
    print "\n";
}

sub obj_fn_weight_lt_gt {
    my ($pretag, $def_lt, $def_gt, $row_head, $more) = @_;
    # <td><a name="${pretag}_INPUT" href="$DOC_URL#${pretag}">$row_head</a>
    print qq{
  <tr>
    <td><a name="${pretag}_INPUT" href="$DOC_URL#generic_penalty_weights">$row_head</a> &nbsp;
    Lt:<input type=text size=$VSMALL_TEXT name="PRIMER_${pretag}_LT" value=$def_lt></td>
    <td>Gt:<input type=text size=$VSMALL_TEXT name="PRIMER_${pretag}_GT" value=$def_gt></td>
    };
    print "<td> $more\n" if $more;
}

sub obj_fn_weight {
    my ($tag, $def, $row_head, $more) = @_;
    #	<td><a name="${tag}_INPUT" href="$DOC_URL#$tag">$row_head</a>

    print qq{
		<tr>
			<td><a name="${tag}_INPUT" href="$DOC_URL#generic_penalty_weights">$row_head</a></td>
			<td><input type=text size=$VSMALL_TEXT name=PRIMER_$tag value=$def></td>
        </tr>
	};

    print "<td> $more\n" if $more;
}

sub print_seq_quality_input {
    print qq{
        <div class="mainpanel">
		<h3> <a name="SEQUENCE_QUALITY"> Sequence Quality </a></h3>
    	<table border=0>
        	<tr>
                <td>
                	<table border=0>
                        <tr>
        					<td> <a name=PRIMER_SEQUENCE_QUALITY_INPUT href="$DOC_URL#PRIMER_SEQUENCE_QUALITY">Sequence Quality </a></td>
                        </tr>
                        <tr>
                        	<td>
								<textarea rows=4 cols=40 name=PRIMER_SEQUENCE_QUALITY value=$PRIMER_SEQUENCE_QUALITY> </textarea>
                            </td>
                        <tr>
                    </table>
                 </td>
            	<td>
                	<table border=0>
                        <tr>
                            <td><a name=PRIMER_MIN_QUALITY_INPUT href="$DOC_URL#PRIMER_MIN_QUALITY">Min Sequence Quality:</a></td>
                            <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_QUALITY value=$PRIMER_MIN_QUALITY></td>
                        </tr>
                        <tr>
                            <td><a name=PRIMER_MIN_END_QUALITY_INPUT href="$DOC_URL#PRIMER_MIN_END_QUALITY">Min End Sequence Quality:</a></td>
                            <td><input type="text" size=$VSMALL_TEXT name=PRIMER_MIN_END_QUALITY value=$PRIMER_MIN_END_QUALITY></td>
                        </tr>
                        <tr>
                            <td><a name=PRIMER_QUALITY_RANGE_MIN_INPUT href="$DOC_URL#PRIMER_QUALITY_RANGE_MIN">Sequence Quality Range Min:</a></td>
                            <td><input type="text" size=$VSMALL_TEXT name=PRIMER_QUALITY_RANGE_MIN value=$PRIMER_QUALITY_RANGE_MIN></td>
                        </tr>
                        <tr>
                            <td><a name=PRIMER_QUALITY_RANGE_MAX_INPUT href="$DOC_URL#PRIMER_QUALITY_RANGE_MAX">Sequence Quality Range Max:</a></td>
                            <td><input type="text" size=$VSMALL_TEXT name=PRIMER_QUALITY_RANGE_MAX value=$PRIMER_QUALITY_RANGE_MAX></td>
                        </tr>
                    </table>
                </td>
             </tr>
    	</table>
        </div>

    };
}

sub read_user_primer_parameters {
    my %hash = ();
    my $file_name = "$PARAM_DIR/$user_id" . "_" . $primer_type . ".txt";
    if (-e $file_name) {
        open (FILE, "<$file_name")
            or die ("<p>Can't not open the file $file_name!  Please check if your parameter directory is correct. Probably you need to change the mode of the directory: chmod 777</p>");
    
        my $is_sequence = 0;
        my $key = "";
        foreach my $line (<FILE>) {
            chomp($line);
            next if $line =~ /^\s*$/;
            if ($line =~ /=/ && $line !~ /^>/) {
                my @name_value = split("=", $line);
                my $name = $name_value[0];
                my $value = $name_value[1];
                if ($value) {
                     if ($line =~ /^SEQUENCE=/) {
                        $key = $name;
                    } else {
                        $key = "";
                    }
                    $hash{$name} = $value;
                }
            }
            elsif ($line !~ /^\s*$/ && $key eq 'SEQUENCE') {
                $hash{$key} = $hash{$key} . "\n" . $line; 
            }
        }
        close (FILE);
    }
    
    return %hash;
}

