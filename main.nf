#!/usr/bin/env nextflow
/*
 * Authors:
 * - Daniel Cook <danielecook@gmail.com>
 *
 */


/*
    Globals
*/

// Define contigs here!
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contigs = Channel.from(CONTIG_LIST)

/*
    Params
*/

date = new Date().format( 'yyyyMMdd' )
params.out = "WI-${date}"
params.debug = false
params.annotation_reference = "WS263"
params.cores = 6
params.tmpdir = "tmp/"
params.email = ""
params.reference = "(required)"
params.manta_path = null
params.tiddit_discord = null
params.snpeff_path="${workflow.workDir}/snpeff"
params.call_sv = false


// Compressed Reference File
File reference = new File("${params.reference}")
if (params.reference != "(required)") {
   reference_handle = reference.getAbsolutePath();
   reference_handle_uncompressed = reference_handle.replace(".gz", "")
} else {
   reference_handle = "(required)"
}

// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.fqs = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bamdir = "${params.out}/bam"
    File fq_file = new File(params.fqs);
    params.fq_file_prefix = "${workflow.projectDir}/test_data"

    // DEBUG Filter thresholds
    min_depth=0
    qual=10
    mq=10
    dv_dp=0.0


} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fq_file_prefix = null;
    params.fqs = "sample_sheet.tsv"

    min_depth=10
    qual=30
    mq=40
    dv_dp=0.5

}

File fq_file = new File(params.fqs);

/*
    ==
    UX
    ==
*/

param_summary = '''


     ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄                         ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄
    ▐░▌       ▐░▌▐░░░░░░░░░░░▌                        ▐░░▌      ▐░▌▐░░░░░░░░░░░▌
    ▐░▌       ▐░▌ ▀▀▀▀█░█▀▀▀▀                        ▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌       ▐░▌     ▐░▌                             ▐░▌▐░▌    ▐░▌▐░▌
    ▐░▌   ▄   ▐░▌     ▐░▌           ▄▄▄▄▄▄▄▄▄▄▄      ▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄
    ▐░▌  ▐░▌  ▐░▌     ▐░▌          ▐░░░░░░░░░░░▌      ▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌
    ▐░▌ ▐░▌░▌ ▐░▌     ▐░▌           ▀▀▀▀▀▀▀▀▀▀▀      ▐░▌   ▐░▌ ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌▐░▌ ▐░▌▐░▌     ▐░▌                             ▐░▌    ▐░▌▐░▌▐░▌
    ▐░▌░▌   ▐░▐░▌ ▄▄▄▄█░█▄▄▄▄                        ▐░▌     ▐░▐░▌▐░▌
    ▐░░▌     ▐░░▌▐░░░░░░░░░░░▌                        ▐░▌      ▐░░▌▐░▌
     ▀▀       ▀▀  ▀▀▀▀▀▀▀▀▀▀▀                         ▀        ▀▀  ▀


''' + """

    parameters              description                    Set/Default
    ==========              ===========                    =======

    --debug                 Set to 'true' to test          ${params.debug}
    --cores                 Regular job cores              ${params.cores}
    --out                   Directory to output results    ${params.out}
    --fqs                   fastq file (see help)          ${params.fqs}
    --fq_file_prefix        fastq prefix                   ${params.fq_file_prefix}
    --reference             Reference Genome (w/ .gz)      ${params.reference}
    --annotation_reference  SnpEff annotation              ${params.annotation_reference}
    --bamdir                Location for bams              ${params.bamdir}
    --tmpdir                A temporary directory          ${params.tmpdir}
    --email                 Email to be sent results       ${params.email}

    HELP: http://andersenlab.org/dry-guide/pipeline-wi/

"""

println param_summary

if (params.reference == "(required)" || params.fqs == "(required)") {

    println """
    The Set/Default column shows what the value is currently set to
    or would be set to if it is not specified (it's default).
    """
    System.exit(1)
}

if (!reference.exists()) {
    println """

    Error: Reference does not exist

    """
    System.exit(1)
}

if (!fq_file.exists()) {
    println """

    Error: fastq sheet does not exist

    """
    System.exit(1)
}


// Read sample sheet
strainFile = new File(params.fqs)

if (params.fq_file_prefix != "") {
    fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
                 .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${params.fq_file_prefix}/${fq1}"), file("${params.fq_file_prefix}/${fq2}"), seq_folder] }
} else {
    fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
                 .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${fq1}"), file("${fq2}"), seq_folder] }
}


fqs.into {
    fqs_kmer
    fqs_align
}

/*
    =============
    Kmer counting
    =============
*/
process kmer_counting {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2, seq_folder from fqs_kmer
    output:
        file("${ID}.kmer.tsv") into kmer_set

    """
        # fqs will have same number of lines
        export OFS="\t"
        fq_wc="`zcat ${fq1} | awk 'NR % 4 == 0' | wc -l`"
        
        zcat ${fq1} ${fq2} | \\
        fastq-kmers -k 6 | \\
        awk -v OFS="\t" -v ID=${ID} -v SM=${SM} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${ID}.kmer.tsv
    """
}


process merge_kmer {

    publishDir "${params.out}/phenotype", mode: "copy"

    input:
        file("kmer*.tsv") from kmer_set.collect()
    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tSM\tID\twc") *.tsv > kmers.tsv
    """

}


/*
    ===============
    Fastq alignment
    ===============

    The output looks strange below,
    but its designed to group like samples together - so leave it!

*/

process perform_alignment {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2, seq_folder from fqs_align
    output:
        set val(SM), file("${ID}.bam"), file("${ID}.bam.bai") into fq_bam_set


    """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference_handle} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
        sambamba index --nthreads=${task.cpus} ${ID}.bam

        if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
            exit 1;
        fi

    """
}

/*
    ===================================
    Merge - Generate isotype-level BAMs
    ===================================
*/

process merge_bam {

    cpus params.cores

    tag { SM }

    input:
        set SM, bam, index from fq_bam_set.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into SM_bam_set
        file("${SM}.picard.sam.markduplicates") into duplicates_set

    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${task.cpus} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${task.cpus} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.picard.sam.markduplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${task.cpus} ${SM}.bam
    """
}

SM_bam_set.into {
                  bam_publish;
                  bam_idxstats;
                  bam_stats;
                  bam_coverage;
                  bam_snp_individual;
                  bam_snp_union;
                  bam_telseq;
                  bam_isotype_stats;
                  bam_manta;
                  bam_tiddit;
                  bam_delly;
                  bam_delly_recall;
}

process bam_isotype_stats {

    cpus params.cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_isotype_stats

    output:
         file("${SM}.samtools.txt") into SM_samtools_stats_set
         file("${SM}.bamtools.txt") into SM_bamtools_stats_set
         file("${SM}_fastqc.zip") into SM_fastqc_stats_set
         file("${SM}.picard.*") into SM_picard_stats_set

    """
        samtools stats --threads=${task.cpus} ${SM}.bam > ${SM}.samtools.txt
        bamtools -in ${SM}.bam > ${SM}.bamtools.txt
        fastqc --threads ${task.cpus} ${SM}.bam
        picard CollectAlignmentSummaryMetrics R=${reference_handle} I=${SM}.bam O=${SM}.picard.alignment_metrics.txt
        picard CollectInsertSizeMetrics I=${SM}.bam O=${SM}.picard.insert_metrics.txt H=${SM}.picard.insert_histogram.txt
    """

}

process bam_publish {

    publishDir "${params.bamdir}/WI/isotype", mode: 'copy', pattern: '*.bam*'

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_publish
    output:
        set file("${SM}.bam"), file("${SM}.bam.bai")

    """
        echo "${SM} saved to publish folder you rockstar."
    """
}

process SM_idx_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_idxstats
    output:
        file("${SM}.bam_idxstats") into bam_idxstats_set
        file("${SM}.bam_idxstats") into bam_idxstats_multiqc

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > ${SM}.bam_idxstats
    """
}

process SM_combine_idx_stats {

    publishDir "${params.out}/alignment", mode: "copy"

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("isotype_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > isotype_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> isotype_bam_idxstats.tsv
    """
}

/*
    =================
    Isotype BAM stats
    =================
*/

process isotype_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        samtools stats ${SM}.bam | \\
        grep ^SN | \\
        cut -f 2- | \\
        awk '{ print "${SM}\t" \$0 }' | \\
        sed 's/://g' > bam_stat
    """
}

process combine_isotype_bam_stats {

    publishDir "${params.out}/alignment", mode: "copy"

    input:
        val stat_files from SM_bam_stat_files.collect()

    output:
        file("isotype_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > isotype_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_stats.tsv
    """
}

/*
    ============
    Coverage BAM
    ============
*/
process coverage_SM {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_coverage

    output:
        val SM into isotype_coverage_sample
        file("${SM}.coverage.tsv") into isotype_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process coverage_SM_merge {

    publishDir "${params.out}/alignment", mode: 'copy'

    input:
        val sm_set from isotype_coverage.toSortedList()

    output:
        file("isotype_coverage.full.tsv") into mt_content
        file("isotype_coverage.tsv") into isotype_coverage_merged

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > isotype_coverage.full.tsv
        cat ${sm_set.join(" ")} >> isotype_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat isotype_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > isotype_coverage.tsv
    """
}

/*
    ==========
    MT content
    ==========
*/

process output_mt_content {

    executor 'local'

    publishDir "${params.out}/phenotype", mode: "copy"

    input:
        file("isotype_coverage.full.tsv") from mt_content

    output:
        file("MT_content.tsv")

    """
        cat <(echo -e 'isotype\\tmt_content') <(cat isotype_coverage.full.tsv | awk '/mt_nuclear_ratio/' | cut -f 1,6) > MT_content.tsv
    """
}

/*
    ======
    telseq
    ======
*/

process call_telseq {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_telseq
    output:
        file("telseq_out.txt") into telseq_results

    """
        telseq -z TTAGGC -H ${SM}.bam > telseq_out.txt
    """
}

process combine_telseq {

    executor 'local'

    publishDir "${params.out}/phenotype", mode: 'copy'

    input:
        file("ind_telseq?.txt") from telseq_results.toSortedList()

    output:
        file("telseq.tsv")

    '''
        telseq -h > telseq.tsv
        cat ind_telseq*.txt | egrep -v '\\[|BAMs' >> telseq.tsv
    '''
}

/*
    ========================
    Call Variants - BCFTools
    ========================    
*/

process call_variants {

    tag { SM }

    cpus params.cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_snp_individual

    output:
        file("${SM}.vcf.gz") into isotype_vcf

    """

    # Subsample high-depth bams
    coverage=`goleft covstats ${SM}.bam | awk 'NR > 1 { printf "%5.0f", \$1 }'`

    if [ \${coverage} -gt 100 ];
    then

        # Add a trap to remove temp files
        function finish {
            rm -f "${SM}.subsample.bam"
            rm -f "${SM}.subsample.bam.bai"
        }
        trap finish EXIT

        echo "Coverage is above 100x; Subsampling to 100x"
        # Calculate fraction of reads to keep
        frac_keep=`echo "100.0 / \${coverage}" | bc -l | awk '{printf "%0.2f", \$0 }'`
        SM_use="${SM}.subsample.bam"
        sambamba view --nthreads=${task.cpus} --show-progress --format=bam --with-header --subsample=\${frac_keep} ${SM}.bam > \${SM_use}
        sambamba index --nthreads ${task.cpus} \${SM_use}
    else
        echo "Coverage is below 100x; No subsampling"
        SM_use="${SM}.bam"
    fi;


    function process_variants {
        bcftools mpileup --redo-BAQ \\
                         --redo-BAQ \\
                         -r \$1 \\
                         --gvcf 1 \\
                         --annotate DP,AD,ADF,ADR,INFO/AD,SP \\
                         --fasta-ref ${reference_handle} \${SM_use} | \\
        bcftools call --multiallelic-caller \\
                      --gvcf 3 \\
                      --multiallelic-caller -O v - | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --mode + --soft-filter quality --include "(QUAL >= ${qual}) || (FORMAT/GT == '0/0') || (TYPE == 'REF')" |  \\
        bcftools filter -O u --mode + --soft-filter min_depth --include "(FORMAT/DP > ${min_depth}) || (TYPE == 'REF')" | \\
        bcftools filter -O u --mode + --soft-filter mapping_quality --include "(INFO/MQ > ${mq}) || (TYPE == 'REF')" | \\
        bcftools filter -O v --mode + --soft-filter dv_dp --include "((FORMAT/AD[*:1])/(FORMAT/DP) >= ${dv_dp}) || (FORMAT/GT == '0/0') || (TYPE == 'REF')" | \\
        awk -v OFS="\\t" '\$0 ~ "^#" { print } \$0 ~ ":AB" { gsub("PASS","", \$7); if (\$7 == "") { \$7 = "het"; } else { \$7 = \$7 ";het"; } } \$0 !~ "^#" { print }' | \\
        awk -v OFS="\\t" '\$0 ~ "^#CHROM" { print "##FILTER=<ID=het,Description=\\"heterozygous_call_after_het_polarization\\">"; print; } \$0 ~ "^#" && \$0 !~ "^#CHROM" { print } \$0 !~ "^#" { print }' | \\
        vk geno transfer-filter - | \\
        bcftools norm -O z --check-ref s --fasta-ref ${reference_handle} > ${SM}.\$1.vcf.gz
    }

    export SM_use;
    export -f process_variants

    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    parallel -j ${task.cpus} --verbose process_variants {} ::: \${contigs}
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".vcf.gz" }'`
    
    # Concatenate and filter
    bcftools concat --threads ${task.cpus-1} \${order} -O z > ${SM}.vcf.gz
    bcftools index --threads ${task.cpus} ${SM}.vcf.gz
    rm \${order}

    """
}


process generate_vcf_list {

    executor 'local'

    cpus 1 

    input:
       val vcf_set from isotype_vcf.toSortedList()

    output:
       file("union_vcfs.txt") into union_vcfs

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
    """
}

union_vcfs_in = union_vcfs.spread(contigs)

process merge_union_vcf_chromosome {

    cpus params.cores

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in

    output:
        set val(chrom), file("${chrom}.merged.vcf.gz"), file("${chrom}.merged.vcf.gz.csi") into raw_vcf

    """
        bcftools merge --gvcf ${reference_handle} \\
                       --regions ${chrom} \\
                       -O z \\
                       -m both \\
                       --file-list ${union_vcfs} > ${chrom}.merged.vcf.gz
        echo "Merging of ${chrom} completed";
        bcftools index ${chrom}.merged.vcf.gz
    """
}


// Generates the initial soft-vcf; but it still
process generate_soft_vcf {

    cpus params.cores

    tag { chrom }

    input:
        set val(chrom), file("${chrom}.merged.vcf.gz"), file("${chrom}.merged.vcf.gz.csi") from raw_vcf

    output:
        set val(chrom), file("${chrom}.soft-filter.vcf.gz"), file("${chrom}.soft-filter.vcf.gz.csi") into soft_filtered_vcf

    """
        bcftools view --threads=${task.cpus-1} ${chrom}.merged.vcf.gz | \\
        bcftools filter -O u --mode=+x --soft-filter="high_missing" --include 'F_MISSING  <= ${params.missing}' - | \\
        bcftools filter -O u --mode=+x --soft-filter="high_heterozygosity" --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' - | \\
        bcftools view -O v --min-af 0.0000000000001 --max-af 0.999999999999 | \\
        vcffixup - | \\
        bcftools norm -O u -m +both --rm-dup both --fasta-ref ${reference_handle} | \\
        bcftools view --threads=${task.cpus-1} -O z - > ${chrom}.soft-filter.vcf.gz
        bcftools index --threads=${task.cpus} -f ${chrom}.soft-filter.vcf.gz
    """
}


/*
    Fetch some necessary datasets
*/

process fetch_ce_gff {

    executor 'local'

    output:
        file("ce.gff3.gz") into ce_gff3
    
    """
        # Download the annotation file
        wget -O ce.gff3.gz ftp://ftp.ensembl.org/pub/release-92/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.92.gff3.gz
    """
}


fix_snpeff_script = file("fix_snpeff_names.py")

process fetch_gene_names {

    executor 'local'

    output:
        file("gene.pkl") into gene_pkl

    """
    fix_snpeff_names.py
    """

}

gene_pkl.into {
                    gene_pkl_snpindel;
                    gene_pkl_manta;
                    gene_pkl_delly;
                    gene_pkl_tiddit;
                    gene_pkl_cnv;
                    gene_pkl_svdb;
              }


process annotate_vcf {

    cpus params.cores

    tag { chrom }

    cache 'deep'

    input:
        set val(chrom), file("${chrom}.soft-filter.vcf.gz"), file("${chrom}.soft-filter.vcf.gz.csi") from soft_filtered_vcf
        file 'snpeff.config' from file("${baseDir}/snpeff_data/snpEff.config")
        file 'snpeff_data' from file("${baseDir}/snpeff_data")
        file 'ce.gff3.gz' from ce_gff3
        file 'gene.pkl' from gene_pkl_snpindel

    output:
        file("${chrom}.soft-annotated.vcf.gz") into soft_annotated_vcf
        file("snpeff_out.csv") into snpeff_multiqc


    """
        # bcftools csq
        bcftools view --threads=${task.cpus-1} -O v ${chrom}.soft-filter.vcf.gz | \\
        bcftools csq -O v --fasta-ref ${reference_handle} \\
                     --gff-annot ce.gff3.gz \\
                     --phase a | \\
        snpEff eff -csvStats snpeff_out.csv \\
                   -no-downstream \\
                   -no-intergenic \\
                   -no-upstream \\
                   -nodownload \\
        -dataDir . \\
        -config snpeff.config \\
        ${params.annotation_reference} | \\
        bcftools view -O v | \\
        fix_snpeff_names.py - | \\
        bcftools view --threads=${task.cpus-1} -O z > ${chrom}.soft-annotated.vcf.gz
        bcftools index --threads=${task.cpus} ${chrom}.soft-annotated.vcf.gz
    """

}


// Generate a list of ordered files.
contig_raw_vcf = CONTIG_LIST*.concat(".soft-annotated.vcf.gz")

process concatenate_union_vcf {

    cpus params.cores

    tag { chrom }

    input:
        val merge_vcf from soft_annotated_vcf.collect()

    output:
        set file("full.soft-filter.vcf.gz"), file("full.soft-filter.vcf.gz.csi") into soft_filtered_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat --threads ${task.cpus-1} -O z ${contig_raw_vcf.join(" ")} > full.soft-filter.vcf.gz
        bcftools index  --threads ${task.cpus} full.soft-filter.vcf.gz
    """
}

/*
    Download annotation tracks
*/

process download_annotation_files {

    executor 'local'

    errorStrategy 'retry'
    maxRetries 5

    output:
        set val("phastcons"), file("elegans.phastcons.wib") into phastcons
        set val("phylop"), file("elegans.phylop.wib") into phylop
        set val("repeatmasker"), file("elegans_repeatmasker.bb") into repeatmasker

    """
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phastcons.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phylop.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans_repeatmasker.bb
    """
}

phastcons.mix(phylop).set { wig }

process wig_to_bed {

    tag { track_name }

    publishDir "${params.out}/tracks", mode: 'copy'

    input:
        set val(track_name), file("track.wib") from wig
    output:
        set file("${track_name}.bed.gz"), file("${track_name}.bed.gz.tbi") into bed_tracks
        
    """
        bigWigToBedGraph track.wib ${track_name}.bed
        bgzip ${track_name}.bed
        tabix ${track_name}.bed.gz
    """

}

process annovar_and_output_soft_filter_vcf {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
        set file("WI.${date}.soft-effect.vcf.gz"), file("WI.${date}.soft-effect.vcf.gz.csi") from soft_filtered_concatenated
        file(track) from bed_tracks.collect()
        file('vcf_anno.conf') from Channel.fromPath("data/vcfanno.conf")

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into soft_filter_vcf_annotated
        set val("soft"), file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into soft_sample_summary
        file("WI.${date}.soft-filter.stats.txt") into soft_filter_stats
        file("WI.${date}.soft-filter.vcf.gz.tbi")

    """
        vcfanno -p ${task.cpus} vcf_anno.conf WI.${date}.soft-effect.vcf.gz | \\
        bcftools view --threads ${task.cpus-1} -O z > WI.${date}.soft-filter.vcf.gz
        bcftools index --threads ${task.cpus} WI.${date}.soft-filter.vcf.gz
        tabix WI.${date}.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt
    """

}


soft_filter_vcf_annotated.into {
                                 soft_filtered_vcf_to_hard;
                                 soft_filtered_vcf_gtcheck;
                                 soft_filter_vcf_strain;
                                 soft_filter_vcf_isotype_list;
                                 soft_filter_vcf_mod_tracks;
                                 soft_filter_vcf_tsv;
                               }


process generate_strain_list {

    executor 'local'

    input:
        set file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi") from soft_filter_vcf_isotype_list

    output:
        file('isotype_list.tsv') into isotype_list

    """
        bcftools query -l WI.${date}.vcf.gz > isotype_list.tsv
    """

}

isotype_list.into {
    sample_files_list;
}


process generate_hard_vcf {

    cpus params.cores

    cache 'deep'

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from soft_filtered_vcf_to_hard

    output:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf
        set val("hard"), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf_summary
        set val("hard"), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_sample_summary
        file("WI.${date}.hard-filter.vcf.gz.tbi")
        file("WI.${date}.hard-filter.stats.txt") into hard_filter_stats


    """
        # Generate hard-filtered VCF
        function generate_hard_filter {
            bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} WI.${date}.soft-filter.vcf.gz | \\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT != "PASS"' - | \\
            bcftools filter -O u --include 'F_MISSING  <= ${params.missing}' - | \\
            bcftools filter -O u --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' - | \\
            bcftools view -O v --min-af 0.0000000000001 --max-af 0.999999999999 | \\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
        }

        export -f generate_hard_filter

        parallel --verbose generate_hard_filter {} ::: I II III IV V X MtDNA

        bcftools concat -O z I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz > WI.${date}.hard-filter.vcf.gz
        bcftools index --threads ${task.cpus} WI.${date}.hard-filter.vcf.gz
        tabix WI.${date}.hard-filter.vcf.gz
        bcftools stats --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt

        # Remove extra files
        rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
    """
}

hard_vcf.into { 
                hard_vcf_to_impute;
                tajima_bed;
                vcf_phylo;
                hard_vcf_variant_accumulation;
                biallelic_snp_vcf
            }


process generate_biallelic_snp_vcf {

    cpus params.cores

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from biallelic_snp_vcf

    output:
        set file("WI.${date}.hard-filter-biallelic-snp.vcf.gz"), file("WI.${date}.hard-filter-biallelic-snp.vcf.gz.csi") into haplotype_vcf
        file("WI.${date}.hard-filter-biallelic-snp.stats.txt") into biallelic_snp_stats

    """
        bcftools view --threads ${task.cpus-1} \\
                      -m 2 \\
                      -M 2 \\
                      --trim-alt-alleles \\
                      --min-af 0.0000000000001 \\
                      --max-af 0.999999999999 \\
                      --include '%TYPE == "SNP"' \\
                      -O z \\
                      WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter-biallelic-snp.vcf.gz

        bcftools index --threads ${task.cpus} WI.${date}.hard-filter-biallelic-snp.vcf.gz
        bcftools stats --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter-biallelic-snp.stats.txt 
    """

}


process calculate_gtcheck {

    publishDir "${params.out}/concordance", mode: 'copy'

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from soft_filtered_vcf_gtcheck

    output:
        file("gtcheck.tsv") into gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 WI.${date}.soft-filter.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """
}


/*
    =================
    Calculate Summary
    =================
*/
process calculate_hard_vcf_summary {

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set val('hard'), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from hard_vcf_summary

    output:
        file("WI.${date}.hard-filter.genotypes.tsv")
        file("WI.${date}.hard-filter.genotypes.frequency.tsv")

    """
        # Calculate singleton freq
        vk calc genotypes WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.genotypes.tsv
        vk calc genotypes --frequency WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.genotypes.frequency.tsv

        # Calculate average discordance; Determine most diverged strains
        awk '\$0 ~ "^CN" { print 1-(\$2/\$3) "\t" \$5 "\n" 1-(\$2/\$3) "\t" \$6 }' | \
        sort -k 2 | \
        datamash mean 1 --group 2 | \
        sort -k2,2n > WI.${date}.hard-filter.avg_concordance.tsv
    """
}


/*
    Variant summary
*/

sample_summary = soft_sample_summary.concat( hard_sample_summary )

process sample_variant_summary {

    tag { summary_vcf }

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set val(summary_vcf), file("out.vcf.gz"), file("out.vcf.gz.csi") from sample_summary

    output:
        set val(summary_vcf), file("${summary_vcf}.variant_summary.json") into sample_summary_out

    """
    bcftools view out.vcf.gz | sample_summary_vcf.py - > ${summary_vcf}.variant_summary.json
    """
}


sample_summary_out.into {
    sample_summary_to_tsv;
    sample_summary_to_split
}

process parse_sample_summary {

    tag { summary_vcf }

    publishDir "${params.out}/variation/sample_summary", mode: 'copy'

    input:
        set val(summary_vcf), file("${summary_vcf}.variant_summary_in.json") from sample_summary_to_tsv

    output:
        file("${summary_vcf}.effect_summary.json")
        file("${summary_vcf}.effect_summary.tsv")
        file("${summary_vcf}.impact_summary.json")
        file("${summary_vcf}.impact_summary.tsv")
        file("${summary_vcf}.biotype_summary.json")
        file("${summary_vcf}.biotype_summary.tsv")
        file("${summary_vcf}.high_impact_variants_summary.json")
        file("${summary_vcf}.high_impact_variants_summary.tsv")
        file("${summary_vcf}.gt_count_summary.json")
        file("${summary_vcf}.gt_count_summary.tsv")
        file("${summary_vcf}.isotype_summary.json")
        file("${summary_vcf}.isotype_summary.tsv")

    """
        # Parse variant summary json
        cat ${summary_vcf}.variant_summary_in.json  | jq 'keys[] as \$parent | {'sample': \$parent} + .[\$parent].ANN.effect' | jq --slurp '.' > ${summary_vcf}.effect_summary.json
        cat ${summary_vcf}.variant_summary_in.json  | jq 'keys[] as \$parent | {'sample': \$parent} + .[\$parent].ANN.impact' | jq --slurp '.' > ${summary_vcf}.impact_summary.json
        cat ${summary_vcf}.variant_summary_in.json  | jq 'keys[] as \$parent | {'sample': \$parent} + .[\$parent].ANN.transcript_biotype' | jq --slurp '.' > ${summary_vcf}.biotype_summary.json
        cat ${summary_vcf}.variant_summary_in.json  | jq 'keys[] as \$parent | {'sample': \$parent} + .[\$parent].ANN.HIGH_impact_genes[]?' | jq --slurp '.' > ${summary_vcf}.high_impact_variants_summary.json
        cat ${summary_vcf}.variant_summary_in.json  | jq 'keys[] as \$parent | 
                                                (.[\$parent].gt_count | keys[]) as \$subcat | 
                                                (.[\$parent].gt_count[\$subcat] | keys[]) as \$n_sample | 
                                                (.[\$parent].gt_count[\$subcat][\$n_sample]) as \$n_count | 
                                                {'sample': \$parent, 'gt': \$subcat, 'n_sample': (\$n_sample | tonumber), 'n_count': \$n_count} ' | jq --slurp '.' > ${summary_vcf}.gt_count_summary.json


        cat ${summary_vcf}.variant_summary_in.json | \
            jq 'keys[] as \$parent |
                (.[\$parent].gt_count.homozygous_alt["1"] + .[\$parent].gt_count.homozygous_ref["1"]) as \$singleton | 
                ([.[\$parent].gt_count.homozygous_alt[]?]  | add) as \$alt_calls |
                ([.[\$parent].gt_count.homozygous_ref[]?] | add) as \$ref_calls |
                ([.[\$parent].gt_count.heterozygous[]?]  | add) as \$het_calls |
                ([.[\$parent].gt_count.missing[]?]  | add) as \$missing_calls |
                (\$alt_calls + \$ref_calls + \$het_calls) as \$n_variants |
                .[\$parent].ANN.impact.HIGH as \$high_impact |
                {"isotype": \$parent,
                 "singletons": \$singleton,
                 "alt_calls": \$alt_calls,
                 "ref_calls": \$ref_calls,
                 "het_calls": \$het_calls,
                 "missing_calls": \$missing_calls,
                 "n_calls": \$n_variants,
                 "high_impact_variants": \$high_impact}' | jq --slurp '.' > ${summary_vcf}.isotype_summary.json

        Rscript `which process_variant_summary.R`
    """

}


/*
    ==============
    Phylo analysis
    ==============
*/
process phylo_analysis {

    publishDir "${params.out}/popgen/trees", mode: "copy"

    tag { contig }

    input:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi"), val(contig) from vcf_phylo.spread(["I", "II", "III", "IV", "V", "X", "MtDNA", "genome"])

    output:
        set val(contig), file("${contig}.tree") into trees

    """
        # genome
        if [ "${contig}" == "genome" ]
        then
            vk phylo tree nj WI.${date}.hard-filter.vcf.gz > genome.tree
            if [[ ! genome.tree ]]; then
                exit 1;
            fi
        else
            vk phylo tree nj WI.${date}.hard-filter.vcf.gz ${contig} > ${contig}.tree
            if [[ ! ${contig}.tree ]]; then
                exit 1;
            fi
        fi
    """
}


process plot_trees {

    publishDir "${params.out}/popgen/trees", mode: "copy"

    tag { contig }

    input:
        set val(contig), file("${contig}.tree") from trees

    output:
        file("${contig}.svg")
        file("${contig}.png")


    """
        Rscript --vanilla `which process_trees.R` ${contig}
    """

}


process tajima_bed {

    publishDir "${params.out}/popgen", mode: 'copy'

    input:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from tajima_bed
    output:
        set file("WI.${date}.tajima.bed.gz"), file("WI.${date}.tajima.bed.gz.tbi") into plot_tajima

    """
        vk tajima --no-header 100000 10000 WI.${date}.hard-filter.vcf.gz | bgzip > WI.${date}.tajima.bed.gz
        tabix WI.${date}.tajima.bed.gz
    """

}


process plot_tajima {

    executor 'local'

    publishDir "${params.out}/popgen", mode: 'copy'

    input:
        set file("tajima.bed.gz"), file("tajima.bed.gz.tbi") from plot_tajima

    output:
        file("tajima_d.png")
        file("tajima_d.thumb.png")

    """
        Rscript --vanilla `which plot_tajima.R`

        # Generate thumbnail
        convert tajima_d.png -density 300 -resize 200 tajima_d.thumb.png
    """

}

/*
    ====================
    variant_accumulation
    ====================
*/

process calc_variant_accumulation {

    publishDir "${params.out}/popgen", mode: 'copy'

    input: 
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from hard_vcf_variant_accumulation

    output:
        file("variant_accumulation.pdf")
        file("variant_accumulation.tsv")
        file("impute_gts.tsv.gz")

    """
    bcftools query -f "[%GT\\t]\n" WI.${date}.hard-filter.vcf.gz  | \\
    awk '{ gsub(":GT", "", \$0); gsub("(# )?[[0-9]+]","",\$0); print \$0 }' | \\
    sed -r 's/([0-9\\.]+)\\/([0-9\\.]+)/\\1/g' | \\
    sed 's/\\./NA/g' | \\
    gzip > impute_gts.tsv.gz

    Rscript --vanilla `which variant_accumulation.R`
    """
}


process imputation {

    cpus params.cores

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from hard_vcf_to_impute
    output:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") into impute_vcf
        file("WI.${date}.impute.stats.txt") into impute_stats
        file("WI.${date}.impute.stats.txt") into filtered_stats
        file("WI.${date}.impute.vcf.gz")
        file("WI.${date}.impute.vcf.gz.tbi")

    """

        function perform_imputation {
            java -jar `which beagle.jar` chrom=\${1} window=8000 overlap=3000 impute=true ne=17500 gt=WI.${date}.hard-filter.vcf.gz out=\${1}
        }

        export -f perform_imputation

        parallel --verbose perform_imputation {} ::: I II III IV V X MtDNA

        bcftools concat I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz > WI.${date}.impute.vcf.gz
        bcftools index --threads=${task.cpus} WI.${date}.impute.vcf.gz
        tabix WI.${date}.impute.vcf.gz
        bcftools stats --verbose WI.${date}.impute.vcf.gz > WI.${date}.impute.stats.txt
    """
}


impute_vcf.into { kinship_vcf;  mapping_vcf; }


process make_kinship {

    publishDir "${params.out}/cegwas", mode: 'copy'

    input:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") from kinship_vcf
    output:
        file("kinship.Rda")

    """
        Rscript -e 'library(cegwas); kinship <- generate_kinship("WI.${date}.impute.vcf.gz"); save(kinship, file = "kinship.Rda");'
    """

}


process make_mapping_rda_file {

    publishDir "${params.out}/cegwas", mode: 'copy'

    input:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") from mapping_vcf
    output:
        file("snps.Rda")

    """
        Rscript -e 'library(cegwas); snps <- generate_mapping("WI.${date}.impute.vcf.gz"); save(snps, file = "snps.Rda");'
    """

}


/*
    Haplotype analysis
*/

process_ibd=file("process_ibd.R")

minalleles = 0.05 // Species the minimum number of samples carrying the minor allele.
r2window = 1500 // Specifies the number of markers in the sliding window used to detect correlated markers.
ibdtrim = 0
r2max = 0.8

process ibdseq {

    publishDir "${params.out}/haplotype", mode: 'copy'

    tag { "ibd" }

    input:
        set file("WI.${date}.hard-filter-biallelic-snp.vcf.gz"), file("WI.${date}.hard-filter-biallelic-snp.vcf.gz.csi") from haplotype_vcf

    output:
        file("haplotype.tsv") into haplotype_analysis

    """
    minalleles=\$(bcftools query --list-samples WI.${date}.hard-filter-biallelic-snp.vcf.gz | wc -l | awk '{ print \$0*${minalleles} }' | awk '{printf("%d\\n", \$0+=\$0<0?0:0.9)}')
    if [[ \${minalleles} -lt 2 ]];
    then
        minalleles=2;
    fi;
    echo "minalleles=${minalleles}"
    for chrom in I II III IV V X; do
        java -jar `which ibdseq.r1206.jar` \\
            gt=WI.${date}.hard-filter-biallelic-snp.vcf.gz \\
            out=haplotype_\${chrom} \\
            ibdtrim=${ibdtrim} \\
            minalleles=\${minalleles} \\
            r2max=${r2max} \\
            nthreads=4 \\
            chrom=\${chrom}
        done;
    cat *.ibd | awk '{ print \$0 "\\t${minalleles}\\t${ibdtrim}\\t${r2window}\\t${r2max}" }' > haplotype.tsv
    """
}

process analyze_ibdseq {

    publishDir "${params.out}/haplotype", mode: 'copy'

    input:
        file("haplotype.tsv") from haplotype_analysis

    output:
        file("processed_haps.Rda")
        file("haplotype_plot_df.Rda") into plot_df


    """
        Rscript --vanilla `which process_ibd.R`
    """
}

process plot_ibdseq {

    publishDir "${params.out}/haplotype", mode: 'copy'

    input:
        file("haplotype_plot_df.Rda") from plot_df

    output:
        file("haplotype_length.png")
        file("max_haplotype_sorted_genome_wide.png")
        file("haplotype.png")
        file("sweep_summary.tsv")

    """
        Rscript --vanilla `which plot_ibd.R`
    """
}


mod_tracks = Channel.from(["LOW", "MODERATE", "HIGH", "MODIFIER"])
soft_filter_vcf_mod_tracks.spread(mod_tracks).set { mod_track_set }


process generate_mod_tracks {

    publishDir "${params.out}/tracks", mode: 'copy'

    tag { severity }

    input:
        set file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), val(severity) from mod_track_set
    output:
        set file("${date}.${severity}.bed.gz"), file("${date}.${severity}.bed.gz.tbi")

    """
        bcftools view --apply-filters PASS WI.${date}.vcf.gz | \
        grep ${severity} | \
        awk '\$0 !~ "^#" { print \$1 "\\t" (\$2 - 1) "\\t" (\$2)  "\\t" \$1 ":" \$2 "\\t0\\t+"  "\\t" \$2 - 1 "\\t" \$2 "\\t0\\t1\\t1\\t0" }' | \\
        bgzip  > ${date}.${severity}.bed.gz
        tabix -p bed ${date}.${severity}.bed.gz
    """
}


sample_files_list.splitText() { it.strip() } .combine(soft_filter_vcf_strain).into { isotype_set_vcf; isotype_set_tsv }


process generate_isotype_vcf {

    publishDir "${params.out}/isotype/vcf", mode: 'copy'

    tag { isotype }

    input:
        set val(isotype), file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi") from isotype_set_vcf

    output:
        set val(isotype), file("${isotype}.${date}.vcf.gz"), file("${isotype}.${date}.vcf.gz.tbi") into isotype_ind_vcf

    """
        bcftools view -O z --samples ${isotype} --exclude-uncalled WI.${date}.vcf.gz  > ${isotype}.${date}.vcf.gz && tabix ${isotype}.${date}.vcf.gz
    """

}


process generate_isotype_tsv {

    publishDir "${params.out}/isotype/tsv", mode: 'copy'

    tag { isotype }

    input:
        set val(isotype), file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi") from isotype_set_tsv

    output:
        set val(isotype), file("${isotype}.${date}.tsv.gz")

    """
        echo 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\tFT\\tGT' > ${isotype}.${date}.tsv
        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\t%FILTER\\t%FT\\t%TGT]\\n' --samples ${isotype} WI.${date}.vcf.gz > ${isotype}.${date}.tsv
        bgzip ${isotype}.${date}.tsv
        tabix -S 1 -s 1 -b 2 -e 2 ${isotype}.${date}.tsv.gz
    """

}

vcf_stats = soft_filter_stats.concat ( hard_filter_stats, impute_stats )

process multiqc_report {

    executor 'local'

    publishDir "${params.out}/report", mode: 'copy'

    input:
        file("stat*") from vcf_stats.toSortedList()
        file(samtools_stats) from SM_samtools_stats_set.toSortedList()
        //file(bamtools_stats) from SM_bamtools_stats_set.toSortedList()
        file(duplicates) from duplicates_set.toSortedList()
        file(fastqc) from SM_fastqc_stats_set.toSortedList()
        file("bam*.idxstats") from bam_idxstats_multiqc.toSortedList()
        file("picard*.stats.txt") from SM_picard_stats_set.collect()
        file("snpeff_out.csv") from snpeff_multiqc
        file("WI.${date}.hard-filter-biallelic-snp.stats.txt") from biallelic_snp_stats
        //file("DELLYsv_snpeff_out.csv") from snpeff_delly_multiqc
        //file("MANTAsv_snpeff_out.csv") from manta_snpeff_multiqc
        //file("WI.${date}.DELLYsv.raw.stats.txt") from delly_bcf_stats
        //file("WI.${date}.MANTAsv.soft-filter.stats.txt") from bcf_manta_stats
        //file("WI.${date}.TIDDITsv.soft-filter.stats.txt") from tiddit_stats
        //file("TIDDITsv_snpeff_out.csv") from snpeff_tiddit_multiqc

    output:
        file("multiqc_data/*.json") into multiqc_json_files
        file("multiqc.html")

    """
        multiqc -k json --filename multiqc.html .
    """

}

process comprehensive_report {


    input:
        file("multiqc_data/*.json") from multiqc_json_files

    """
        echo "great"
        exit 0
    """

}

workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }

    // mail summary
    if (params.email) {
        ['mail', '-s', 'wi-nf', params.email].execute() << summary
    }


}