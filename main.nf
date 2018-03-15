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
params.annotation_reference = "WS261"
params.cores = 6
params.tmpdir = "tmp/"
params.email = ""
params.reference = "(required)"
params.manta_path = null
params.tiddit_discord = null
params.snpeff_path="${workflow.workDir}/snpeff"


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

} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fq_file_prefix = null;
    params.fqs = "sample_sheet.tsv"
}

File fq_file = new File(params.fqs);

/*
    ==
    UX
    ==
*/

param_summary = '''


     ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄                         ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄
    ▐░▌       ▐░▌▐░░░░░░░░░░░▌                       ▐░░▌      ▐░▌▐░░░░░░░░░░░▌
    ▐░▌       ▐░▌ ▀▀▀▀█░█▀▀▀▀                        ▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌       ▐░▌     ▐░▌                            ▐░▌▐░▌    ▐░▌▐░▌
    ▐░▌   ▄   ▐░▌     ▐░▌           ▄▄▄▄▄▄▄▄▄▄▄      ▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄
    ▐░▌  ▐░▌  ▐░▌     ▐░▌          ▐░░░░░░░░░░░▌     ▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌
    ▐░▌ ▐░▌░▌ ▐░▌     ▐░▌           ▀▀▀▀▀▀▀▀▀▀▀      ▐░▌   ▐░▌ ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌▐░▌ ▐░▌▐░▌     ▐░▌                            ▐░▌    ▐░▌▐░▌▐░▌
    ▐░▌░▌   ▐░▐░▌ ▄▄▄▄█░█▄▄▄▄                        ▐░▌     ▐░▐░▌▐░▌
    ▐░░▌     ▐░░▌▐░░░░░░░░░░░▌                       ▐░▌      ▐░░▌▐░▌
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

    publishDir params.out + "/phenotype", mode: 'copy'

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

    publishDir params.out + "/alignment", mode: 'copy'

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

    publishDir params.out + "/alignment", mode: 'copy'

    input:
        val stat_files from SM_bam_stat_files.toSortedList()

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

    publishDir params.out + "/alignment", mode: 'copy'

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

    publishDir params.out + "/phenotype", mode: 'copy'

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

    publishDir params.out + "/phenotype", mode: 'copy'

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
    ====================
    Call Variants - GATK
    ====================
    
    Generate SM gVCFs
*/
process call_variants_individual {

    tag { SM }

    cpus params.cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_snp_individual

    output:
        file("${SM}.g.vcf") into individual_sites
        file("${SM}.g.vcf.idx") into individual_sites_index

    """

    function split_gatk() {
      gatk-launch HaplotypeCaller \\
            -R ${reference_handle_uncompressed} \\
            -I ${SM}.bam \\
            --emit-ref-confidence GVCF \\
            --sample-ploidy 1 \\
            --genotyping-mode DISCOVERY \\
            --max-genotype-count 3000 \\
            --max-alternate-alleles 100 \\
            --annotation DepthPerAlleleBySample \\
            --annotation Coverage \\
            --annotation GenotypeSummaries \\
            --annotation TandemRepeat \\
            --annotation StrandBiasBySample \\
            --annotation ChromosomeCounts \\
            --annotation AS_QualByDepth \\
            --annotation AS_StrandOddsRatio \\
            --annotation AS_MappingQualityRankSumTest \\
            --annotation DepthPerSampleHC \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --annotation-group StandardHCAnnotation \\
            -L \${1} \\
            -O ${SM}_\${1}.g.vcf            
    }

    export -f split_gatk

    parallel --verbose split_gatk {} ::: I II III IV V X MtDNA

    gatk-launch CombineGVCFs \\
        -R ${reference_handle_uncompressed} \\
        --variant ${SM}_I.g.vcf \\
        --variant ${SM}_II.g.vcf \\
        --variant ${SM}_III.g.vcf \\
        --variant ${SM}_IV.g.vcf \\
        --variant ${SM}_V.g.vcf \\
        --variant ${SM}_X.g.vcf \\
        --variant ${SM}_MtDNA.g.vcf \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        -O ${SM}.g.vcf 

    gatk-launch IndexFeatureFile \\
        -F ${SM}.g.vcf

    """
}

// Merge gVCFs

process merge_gvcfs {

    cpus params.cores

    input:
        file gvcfs from individual_sites.toSortedList()
        file indices from individual_sites_index.toSortedList()

    output:
        set file("wild_isolate.vcf.gz"), file("wild_isolate.vcf.gz.tbi") into raw_wild_isolate_vcf


    """
    find . -name '*.g.vcf' > input.list

    function combine_gvcfs() {
    gatk-launch GenomicsDBImport \\
        -R ${reference_handle_uncompressed} \\
            -V input.list \\
            --genomicsdb-workspace-path my_database_\${1} \\
            -L \${1}

    gatk-launch GenotypeGVCFs \\
             -R ${reference_handle_uncompressed} \\
             -new-qual \\
             -stand-call-conf 0 \\
            --max-genotype-count 3000 \\
            --max-alternate-alleles 100 \\
             --annotation DepthPerAlleleBySample \\
             --annotation Coverage \\
             --annotation GenotypeSummaries \\
             --annotation TandemRepeat \\
             --annotation StrandBiasBySample \\
             --annotation ChromosomeCounts \\
             --annotation AS_QualByDepth \\
             --annotation AS_StrandOddsRatio \\
             --annotation AS_MappingQualityRankSumTest \\
             --annotation DepthPerSampleHC \\
             --variant gendb://my_database_\${1} \\
             --annotation-group StandardAnnotation \\
             --annotation-group AS_StandardAnnotation \\
             -L \${1} \\
             -O wild_isolate_\${1}.vcf
    }

    export -f combine_gvcfs

    parallel --verbose combine_gvcfs {} ::: I II III IV V X MtDNA

    gatk-launch GatherVcfs \\
        -R ${reference_handle_uncompressed} \\
        -I wild_isolate_I.vcf \\
        -I wild_isolate_II.vcf \\
        -I wild_isolate_III.vcf \\
        -I wild_isolate_IV.vcf \\
        -I wild_isolate_V.vcf \\
        -I wild_isolate_X.vcf \\
        -I wild_isolate_MtDNA.vcf \\
        -O wild_isolate.vcf 

    bgzip --threads=${task.cpus} -c wild_isolate.vcf > wild_isolate.vcf.gz
    tabix -p vcf wild_isolate.vcf.gz

    """
}

process gatk_to_diploid {

    publishDir "${params.out}/variation", mode: 'copy'

    input:
      set file("wild_isolate.vcf.gz"), file("wild_isolate.vcf.gz.tbi") from raw_wild_isolate_vcf

    output:
      set file("WI.${date}.raw.vcf.gz"), file("WI.${date}.raw.vcf.gz.tbi") into raw_wild_isolate_diploid_vcf


    """
        bcftools view wild_isolate.vcf.gz | \\
        sed -E 's/\\t([0-9\\.]+):/\\t\\1\\/\\1:/g' | \\
        bgzip -c --threads ${task.cpus} > WI.${date}.raw.vcf.gz
        tabix -p vcf WI.${date}.raw.vcf.gz
    """
}

// To save some steps in this process we can incorporate TYPE in info field, this will prevent the need to split and apply indel filters 

process apply_filters {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
        set file(unionvcf), file(unionvcfindex) from raw_wild_isolate_diploid_vcf

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf
        set val('filtered'), file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf_stat
        file("WI.${date}.soft-filter.stats.txt") into soft_filter_stats

    """
        # Implement bash-trap; This removes files whether the process succeeds or fails.
        function finish {
            rm -f wi_norm.vcf.gz
            rm -f indel_soft_filters.vcf.gz
            rm -f snps.vcf.gz
            rm -f snp_indel_soft_filters.vcf
            rm -f wi_qual.vcf
        }
        trap finish EXIT


        gatk-launch VariantFiltration \\
            -R ${reference_handle_uncompressed} \\
            --variant ${unionvcf} \\
            --genotype-filter-expression "DP < ${params.min_depth}" \\
            --genotype-filter-name "depth" \\
            -O wi_dp.vcf

        bcftools norm -m -any --threads ${task.cpus} -O z wi_dp.vcf > wi_norm.vcf.gz 

        # Output indels
        bcftools view -v indels wi_norm.vcf.gz | \\
        bcftools filter -O v --threads ${task.cpus-1} --mode + --soft-filter indelsor --include "INFO/SOR > ${params.strand_odds_ratio}" | \\
        bcftools filter -O z --threads ${task.cpus-1} --mode + --soft-filter indelqd --include "INFO/QD > ${params.quality_by_depth}" > indel_soft_filters.vcf.gz 
        bcftools index --threads=${task.cpus} indel_soft_filters.vcf.gz

        # Output snps
        bcftools view -v snps -O z wi_norm.vcf.gz > snps.vcf.gz
        bcftools index --threads=${task.cpus} snps.vcf.gz

        bcftools concat --threads ${task.cpus-1} \\
                        --allow-overlaps \\
                        indel_soft_filters.vcf.gz \\
                        snps.vcf.gz | \\
        bcftools filter -O v \\
                        --mode + \\
                        --soft-filter mapping_quality \\
                        --include "INFO/MQ > ${params.mapping_quality}" > snp_indel_soft_filters.vcf

        gatk-launch VariantFiltration \\
            -R ${reference_handle_uncompressed} \\
            --variant snp_indel_soft_filters.vcf \\
            --genotype-filter-expression "( AD[1] / (AD[0] + AD[0]) ) < ${params.dv_dp}" \\
            --genotype-filter-name "dv_dp" \\
            --genotype-filter-expression "QD < 10.0 && AD[1] / (AD[1] + AD[0]) < ${params.dv_dp} && ReadPosRankSum < 0.0" \\
            --genotype-filter-name "dv_dp_qd_readposranksum" \\
            -O wi_dv_dp.vcf

        bcftools norm -m +any -O v -o wi_dv_dp_norm.vcf wi_dv_dp.vcf

        gatk-launch VariantFiltration \\
            -R ${reference_handle_uncompressed} \\
            --variant wi_dv_dp_norm.vcf \\
            --filter-expression "QUAL < ${params.qual}" \\
            --filter-name "quality" \\
            -O wi_qual.vcf

        bcftools filter -O z --threads ${task.cpus-1} --mode + --soft-filter high_missing --include "F_MISSING<=${params.missing}" wi_qual.vcf > WI.${date}.soft-filter.vcf.gz
        bcftools index --threads ${task.cpus} -f WI.${date}.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt
    """
}


filtered_vcf.into {
                    filtered_vcf_snpeff;
                    filtered_vcf_to_hard;
                    filtered_vcf_gtcheck;
                    filtered_vcf_primer;
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

process annotate_vcf_snpeff {

    cpus params.cores

    errorStrategy 'retry'
    maxRetries 2

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_snpeff
        file("gene.pkl") from gene_pkl_snpindel

    output:
        set file("WI.${date}.snpeff.vcf.gz"), file("WI.${date}.snpeff.vcf.gz.csi") into snpeff_vcf
        file("snpeff_out.csv") into snpeff_multiqc

    script:
        """
            bcftools view --threads=${params.cores-1} -O v WI.${date}.soft-filter.vcf.gz | \\
            snpEff eff -csvStats snpeff_out.csv \\
            -no-downstream -no-intergenic -no-upstream \\
            -dataDir ${params.snpeff_path} \\
            -config ${params.snpeff_path}/snpEff.config \\
            ${params.annotation_reference} | \\
            bcftools view -O v | \\
            python `which fix_snpeff_names.py` - | \\
            bcftools view --threads=${task.cpus-1} -O z > WI.${date}.snpeff.vcf.gz
            bcftools index --threads=${task.cpus} WI.${date}.snpeff.vcf.gz
        """

}

process generate_hard_vcf {

    cpus params.cores

    publishDir params.out + "/variation", mode: 'copy'

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_to_hard

    output:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf_to_impute
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into tajima_bed
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into vcf_phylo
        set val('clean'), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf
        file("WI.${date}.hard-filter.vcf.gz.tbi")
        file("WI.${date}.hard-filter.stats.txt") into hard_filter_stats


    """
        # Generate hard-filtered (clean) vcf
        bcftools view WI.${date}.soft-filter.vcf.gz | \\
        bcftools filter --set-GTs . --exclude 'FORMAT/FT != "PASS"' | \\
        vk filter MISSING --max=${params.missing} - | \\
        vk filter HET --max=0.10 - | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        vcffixup - | \\
        bcftools view --trim-alt-alleles -O z > WI.${date}.hard-filter.vcf.gz
        bcftools index -f WI.${date}.hard-filter.vcf.gz
        tabix WI.${date}.hard-filter.vcf.gz
        bcftools stats --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt
    """
}

hard_vcf.set { hard_vcf_summary }

process calculate_gtcheck {

    publishDir params.out + "/concordance", mode: 'copy'

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_gtcheck

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

    publishDir params.out + "/variation", mode: 'copy'

    input:
        set val('clean'), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from hard_vcf_summary

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
        file("${contig}.pdf")
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
        set file("WI.${date}.tajima.bed.gz"), file("WI.${date}.tajima.bed.gz.tbi")

    """
        vk tajima --no-header 100000 10000 WI.${date}.hard-filter.vcf.gz | bgzip > WI.${date}.tajima.bed.gz
        tabix WI.${date}.tajima.bed.gz
    """

}


process imputation {

    cpus params.cores

    publishDir params.out + "/variation", mode: 'copy'


    input:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") from hard_vcf_to_impute
    output:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") into impute_vcf
        file("WI.${date}.impute.stats.txt") into impute_stats
        file("WI.${date}.impute.stats.txt") into filtered_stats
        file("WI.${date}.impute.vcf.gz")
        file("WI.${date}.impute.vcf.gz.tbi")

    """
        beagle nthreads=${task.cpus} window=8000 overlap=3000 impute=true ne=17500 gt=WI.${date}.hard-filter.vcf.gz out=WI.${date}.impute
        bcftools index --threads=${task.cpus} WI.${date}.impute.vcf.gz
        tabix WI.${date}.impute.vcf.gz
        bcftools stats --verbose WI.${date}.impute.vcf.gz > WI.${date}.impute.stats.txt
    """
}


impute_vcf.into { kinship_vcf;  mapping_vcf; haplotype_vcf }


process make_kinship {

    publishDir params.out + "/cegwas", mode: 'copy'

    input:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") from kinship_vcf
    output:
        file("kinship.Rda")

    """
        Rscript -e 'library(cegwas); kinship <- generate_kinship("WI.${date}.impute.vcf.gz"); save(kinship, file = "kinship.Rda");'
    """

}


process make_mapping_rda_file {

    publishDir params.out + "/cegwas", mode: 'copy'

    input:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") from mapping_vcf
    output:
        file("snps.Rda")

    """
        Rscript -e 'library(cegwas); snps <- generate_mapping("WI.${date}.impute.vcf.gz"); save(snps, file = "snps.Rda");'
    """

}


process download_annotation_files {

    executor 'local'

    errorStrategy 'retry'
    maxRetries 5

    output:
        set val("phastcons"), file("elegans.phastcons.wib") into phastcons
        set val("phylop"), file("elegans.phylop.wib") into phylop
        set val("repeatmasker"), file("elegans_repeatmasker.bb") into repeatmasker
        file("Caenorhabditis_elegans.WBcel235.91.gff3.gz") into ensembl_gff3_csq

    """
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phastcons.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans.phylop.wib
        wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS258/MULTI_SPECIES/hub/elegans/elegans_repeatmasker.bb
        wget ftp://ftp.ensembl.org/pub/current_gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.91.gff3.gz
    """
}

phastcons.mix(phylop).set { wig }

process wig_to_bed {

    tag { track_name }

    publishDir params.out + '/tracks', mode: 'copy'

    input:
        set val(track_name), file("track.wib") from wig
    output:
        file("${track_name}.bed.gz") into bed_tracks
        file("${track_name}.bed.gz.tbi") into bed_indices

    """
        bigWigToBedGraph track.wib ${track_name}.bed
        bgzip ${track_name}.bed
        tabix ${track_name}.bed.gz
    """

}


process annovar_and_output_soft_filter_vcf {

    publishDir params.out + "/variation", mode: 'copy'

    cpus params.cores

    input:
        set file("WI.${date}.snpeff.vcf.gz"), file("WI.${date}.snpeff.vcf.gz.csi") from snpeff_vcf
        file(track) from bed_tracks.toSortedList()
        file(track) from bed_indices.toSortedList()
        file('vcf_anno.conf') from Channel.fromPath("vcfanno.conf")
        file("ensembl.gff3.gz") from ensembl_gff3_csq

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi"), file("WI.${date}.soft-filter.vcf.gz.tbi") into soft_filter_vcf
        file("WI.${date}.soft-filter.stats.txt")

    """
        vcfanno -p ${task.cpus} vcf_anno.conf WI.${date}.snpeff.vcf.gz | \\
        bcftools view -O z > WI.${date}.soft-filter.vcf.gz
        bcftools index WI.${date}.soft-filter.vcf.gz
        tabix WI.${date}.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt
    """

}

soft_filter_vcf.into {
                        soft_filter_vcf_strain;
                        soft_filter_vcf_isotype_list;
                        soft_filter_vcf_mod_tracks;
                        soft_filter_vcf_tsv
                     }


mod_tracks = Channel.from(["LOW", "MODERATE", "HIGH", "MODIFIER"])
soft_filter_vcf_mod_tracks.spread(mod_tracks).set { mod_track_set }


process generate_mod_tracks {

    publishDir params.out + '/tracks', mode: 'copy'

    tag { severity }

    input:
        set file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), file("WI.${date}.vcf.gz.tbi"), val(severity) from mod_track_set
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

process generate_strain_list {

    executor 'local'

    input:
        set file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), file("WI.${date}.vcf.gz.tbi") from soft_filter_vcf_isotype_list

    output:
        file('isotype_list.tsv') into isotype_list

    """
        bcftools query -l WI.${date}.vcf.gz > isotype_list.tsv
    """

}


isotype_list.splitText() { it.strip() } .spread(soft_filter_vcf_strain).into { isotype_set_vcf; isotype_set_tsv }


process generate_isotype_vcf {

    publishDir params.out + '/isotype/vcf', mode: 'copy'

    tag { isotype }

    input:
        set val(isotype), file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), file("WI.${date}.vcf.gz.tbi") from isotype_set_vcf

    output:
        set val(isotype), file("${isotype}.${date}.vcf.gz"), file("${isotype}.${date}.vcf.gz.tbi") into isotype_ind_vcf

    """
        bcftools view -O z --samples ${isotype} --exclude-uncalled WI.${date}.vcf.gz  > ${isotype}.${date}.vcf.gz && tabix ${isotype}.${date}.vcf.gz
    """

}


process generate_isotype_tsv {

    publishDir params.out + '/isotype/tsv', mode: 'copy'

    tag { isotype }

    input:
        set val(isotype), file("WI.${date}.vcf.gz"), file("WI.${date}.vcf.gz.csi"), file("WI.${date}.vcf.gz.tbi") from isotype_set_tsv

    output:
        set val(isotype), file("${isotype}.${date}.tsv.gz")

    """
        echo 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\tFT\\tGT' > ${isotype}.${date}.tsv
        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\t%FILTER\\t%FT\\t%TGT]\\n' --samples ${isotype} WI.${date}.vcf.gz > ${isotype}.${date}.tsv
        bgzip ${isotype}.${date}.tsv
        tabix -S 1 -s 1 -b 2 -e 2 ${isotype}.${date}.tsv.gz
    """

}

vcf_stats = soft_filter_stats.concat( hard_filter_stats, impute_stats )


/*
    Manta-sv
*/

process manta_call {
    
    tag { SM }

    
    when:
        params.manta_path

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_manta

    output:
        file "*.vcf.gz" into individual_output_vcf_zipped
        file "*.vcf.gz.tbi" into individual_output_index
        set val(SM), file("${SM}_manta.vcf") into manta_to_db


    """
        configManta.py \\
        --bam ${SM}.bam \\
        --referenceFasta ${reference_handle_uncompressed} \\
        --outputContig \\
        --runDir results/

        python results/runWorkflow.py -m local -j 8

        cp results/results/variants/diploidSV.vcf.gz .
        cp results/results/variants/diploidSV.vcf.gz.tbi .

        mv diploidSV.vcf.gz ${SM}_manta.vcf.gz 
        mv diploidSV.vcf.gz.tbi ${SM}_manta.vcf.gz.tbi

        bcftools view -Ov -o ${SM}_manta.vcf ${SM}_manta.vcf.gz 
    """

}

individual_output_vcf_zipped
  .toSortedList()
  .set { merged_deletion_vcf }


individual_output_index
  .toSortedList()
  .set { merged_vcf_index }


process merge_manta_vcf {

    publishDir params.out + "/variation", mode: 'copy'

    when:
        params.manta_path

    input:
      file merged_deletion_vcf
      file merged_vcf_index

    output:
      set file("WI.${date}.MANTAsv.soft-filter.vcf.gz"), file("WI.${date}.MANTAsv.soft-filter.vcf.gz.csi") into processed_manta_vcf
      file("WI.${date}.MANTAsv.soft-filter.stats.txt") into bcf_manta_stats

    """
        bcftools merge -m all --threads ${task.cpus-1} -o WI.${date}.MANTAsv.soft-filter.vcf.gz -Oz ${merged_deletion_vcf}
        bcftools index --threads ${task.cpus} -f WI.${date}.MANTAsv.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.MANTAsv.soft-filter.vcf.gz > WI.${date}.MANTAsv.soft-filter.stats.txt
    """

}

process prune_manta {
    
    publishDir params.out + "/variation", mode: 'copy'

    input:
        set file(mantavcf), file(mantaindex) from processed_manta_vcf
        file("gene.pkl") from gene_pkl_manta

    output:
        set file("WI.${date}.MANTAsv.LargeRemoved.snpeff.vcf.gz"), file("WI.${date}.MANTAsv.LargeRemoved.snpeff.vcf.gz.csi") into snpeff_manta_vcf
        file("WI.${date}.MANTAsv.CONTIGS.tsv.gz") into manta_contigs
        file("MANTAsv_snpeff_out.csv") into manta_snpeff_multiqc

    """
        bcftools plugin setGT -Oz -o manta_gt_filled.vcf.gz -- ${mantavcf} -t . -n 0
        bcftools query -l manta_gt_filled.vcf.gz | sort > sample_names.txt
        bcftools view --samples-file=sample_names.txt -Oz -o manta_gt_filled_sorted.vcf.gz manta_gt_filled.vcf.gz

        bcftools view manta_gt_filled_sorted.vcf.gz | \\
        bcftools filter -e 'INFO/SVLEN>100000' | \\
        bcftools filter -e 'INFO/SVLEN<-100000' | \\
        bcftools view -Oz -o manta_gt_filled_sorted_largeRemoved.vcf.gz

        bcftools view -O v manta_gt_filled_sorted_largeRemoved.vcf.gz | \\
        snpEff eff -csvStats MANTAsv_snpeff_out.csv \\
        -no-downstream -no-intergenic -no-upstream \\
        -dataDir ${params.snpeff_path} \\
        -config ${params.snpeff_path}/snpEff.config \\
        ${params.annotation_reference} | \\
        bcftools view -O v | \\
        python `which fix_snpeff_names.py` - | \\
        bcftools view -O z > WI.${date}.MANTAsv.LargeRemoved.snpeff.vcf.gz

        bcftools index -f WI.${date}.MANTAsv.LargeRemoved.snpeff.vcf.gz

        bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVTYPE\\t%SVLEN\\t%CONTIG[\\t%GT]\\n' WI.${date}.MANTAsv.LargeRemoved.snpeff.vcf.gz > WI.${date}.MANTAsv.CONTIGS.tsv

        bgzip  WI.${date}.MANTAsv.CONTIGS.tsv
    """

}

/*
    ========
    Delly-sv
    ========
*/

process delly_sv {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_delly

    output:
        file "*.bcf" into dellybcf


    """
        delly call ${SM}.bam -g ${reference_handle_uncompressed} -o ${SM}.bcf
    """

}

dellybcf
    .toSortedList()
    .into{ deletion_bcf }

process combine_delly {
      
  input:
    file deletion_bcf

  output:
    file "*.bcf" into combined_delly_bcf


  script:
    """
        delly merge ${deletion_bcf} -m 100 -n 100000 -b 500 -r 0.5 -o WI.delly.first.bcf
    """

}

process recall_deletions {
        
    tag { SM }
    
    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_delly_recall
        file combined_delly_bcf

    output:
        file "${SM}_second.bcf" into recalled_delly_sv
        file "${SM}_second.bcf.csi" into recalled_delly_sv_index
        set val(SM), file("${SM}_second.vcf") into delly_to_db


    """
        delly call ${SM}.bam -v ${combined_delly_bcf} -g ${reference_handle_uncompressed} -o ${SM}_second.bcf
        bcftools view -Ov -o ${SM}_second.vcf ${SM}_second.bcf
    """

}

recalled_delly_sv
  .toSortedList()
  .into{ recalled_delly_sv_bcf }

recalled_delly_sv_index
  .toSortedList()
  .into{ recalled_delly_sv_bcf_index }

process combine_second_deletions {
        
    publishDir params.out + "/variation", mode: 'copy'

    input:
        file recalled_delly_sv_bcf
        file recalled_delly_sv_bcf_index

    output:
        set file("WI.${date}.DELLYsv.raw.vcf.gz"), file("WI.${date}.DELLYsv.raw.vcf.gz.csi") into raw_recalled_wi_dell_sv
        set file("WI.${date}.DELLYsv.germline-filter.vcf.gz"), file("WI.${date}.DELLYsv.germline-filter.vcf.gz.csi") into germline_recalled_wi_dell_sv
        file "WI.${date}.DELLYsv.raw.stats.txt" into delly_bcf_stats

    script:
        """
            bcftools merge -m id -O b -o WI.${date}.DELLYsv.raw.bcf ${recalled_delly_sv_bcf}
            bcftools index -f WI.${date}.DELLYsv.raw.bcf
            bcftools query -l WI.${date}.DELLYsv.raw.bcf | sort > sample_names.txt
            bcftools view --samples-file=sample_names.txt -Oz -o WI.${date}.DELLYsv.raw.vcf.gz WI.${date}.DELLYsv.raw.bcf
            bcftools index -f WI.${date}.DELLYsv.raw.vcf.gz

            delly filter -f germline WI.${date}.DELLYsv.raw.bcf -o WI.${date}.DELLYsv.germline-filter.bcf

            bcftools index -f WI.${date}.DELLYsv.germline-filter.bcf
            bcftools view --samples-file=sample_names.txt -Oz -o WI.${date}.DELLYsv.germline-filter.vcf.gz WI.${date}.DELLYsv.germline-filter.bcf
            bcftools index -f WI.${date}.DELLYsv.germline-filter.vcf.gz

            bcftools stats --verbose WI.${date}.DELLYsv.raw.vcf.gz > WI.${date}.DELLYsv.raw.stats.txt     
        """
}

process delly_snpeff {
    
    publishDir params.out + "/variation", mode: 'copy'

    input:
        set file(dellysv), file(dellysvindex) from germline_recalled_wi_dell_sv
        file("gene.pkl") from gene_pkl_manta

    output:
        set file("WI.${date}.DELLYsv.snpEff.vcf.gz"), file("WI.${date}.DELLYsv.snpEff.vcf.gz.csi") into snpeff_delly_vcf
        file("DELLYsv_snpeff_out.csv") into snpeff_delly_multiqc

    when:
        params.tiddit_discord


    script:
      """
        bcftools view --threads ${task.cpus-1} ${dellysv} | \\
        snpEff eff -csvStats DELLYsv_snpeff_out.csv \\
        -no-downstream \\
        -no-intergenic \\
        -no-upstream \\
        -dataDir ${params.snpeff_path} \\
        -config ${params.snpeff_path}/snpEff.config \\
        ${params.annotation_reference} | \\
        bcftools view -O v | \\
        python `which fix_snpeff_names.py` - | \\
        bcftools view -O z > WI.${date}.DELLYsv.snpEff.vcf.gz

        bcftools index --threads ${task.cpus} -f WI.${date}.DELLYsv.snpEff.vcf.gz
      """
}

/*
    =========
    TIDDIT-sv
    =========
*/

process tiddit_call_sv {
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_tiddit

    output:
        file "${SM}_tiddit.vcf.gz" into tiddit_vcf
        file "${SM}_tiddit.vcf.gz.csi" into tiddit_index
        file "*.signals.tab" into tiddit_coverage
        set val(SM), file("${SM}_tiddit.vcf") into tiddit_to_db

    when:
        params.tiddit_discord

    """
        python2 ${params.tiddit} \\
        --sv \\
        -o ${SM}_tiddit \\
        -p ${params.tiddit_discord} \\
        -r ${params.tiddit_discord} \\
        --bam ${SM}.bam \\
        --ref ${reference_handle_uncompressed}

        bcftools view ${SM}_tiddit.vcf | \\
        vk geno transfer-filter - | \\
        bcftools view --threads ${task.cpus-1} -O z > ${SM}_tiddit.vcf.gz

        bcftools index --threads ${task.cpus} -f ${SM}_tiddit.vcf.gz
    """

}

tiddit_vcf
    .toSortedList()
    .set { tiddit_sample_vcfs }

tiddit_index
    .toSortedList()
    .set { tiddit_sample_indices }

tiddit_coverage
    .toSortedList()
    .set { tiddit_sample_coverage }


process merge_tiddit_vcf {
    
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file tiddit_sample_vcfs
        file tiddit_sample_coverage
        file tiddit_sample_indices

    output:
        set file("WI.${date}.TIDDITsv.soft-filter.vcf.gz"), file("WI.${date}.TIDDITsv.soft-filter.vcf.gz.csi") into softfilter_tiddit_vcf
        file("WI.${date}.TIDDITsv.soft-filter.stats.txt") into tiddit_stats

    when:
        params.tiddit_discord

    """
        bcftools merge -m all --threads ${task.cpus-1} -Ov ${tiddit_sample_vcfs} | \\
        bcftools filter --threads ${task.cpus-1} --set-GTs . --exclude 'FORMAT/FT != "PASS"' -O z | \\
        bcftools view -O z \\
                --samples-file=<(bcftools query -l WI.${date}.TIDDITsv.soft-filter_unsorted.vcf.gz | sort) \\
                 WI.${date}.TIDDITsv.soft-filter_unsorted.vcf.gz > WI.${date}.TIDDITsv.soft-filter.vcf.gz

        bcftools index --threads ${task.cpus} -f WI.${date}.TIDDITsv.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.TIDDITsv.soft-filter.vcf.gz > WI.${date}.TIDDITsv.soft-filter.stats.txt
    """

}

process tiddit_snpeff {
        
    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file(tiddit_joint_vcf), file(tiddit_joint_index) from softfilter_tiddit_vcf
        file("gene.pkl") from gene_pkl_tiddit

    output:
        set file("WI.${date}.TIDDITsv.snpEff.vcf.gz"), file("WI.${date}.TIDDITsv.snpEff.vcf.gz.csi") into snpeff_tiddit_vcf
        file("TIDDITsv_snpeff_out.csv") into snpeff_tiddit_multiqc
    
    when:
        params.tiddit_discord

    """
        bcftools view ${tiddit_joint_vcf} | \\
        snpEff eff -csvStats TIDDITsv_snpeff_out.csv \\
        -no-downstream -no-intergenic -no-upstream \\
        -dataDir ${params.snpeff_path} \\
        -config ${params.snpeff_path}/snpEff.config \\
        ${params.annotation_reference} | \\
        bcftools view -O v | \\
        python `which fix_snpeff_names.py` - | \\
        bcftools view -O z > WI.${date}.TIDDITsv.snpEff.vcf.gz

        bcftools index -f WI.${date}.TIDDITsv.snpEff.vcf.gz
    """
}


manta_to_db
    .join(delly_to_db)
    .join(tiddit_to_db)
    .set { variant_db }


process merge_sv_callers {

    input:
        set val(SM), file(mantasv), file(dellysv), file(tidditsv) from variant_db

    output:
        file("${SM}_merged_caller.vcf.gz") into sample_merged_svcaller_vcf
        file("${SM}_merged_caller.vcf.gz.csi") into sample_merged_svcaller_index

    when:
        params.tiddit_discord

    script:
      """
        echo ${SM} > samplename.txt

        svdb --merge --pass_only --no_var --vcf ${dellysv}:delly ${mantasv}:manta ${tidditsv}:tiddit --priority delly,manta,tiddit| \\
        bcftools reheader -s samplename.txt | \\
        awk '\$0 ~ "#" {print} !seen[\$1"\t"\$2]++ {print}' | \\
        bcftools view -Oz -o ${SM}_merged_caller.vcf.gz

        bcftools index -f ${SM}_merged_caller.vcf.gz
      """
}

sample_merged_svcaller_vcf
    .toSortedList()
    .set { merged_sample_sv_vcf }

sample_merged_svcaller_index
    .toSortedList()
    .set { merged_sample_sv_index }


process merge_wi_sv_callers {

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file(mergedSVvcf) from merged_sample_sv_vcf
        file(mergedSVindex) from merged_sample_sv_index
        file("gene.pkl") from gene_pkl_svdb

    output:
        file "WI.${date}.MERGEDsv.snpEff.vcf.gz" into wi_mergedsv

    when:
        params.tiddit_discord

    """
        bcftools merge -m all --threads ${task.cpus} -O z ${mergedSVvcf} > temp_merged.vcf.gz

        bcftools view -Ov temp_merged.vcf.gz | grep '^#'  > db_merged_temp_sorted.vcf
        bcftools view -Ov temp_merged.vcf.gz | grep -v -E '^X|^MtDNA|^#' | sort -k1,1d -k2,2n >> db_merged_temp_sorted.vcf
        bcftools view -Ov temp_merged.vcf.gz | grep -E '^X' | sort -k1,1d -k2,2n >> db_merged_temp_sorted.vcf
        bcftools view -Ov temp_merged.vcf.gz | grep -E '^MtDNA' | sort -k1,1d -k2,2n >> db_merged_temp_sorted.vcf

        bcftools query -l db_merged_temp_sorted.vcf | sort > sample_names.txt

        bcftools view --samples-file=sample_names.txt -Ov db_merged_temp_sorted.vcf | \\
        snpEff eff -csvStats MERGEDsv_snpeff_out.csv \\
        -no-downstream -no-intergenic -no-upstream \\
        -dataDir ${params.snpeff_path} \\
        -config ${params.snpeff_path}/snpEff.config \\
        ${params.annotation_reference} | \\
        bcftools view -O v | \\
        python `which fix_snpeff_names.py` - | \\
        bcftools view -O z > WI.${date}.MERGEDsv.snpEff.vcf.gz

        bcftools index -f WI.${date}.MERGEDsv.snpEff.vcf.gz
    """
}

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
        file("DELLYsv_snpeff_out.csv") from snpeff_delly_multiqc
        file("MANTAsv_snpeff_out.csv") from manta_snpeff_multiqc
        file("WI.${date}.DELLYsv.raw.stats.txt") from delly_bcf_stats
        file("WI.${date}.MANTAsv.soft-filter.stats.txt") from bcf_manta_stats
        file("WI.${date}.TIDDITsv.soft-filter.stats.txt") from tiddit_stats
        file("TIDDITsv_snpeff_out.csv") from snpeff_tiddit_multiqc

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