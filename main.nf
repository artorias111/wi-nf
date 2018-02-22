#!/usr/bin/env nextflow
/* 
 * Authors: 
 * - Daniel Cook <danielecook@gmail.com>
 *  
 */

/*
    Filtering configuration
*/

min_depth=10
qual=30
mq=40
dv_dp=0.5

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
params.cores = 4
params.tmpdir = "tmp/"
params.email = ""
File reference = new File("${params.reference}")
if (params.reference != "(required)") {
   reference_handle = reference.getAbsolutePath();
} else {
   reference_handle = "(required)"
}

// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.fqs = "${workflow.projectDir}/test_data/SM_sample_sheet.tsv"
    params.bamdir = "${params.out}/bam"
    File fq_file = new File(params.fqs);
    params.fq_file_prefix = "${workflow.projectDir}/test_data"

    // lower filter thresholds
    min_depth=0
    qual=10
    mq=10
    dv_dp=0.0

} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fq_file_prefix = null;
    params.fqs = "SM_sample_sheet.tsv"
}
File fq_file = new File(params.fqs);




/*
    UX
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
    --reference             Reference Genome               ${params.reference}
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
    Kmer counting
*/
process kmer_counting {

    container null


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
    zcat ${fq1} ${fq2} | fastq-kmers -k 6 | awk -v OFS="\t" -v ID=${ID} -v SM=${SM} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${ID}.kmer.tsv
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
    Fastq alignment
*/
// The output looks strange below,
// but its designed to group like samples together - so leave it!

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
  Merge - Generate SM Bam
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
                  bam_SM_stats;
}

process bam_SM_stats {

    cpus params.cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_SM_stats

    output:
         file("${SM}.samtools.txt") into SM_samtools_stats_set
         file("${SM}.bamtools.txt") into SM_bamtools_stats_set
         file("${SM}_fastqc.zip") into SM_fastqc_stats_set
         file("${SM}.picard.*") into SM_picard_stats_set

    """
        samtools stats ${SM}.bam > ${SM}.samtools.txt
        bamtools stats -in ${SM}.bam > ${SM}.bamtools.txt
        fastqc --threads ${task.cpus} ${SM}.bam
        picard CollectAlignmentSummaryMetrics R=${reference_handle} I=${SM}.bam O=${SM}.picard.alignment_metrics.txt
        picard CollectInsertSizeMetrics I=${SM}.bam O=${SM}.picard.insert_metrics.txt H=${SM}.picard.insert_histogram.txt
    """

}

process bam_publish {

    publishDir params.bamdir + "/WI/isotype", mode: 'copy', pattern: '*.bam*'

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_publish
    output:
        set file("${SM}.bam"), file("${SM}.bam.bai")

    """
        echo "${SM} saved to publish folder"
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
    SM bam stats
*/

process isotype_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
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
    Coverage Bam
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
    telseq
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
    Call Variants
*/

process call_variants_individual {

    tag { SM }

    cpus params.cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_snp_individual

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${task.cpus} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference_handle} ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.individual.vcf.gz
    bcftools index ${SM}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${SM}.individual.vcf.gz | \\
    bcftools filter --include 'DP > 3' | \\
    egrep '(^#|1/1)' | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' > ${SM}.individual.sites.tsv

    """
}

process merge_variant_list {

    publishDir params.out + "/variation", mode: 'copy'
    
    input:
        val sites from individual_sites.toSortedList()

    output:
        file("sitelist.tsv.gz") into gz_sitelist
        file("sitelist.tsv.gz") into sitelist_stat
        file("sitelist.tsv.gz.tbi") into gz_sitelist_index


    """
        echo ${sites}
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > sitelist.tsv
        bgzip sitelist.tsv -c > sitelist.tsv.gz && tabix -s1 -b2 -e2 sitelist.tsv.gz
    """
}

union_vcf_set = bam_snp_union.combine(gz_sitelist).combine(gz_sitelist_index)

process call_variants_union {

    tag { SM }

    cpus params.cores

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_set

    output:
        val SM into union_vcf_SM
        file("${SM}.union.vcf.gz") into union_vcf_list

    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${task.cpus} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference_handle} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Concatenate and filter
        bcftools concat \${order} -O v | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --threads ${task.cpus} --mode + --soft-filter quality --include "QUAL >= ${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads ${task.cpus} --mode + --soft-filter min_depth --include "FORMAT/DP > ${min_depth}" | \\
        bcftools filter -O u --threads ${task.cpus} --mode + --soft-filter mapping_quality --include "INFO/MQ > ${mq}" | \\
        bcftools filter -O v --threads ${task.cpus} --mode + --soft-filter dv_dp --include "(FORMAT/AD[1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
        awk -v OFS="\t" '\$0 ~ "^#" { print } \$0 ~ ":AB" { gsub("PASS","", \$7); if (\$7 == "") { \$7 = "het"; } else { \$7 = \$7 ";het"; } } \$0 !~ "^#" { print }' | \\
        awk -v OFS="\t" '\$0 ~ "^#CHROM" { print "##FILTER=<ID=het,Description=\\"heterozygous_call_after_het_polarization\\">"; print; } \$0 ~ "^#" && \$0 !~ "^#CHROM" { print } \$0 !~ "^#" { print }' | \\
        vk geno transfer-filter - | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """
}

process generate_union_vcf_list {

    executor 'local'

    cpus 1 

    input:
       val vcf_set from union_vcf_list.toSortedList()

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
        val(chrom) into contigs_list_in
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        bcftools merge --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}

// Generate a list of ordered files.
contig_raw_vcf = CONTIG_LIST*.concat(".merged.raw.vcf.gz")

process concatenate_union_vcf {

    echo true

    cpus params.cores

    input:
        val merge_vcf from raw_vcf.toSortedList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat --threads ${task.cpus} -O z ${contig_raw_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
    """
}

// Generates the initial soft-vcf; but it still
// needs to be annotated with snpeff and annovar.
process generate_soft_vcf {

    cpus params.cores

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf
        set val('filtered'), file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf_stat
        file("WI.${date}.soft-filter.stats.txt") into soft_filter_stats

    """
        bcftools view merged.raw.vcf.gz | \\
        vk filter MISSING --max=0.90 --soft-filter="high_missing" --mode=x - | \
        vk filter HET --max=0.10 --soft-filter="high_heterozygosity" --mode=+ - | \
        vk filter REF --min=1 - | \
        vk filter ALT --min=1 - | \
        vcffixup - | \\
        bcftools view -O z - > WI.${date}.soft-filter.vcf.gz
        bcftools index -f WI.${date}.soft-filter.vcf.gz
        bcftools stats --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt
    """
}

filtered_vcf.into { 
                    filtered_vcf_snpeff;
                    filtered_vcf_to_clean;
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

process annotate_vcf_snpeff {

    cpus params.cores

    errorStrategy 'retry'
    maxRetries 2

    input:
        set file("merged.WI.${date}.soft-filter.vcf.gz"), file("merged.WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_snpeff
        file("gene.pkl") from gene_pkl

    output:
        set file("WI.${date}.snpeff.vcf.gz"), file("WI.${date}.snpeff.vcf.gz.csi") into snpeff_vcf
        file("snpeff_out.csv") into snpeff_multiqc

    script:
        using_container = !(task.container == null)

        """
            # First run generates the list of gene identifiers
            # If running a docker container, the snpeff database must be built.
            if [ "${using_container}" == "true" ]; then
                setup_annotation_db.sh ${params.annotation_reference}
            fi;

            bcftools view -O v merged.WI.${date}.soft-filter.vcf.gz | \\
            snpEff eff -csvStats snpeff_out.csv -noInteraction -no-downstream -no-intergenic -no-upstream ${params.annotation_reference} | \\
            bcftools view -O v | \\
            python `which fix_snpeff_names.py` - | \\
            bcftools view -O z > WI.${date}.snpeff.vcf.gz
            bcftools index WI.${date}.snpeff.vcf.gz
        """

}


process generate_hard_vcf {

    cpus params.cores

    publishDir params.out + "/variation", mode: 'copy'

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_to_clean

    output:
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf_to_impute
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into tajima_bed
        set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into vcf_phylo
        set val('clean'), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.csi") into hard_vcf
        file("WI.${date}.hard-filter.vcf.gz.tbi")
        file("WI.${date}.hard-filter.stats.txt") into hard_filter_stats


    """
        # Generate hard-filtered (clean) vcf
        bcftools view --types snps WI.${date}.soft-filter.vcf.gz | \\

        bcftools filter --set-GTs . --exclude 'FORMAT/FT != "PASS"' | \\
        vk filter MISSING --max=0.90 - | \\
        vk filter HET --max=0.10 - | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        bcftools filter --threads ${task.cpus} --set-GTs . --exclude "GT != '0/0' && GT != '1/1'" | \\
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
        set file("merged.filtered.snp.vcf.gz"), file("merged.filtered.snp.vcf.gz.csi") from filtered_vcf_gtcheck

    output:
        file("gtcheck.tsv") into gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 merged.filtered.snp.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """

}


process generate_primers {

    input:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf

    output:
        file('primers.tsv')

    """
        vk primer snip --ref=WS245 WI.${date}.soft-filter.vcf.gz > primers.tsv
    """


}

/*
    Calculate Singletons
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
    Phylo analysis
*/
process phylo_analysis {

    publishDir params.out + "/popgen/trees", mode: "copy"

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

    publishDir params.out + "/popgen/trees", mode: "copy"

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

    publishDir params.out + "/popgen", mode: 'copy'

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
        java -jar `which beagle.jar` nthreads=${task.cpus} window=8000 overlap=3000 impute=true ne=17500 gt=WI.${date}.hard-filter.vcf.gz out=WI.${date}.impute
        bcftools index WI.${date}.impute.vcf.gz
        tabix WI.${date}.impute.vcf.gz
        bcftools stats --verbose WI.${date}.impute.vcf.gz > WI.${date}.impute.stats.txt
    """
}


impute_vcf.into { kinship_vcf;  mapping_vcf }

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

process make_mapping {

    publishDir params.out + "/cegwas", mode: 'copy'

    input:
        set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") from mapping_vcf
    output:
        file("snps.Rda")

    """
        Rscript -e 'library(cegwas); snps <- generate_mapping("WI.${date}.impute.vcf.gz"); save(snps, file = "snps.Rda");'
    """

}


/*
    Perform concordance analysis

process process_concordance_results {


    publishDir params.out + "/concordance", mode: "copy"

    input:
        file "gtcheck.tsv" from gtcheck
        file "filtered.stats.txt" from filtered_stats
        file "isotype_coverage.tsv" from isotype_coverage_merged

    output:
        file("concordance.pdf")
        file("concordance.png")
        file("concordance_99.pdf")
        file("concordance_99.png")
        file("isotype_groups.tsv")
        file("gtcheck.tsv")

    """
    Rscript --vanilla `which process_concordance.R`
    """

}
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


process multiqc_report {

    executor 'local'

    publishDir params.out + "/report", mode: 'copy'

    input:
        file("stat*") from vcf_stats.toSortedList()
        file(samtools_stats) from SM_samtools_stats_set.toSortedList()
        //file(bamtools_stats) from SM_bamtools_stats_set.toSortedList()
        file(duplicates) from duplicates_set.toSortedList()
        file(fastqc) from SM_fastqc_stats_set.toSortedList()
        file("bam*.idxstats") from bam_idxstats_multiqc.toSortedList()
        file("picard*.stats.txt") from SM_picard_stats_set.collect()
        file("snpeff_out.csv") from snpeff_multiqc

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
        file("sitelist.tsv.gz") from sitelist_stat

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

    // mail summary
    ['mail', '-s', 'wi-nf', params.email].execute() << summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }



}
