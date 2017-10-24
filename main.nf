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

date = new Date().format( 'yyyy-MM-dd' )
params.out = "WI-${date}"
params.debug = false
params.annotation_reference = "WS261"
params.cores = 4
params.tmpdir = "tmp/"
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

    // lower filter thresholds
    min_depth=0
    qual=10
    mq=10
    dv_dp=0.0

} else {
    // The SM sheet that is used is located in the root of the git repo
    params.fqs = "${workflow.projectDir}/SM_sample_sheet.tsv"
    params.bamdir = "(required)"
}

File fq_file = new File("${params.fqs}");
params.fq_file_prefix = fq_file.getParentFile().getAbsolutePath();


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
fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
             .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${params.fq_file_prefix}/${fq1}"), file("${params.fq_file_prefix}/${fq2}"), seq_folder] }



/*
    Fastq alignment
*/
process perform_alignment {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2, seq_folder from fqs
    output:
        set val(SM), file("${ID}.bam"), file("${ID}.bam.bai") into fq_bam_set

    
    """
        bwa mem -t ${params.cores} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference_handle} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${params.cores} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${params.cores} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
        sambamba index --nthreads=${params.cores} ${ID}.bam

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
        file("${SM}.duplicates.txt") into duplicates_file
        
    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${params.cores} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${params.cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${params.cores} ${SM}.bam
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
                  bam_to_cram
}

process bam_publish {

    publishDir params.bamdir + "/WI/isotype", mode: 'copy', pattern: '*.bam*'

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_publish
    output:
        set file("${SM}.bam"), file("${SM}.bam.bai")

    """
        echo "${SM} saved to publish folder"
    """
}

process convert_to_cram {

    cpus params.cores

    tag { SM }

    publishDir params.bamdir + "/WI/isotype_cram", mode: 'copy', pattern: '*.cram*'

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_to_cram
    output:
        set file("${SM}.cram"), file("${SM}.cram.crai")

    """
        sambamba view --nthreads=${params.cores} --format=cram --ref-filename=${reference_handle} --output-filename ${SM}.cram ${SM}.bam 
        sambamba index --cram-input ${SM}.cram
    """
}

process SM_idx_stats {
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_idxstats
    output:
        file bam_idxstats into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats
    """
}

process SM_combine_idx_stats {

    publishDir params.out + "/isotype", mode: 'copy'

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
    """
}

/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    publishDir params.out + "/SM", mode: 'copy'

    input:
        val stat_files from SM_bam_stat_files.toSortedList()

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}

process format_duplicates {

    publishDir params.out + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv")


    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
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
        val SM into SM_coverage_sample
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process coverage_SM_merge {

    publishDir params.out + "/isotype", mode: 'copy'

    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv") into SM_coverage_merged

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
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

    publishDir params.out + "/SM", mode: 'copy'

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
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${params.cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference_handle} ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
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

    publishDir params.out + "/sitelist", mode: 'copy'
    
    input:
        val sites from individual_sites.toSortedList()

    output:
        file("sitelist.tsv.gz") into gz_sitelist
        file("sitelist.tsv.gz.tbi") into gz_sitelist_index
        file("sitelist.tsv") into sitelist
        file("sitelist.count.txt")


    """
        echo ${sites}
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > sitelist.tsv
        cat sitelist.tsv | wc -l > sitelist.count.txt
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
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${params.cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference_handle} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Concatenate and filter
        bcftools concat \${order} -O v | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --threads ${params.cores} --mode + --soft-filter quality --include "QUAL >= ${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads ${params.cores} --mode + --soft-filter min_depth --include "FORMAT/DP > ${min_depth}" | \\
        bcftools filter -O u --threads ${params.cores} --mode + --soft-filter mapping_quality --include "INFO/MQ > ${mq}" | \\
        bcftools filter -O v --threads ${params.cores} --mode + --soft-filter dv_dp --include "(FORMAT/AD[1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
        awk -v OFS="\t" '\$0 ~ "^#" { print } \$0 ~ ":AB" { gsub("PASS","", \$7); if (\$7 == "") { \$7 = "het"; } else { \$7 = \$7 ";het"; } } \$0 !~ "^#" { print }' | \\
        awk -v OFS="\t" '\$0 ~ "^#CHROM" { print "##FILTER=<ID=het,Description=\\"heterozygous_call_after_het_polarization\\">"; print; } \$0 ~ "^#" && \$0 !~ "^#CHROM" { print } \$0 !~ "^#" { print }' | \\
        vk geno transfer-filter - | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """
}

process generate_union_vcf_list {

    cpus 1 

    publishDir params.out + "/vcf", mode: 'copy'

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

    //publishDir params.out + "/vcf", mode: 'copy'

    input:
        //val chrom from contigs_list_in
        val merge_vcf from raw_vcf.toSortedList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat -O z ${contig_raw_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
    """
}

process generate_soft_vcf {

    publishDir params.out + "/vcf", mode: 'copy'

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf
        set val('filtered'), file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi") into filtered_vcf_stat

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

filtered_vcf.into { filtered_vcf_snpeff; filtered_vcf_to_clean; filtered_vcf_gtcheck }

fix_snpeff_script = file("fix_snpeff_names.py")

process annotate_vcf_snpeff {

    publishDir params.out + "/vcf", mode: 'copy'

    input:
        set file("merged.WI.${date}.soft-filter.vcf.gz"), file("merged.WI.${date}.soft-filter.vcf.gz.csi") from filtered_vcf_snpeff

    output:
        set file("WI.${date}.snpeff.vcf.gz"), file("WI.${date}.snpeff.vcf.gz.csi") into snpeff_vcf

    """
        # First run generates the list of gene identifiers
        setup_annotation_db.sh ${params.annotation_reference}
        fix_snpeff_names.py

        bcftools view -O v merged.WI.${date}.soft-filter.vcf.gz | \\
        snpEff eff -noInteraction -no-downstream -no-intergenic -no-upstream ${params.annotation_reference} | \\
        bcftools view -O v | \\
        python `which fix_snpeff_names.py` - | \\
        bcftools view -O z > WI.${date}.snpeff.vcf.gz
        bcftools index WI.${date}.snpeff.vcf.gz
        bcftools stats --verbose WI.${date}.snpeff.vcf.gz > WI.${date}.snpeff.stats.txt
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

    """
        # Generate hard-filtered (clean) vcf
        bcftools view --types snps WI.${date}.soft-filter.vcf.gz | \\
        bcftools filter --set-GTs . --exclude 'FORMAT/FT != "PASS"' | \\
        vk filter MISSING --max=0.90 - | \\
        vk filter HET --max=0.10 - | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        bcftools filter --threads ${params.cores} --set-GTs . --exclude "GT != '0/0' && GT != '1/1'" | \\
        vcffixup - | \\
        bcftools view -O z > WI.${date}.hard-filter.vcf.gz
        bcftools index -f WI.${date}.hard-filter.vcf.gz
        bcftools stats --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt
    """
}

hard_vcf.set { hard_vcf_summary }


/*
    Calculate Singletons
*/

process calculate_hard_vcf_summary {

    publishDir params.out + "/summary", mode: 'copy'

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

    publishDir params.out + "/phylo", mode: "copy"

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

    publishDir params.out + "/phylo", mode: "copy"

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
        tajima --no-header 100000 10000 WI.${date}.hard-filter.vcf.gz | bgzip > WI.${date}.tajima.bed.gz
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
        set val('impute'), file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.csi") into impute_vcf_stat

    """
        java -jar `which beagle.jar` nthreads=${params.cores} window=8000 overlap=3000 impute=true ne=17500 gt=WI.${date}.hard-filter.vcf.gz out=WI.${date}.impute
        bcftools index -f WI.${date}.impute.vcf.gz
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

/* 
    Stat VCFs
*/

vcf_stat_set = filtered_vcf_stat.concat( impute_vcf_stat)

/*
    Stat IMPUTE and CLEAN Here; Final VCF statd below
*/

process stat_tsv {

    tag { vcf_stat }

    publishDir params.out + "/vcf", mode: 'copy'

    input:
        set val(vcf_stat), file("${vcf_stat}.vcf.gz"), file("${vcf_stat}.vcf.gz.csi") from vcf_stat_set

    output:
        file("WI.${date}.${vcf_stat}.stats.txt") into filtered_stats

    """
        bcftools stats --verbose ${vcf_stat}.vcf.gz > WI.${date}.${vcf_stat}.stats.txt
    """

}

/*
    Perform concordance analysis
*/

process process_concordance_results {


    publishDir params.out + "/concordance", mode: "copy"

    input:
        file "gtcheck.tsv" from gtcheck
        file "filtered.stats.txt" from filtered_stats
        file "SM_coverage.tsv" from SM_coverage_merged

    output:
        file("concordance.pdf")
        file("concordance.png")
        file("xconcordance.pdf")
        file("xconcordance.png")
        file("isotype_groups.tsv")
        file("gtcheck.tsv")

    """
    Rscript --vanilla `which process_concordance.R`
    """

}

process download_annotation_files {
    executor 'local'

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


process soft_filter_vcf {

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
        vcfanno -p ${params.cores} vcf_anno.conf WI.${date}.snpeff.vcf.gz | \\
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
